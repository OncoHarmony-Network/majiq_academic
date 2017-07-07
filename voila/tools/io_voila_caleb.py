from voila.tools import Tool
from voila.tools import find_voila_files
import os
import pandas as pa
import pdb
import copy


class ThisisFindBinaryLSVs(Tool):
    help = 'Given a list of voila dPSI txt files, return LSVs that are binary-like'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.1)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.0)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="tsv",
                            type=str,
                            help=help_mes)
        help_mes = "Flag: don't consider IR LSVs"
        parser.add_argument('--no_ir',
                            action='store_true',
                            help=help_mes,
                            default=False)

        return parser

    def run(self, args):
        # this is for code readability, not efficiency
        consider_ir = True
        if args.no_ir:
            consider_ir = False
        imported = quick_import(dir=args.directory,
                                cutoff_d_psi=0,
                                cutoff_prob=0,
                                pattern=args.pattern,
                                keep_ir=consider_ir)


def quick_import(dir,
                 cutoff_d_psi=0.2,
                 cutoff_prob=0.95,
                 keep_ir=True,
                 cutoff_sum=False,
                 return_funky_ids=False,
                 pattern="tsv",
                 deseq_dir=False,
                 deseq_prefix=False,
                 deseq_pat="*csv",
                 deseq_log2fc=1.2,
                 deseq_pv=0.05,
                 deseq_sep="\t",
                 deseq_id_colname="ensembl_gene_id",
                 just_one=False,
                 stop_at=False):
    """
    Given a directory with '*_quantify_deltapsi' files, import all dPSI
        text files, throwing away intron events, return dictionary as follows.
        Key : Value
        NAME_OF_COMPARISON : Non-intron containing imported LSV dictionary

    Arguments:
        cutoff_d_psi: LSVs must have at least 1 junction dPSI >= this.
            (and those junctions must meet the Cutoff_prob)
            0<=X<=1
        cutoff_prob: All LSVs must have at least 1 junction prob >= this.
            (and those junctions must meet the Cutoff_dPSI)
            0<=X<=1
        keep_ir: Boolean. Keep LSVs with significantly changing introns?
            If Intron is significantly changing in any of the comparisons,
            the LSV is considered to be an intron retention event.
        cutoff_sum: if True, dPSI cutoff works by:
            sum(dPSI>0)>=Cutoff
            OR
            abs(sum(dPSI<0))>=Cutoff
        return_funky_ids: Return LSV IDs for those LSVs that I think
            are odd and will break some of my other programs?
        pattern: What does the base voila deltapsi output file format look like?\
            Default (pre voila 1.0.0) looks like:
                .deltapsi_quantify_deltapsi.txt
            voila 1.0.0 looks like:
                .deltapsi_deltapsi.tsv
        DESeq*:
            Optional arguments to extract DESeq results, and use them to remove LSVs
                in genes exhibiting significant changes in gene expression.
                See get_deseq_diff_expr_genes() for details.
        just_one: only import and return one txt path (only really need this for lookup.lookup_everywhere()...)
        stop_at: if provided, stop reading voila file when you reach this LSV ID

    Assumptions:
        Directory contains files that end in ".deltapsi_quantify_deltapsi.txt"
            or DIrectory points directly at a deltapsi.txt file.

        Expected format of input files:
            ~/path/NAME_OF_COMPARISON.deltapsi_quantify_deltapsi.txt

    Returns all voila dPSI results in a dictionary format.
    """
    if not (os.path.isdir(dir)):
        if _is_deltapsi_file(dir, Pattern=pattern):
            basename = os.path.basename(dir)
            dpsi_comparison_name = [_get_deltapsi_txt_file_comparison(basename)]
            dpsi_files = [dir]
        else:
            raise ValueError(dir + " not found.")
    else:
        print("Searching for %s files ..." % pattern)
        dpsi_comparison_name, dpsi_files = find_voila_files.get_voila_files(dir,
                                                                            pattern=pattern,
                                                                            get_comp_names=True)
        if len(dpsi_comparison_name) != len(dpsi_files):
            raise ValueError("Something is probably screwy with the names "
                             "of the dPSI text files...")
    if len(dpsi_files) == 0:
        raise RuntimeError("Didn't find any voila txt files...")
    print("Found " + str(len(dpsi_files)) +
          " dPSI text files to import ...")
    imported_files = dict()
    funky_ids = list()

    if just_one:
        dpsi_files = [dpsi_files[0]]
    for f, comparison_name in zip(dpsi_files, dpsi_comparison_name):
        if return_funky_ids:
            imported_file, funk = import_dpsi(f,
                                              cutoff_d_psi,
                                              cutoff_prob,
                                              return_funky_ids=return_funky_ids,
                                              stop_at=stop_at)
            funky_ids.extend(funk)
        else:
            imported_file = import_dpsi(f, cutoff_d_psi, cutoff_prob,
                                        return_funky_ids=return_funky_ids,
                                        stop_at=stop_at)
        imported_files[comparison_name] = imported_file
    if cutoff_d_psi != 0 or cutoff_prob != 0 or not keep_ir:
        imported_files = subset_significant(imported_files,
                                            cutoff_dpsi=cutoff_d_psi,
                                            cutoff_prob=cutoff_prob,
                                            keep_introns=keep_ir,
                                            cutoff_sum=cutoff_sum)
    if return_funky_ids:
        funky_ids = list(set(funky_ids))
        return imported_files, funky_ids
    lsvs_length(imported_files)
    if deseq_dir:
        deseqres = get_deseq_diff_expr_genes(deseq_dir=deseq_dir,
                                             deseq_res_pref_pattern=deseq_prefix,
                                             deseq_fname_pattern=deseq_pat,
                                             log2_foldchange=deseq_log2fc,
                                             pval_adj=deseq_pv,
                                             deseq_delim=deseq_sep,
                                             deseq_geneid_colname=deseq_id_colname,
                                             recursive=False)
        remove_genes(imported_files, deseqres)
    return imported_files


def import_dpsi(fp,
                cutoff_d_psi=0,
                cutoff_prob=0,
                keep_cutoff_info=False,
                return_funky_ids=False,
                stop_at=False):
    """
    Given a file path pointing at voila .txt output,
    return a dictionary with LSV_ID->data about that LSV.

    Funky_ids are LSVs with fewer exon coordinates than there are
        junctions... these break downstream code of mine

    Arguments:
        fp: file path for a voila d_psi quantify text file
        cutoff_d_psi: LSVs must have at least 1 junction d_psi >= this.
            (and those junctions must meet the Cutoff_prob)
            0<=X<=1
        cutoff_prob: All LSVs must have at least 1 junction prob >= this.
            (and those junctions must meet the Cutoff_dPSI)
            0<=X<=1
        keep_cutoff_info: depcreated, don't use
        return_funky_ids: deprecated, don't use
        stop_at: if provided, stop reading voila file when you reach this LSV ID


    """
    if not isinstance(fp, str):
        raise TypeError("Expected file path to be string, instead it was %s" % type(fp))
    lsv_dictionary = dict()
    funky_lsvs = list()
    pre_voila_1_0_0 = True
    has_voila = True
    file_headers = list()
    with open(fp, "r") as handle:
        line_i = 0
        found_stop_at = False
        can_stop = False
        for line in handle:
            if isinstance(stop_at, str):
                if stop_at in line:
                    found_stop_at = True
                    can_stop = True
                else:
                    found_stop_at = False
            line_split = line.rstrip("\r\n").split("\t")
            if line_i == 0:
                # Fix pound sign silliness
                line_split[line_split.index("#Gene Name")] = "Gene Name"
                file_headers.extend(line_split)
                condition_1_name = line_split[5].split(" ")[0]
                condition_2_name = line_split[6].split(" ")[0]
                print("Importing %s vs %s deltapsi data ..." % (condition_1_name, condition_2_name))
                if "Voila Link" in file_headers:
                    has_voila = True
                else:
                    has_voila = False
                # p_thresh = line_split[4]
                # if DPSI_HEADER[4] == "TBD":
                #     left = p_thresh.find("=")
                #     right = p_thresh.find(") ")
                #     p_thresh_float = float(p_thresh[left+1:right])
                #     DPSI_HEADER[4] = p_thresh
                #     print "Voila results were generated with threshold d_psi of %s"%(p_thresh_float)
                # elif DPSI_HEADER[4] != p_thresh:
                #     ERROR="The deltapsi txtfiles were generated with different thresholds."
                #     ERROR+=" That isn't currently supported (and you probably don't want to"
                #     ERROR+=" compare these results ... talk to Caleb."
                #     raise RuntimeError(ERROR)
                line_i += 1
                continue

            gene_name = str(line_split[0])

            gene_id = str(line_split[1])

            d_psi_floated = [float(x) for x in line_split[3].split(';')]

            prob_d_psis_floated = [float(x) for x in line_split[4].split(';')]

            # Saves the junction # and directionality of change
            over_cutoff = list()

            # Just saves which junction index for those that are significant
            sig_junctions = list()

            if keep_cutoff_info:
                # Junction index
                i = 0
                for d_psi, prob_dPSI in zip(d_psi_floated, prob_d_psis_floated):
                    # If the d_psi is > than 0 and > the cutoff d_psi and prob.
                    if ((d_psi > cutoff_d_psi) and (prob_dPSI > cutoff_prob)):
                        # junctionIndex_over
                        over_cutoff.append(str(i) + "_over")
                        sig_junctions.append(i)
                    # Else if the d_psi is < than 0 and > the cutoff d_psi and prob.
                    elif ((d_psi < -cutoff_d_psi) and (prob_dPSI > cutoff_prob)):
                        # junctionIndex_under
                        over_cutoff.append(str(i) + "_under")
                        sig_junctions.append(i)
                    if abs(d_psi) < 0.2 and prob_dPSI >= 0.95:
                        pdb.set_trace()
                    i += 1

            LSV_ID = str(line_split[2])
            if "target" in LSV_ID:
                ref_type = "target"
            elif "source" in LSV_ID:
                ref_type = "source"

            psi_1_floated = [float(x) for x in line_split[5].split(';')]

            psi_2_floated = [float(x) for x in line_split[6].split(';')]

            lsv_type = str(line_split[7])

            a5ss = bool(line_split[8])

            a3ss = bool(line_split[9])

            es = bool(line_split[10])

            n_junctions = int(line_split[11])

            n_exons = int(line_split[12])

            if line_i == 1:
                if "True" in str(line_split[13]):
                    pre_voila_1_0_0 = False
                if "False" in str(line_split[13]):
                    pre_voila_1_0_0 = False
            if pre_voila_1_0_0:
                de_novo_junct = int(line_split[13])
            else:  # Else it is boolean
                de_novo_junct = str(line_split[13])

            chrom = str(line_split[14])

            strand = str(line_split[15])

            junct_coord = str(line_split[16]).split(";")

            exon_coord = str(line_split[17]).split(";")

            if len(d_psi_floated) != len(exon_coord) - 1:
                funky_lsvs.append(LSV_ID)
                # line_i+=1
                # continue

            exon_alt_start = str(line_split[18])

            exon_alt_end = str(line_split[19])

            ir_coords = str(line_split[20])

            if has_voila:
                voila_link = str(line_split[21])

            lsv_dictionary[LSV_ID] = dict({
                file_headers[0]: gene_name,
                file_headers[1]: gene_id,
                file_headers[2]: LSV_ID,
                file_headers[3]: d_psi_floated,
                file_headers[4]: prob_d_psis_floated,
                file_headers[5]: psi_1_floated,
                file_headers[6]: psi_2_floated,
                file_headers[7]: lsv_type,
                file_headers[8]: a5ss,
                file_headers[9]: a3ss,
                file_headers[10]: es,
                file_headers[11]: n_junctions,
                file_headers[12]: n_exons,
                file_headers[13]: de_novo_junct,
                file_headers[14]: chrom,
                file_headers[15]: strand,
                file_headers[16]: junct_coord,
                file_headers[17]: exon_coord,
                file_headers[18]: exon_alt_start,
                file_headers[19]: exon_alt_end,
                file_headers[20]: ir_coords})
            if has_voila:
                lsv_dictionary["Voila Link"] = voila_link

            if keep_cutoff_info:
                lsv_dictionary[LSV_ID]["dPSI_over_cutoff"] = over_cutoff
                lsv_dictionary[LSV_ID]["sig_junctions"] = sig_junctions

            lsv_dictionary[LSV_ID]["Reference_Type"] = ref_type

            line_i += 1
            if can_stop:
                if not found_stop_at:
                    break

        lsv_dictionary["condition_1_name"] = condition_1_name
        lsv_dictionary["condition_2_name"] = condition_2_name

        n_lsvs = str(line_i - 1)
        n_funky = str(len(funky_lsvs))
        fn = str(os.path.basename(fp))
        # print "%s LSVs (%s funky ones discarded) extracted from %s"%(n_lsvs,n_funky,fn)
        if return_funky_ids:
            return lsv_dictionary, funky_lsvs
        return lsv_dictionary


def subset_significant(data,
                       cutoff_dpsi=0.2,
                       cutoff_prob=0.95,
                       keep_introns=False,
                       cutoff_sum=False,
                       intron_dpsi_thresh=0.05):
    """
    Given a quick_import dictionary, copy each LSV dicitonary to a new dictionary
        whereby only LSVs that meet provided cutoff settings are retained.

        See: help(get_sig_lsv_ids) for argument details.
            Briefly: Cutoff_sum, when True, means LSV is considered
            significant if the sum of +dPSIs>=cutoff or abs(sum(-dPSI))>=cutoff.
            Introns_across_comparisons: if Keep_introns=False, and if this arg is True,
            if intron retention changing in any comparison, that LSV is considered
            intronic and is discarded from all
    """
    if check_is_quick_import(data, the_bool=True):
        new_dict = dict()
        lsv_dict_names = data.keys()
        for LSV_dict_name in lsv_dict_names:
            lsv_dict = data[LSV_dict_name]
            new_dict[LSV_dict_name] = subset_significant(data=lsv_dict,
                                                         cutoff_dpsi=cutoff_dpsi,
                                                         cutoff_prob=cutoff_prob,
                                                         keep_introns=keep_introns,
                                                         cutoff_sum=cutoff_sum,
                                                         intron_dpsi_thresh=intron_dpsi_thresh)
        return new_dict
    elif check_is_lsv_dict(data, da_bool=True):
        over_cutoff_ids = get_sig_lsv_ids(data, cutoff_dpsi, cutoff_prob, cutoff_sum)
        if keep_introns:
            ids_to_keep = over_cutoff_ids
        else:
            non_intronic_ids = no_intron_retention(data,
                                                   Return_intronic_ids=False,
                                                   Intron_dPSI_thresh=intron_dpsi_thresh)
            non_intronic_ids = set(non_intronic_ids)
            over_cutoff_ids = set(over_cutoff_ids)
            ids_to_keep = list(over_cutoff_ids & non_intronic_ids)
        over_cutoff_dict = lsv_dict_subset(data, ids_to_keep, save_LSV_data=True)
        return over_cutoff_dict
    else:
        raise ValueError("subset_significant only takes quick_import-style or LSV dictionaries.")


def _get_deltapsi_txt_file_comparison(File):
    if os.path.exists(File):
        basename = os.path.basename(File)
    else:
        basename = File
    split_file_name = basename.split(".")
    comparison_name = split_file_name[0]
    return comparison_name


def _is_deltapsi_file(File, Pattern=".deltapsi_quantify_deltapsi.txt"):
    if not os.path.exists(File):
        raise ValueError("Suposed deltapsi txt file doesn't exist.")
    if Pattern in File:
        return True
    else:
        return False


def _recursive_dirs(Directory):
    """
    Given a directory, return all sub directories, recursively.
    """
    is_dir = list()
    for thing in os.listdir(Directory):
        thing = os.path.join(Directory, thing)
        if os.path.isdir(thing):
            is_dir.append(thing)
            is_dir.extend(_recursive_dirs(thing))
    return is_dir


def no_intron_retention(LSV_dict,
                        Intron_dPSI_thresh=0.05,  # TODO: account for PSI?
                        Return_intronic_ids=False,
                        Consider_other_junctions=False):
    """
    Given LSV dictionary, return LSV IDs corresponding to
        non-significantly changing intron-retention junctions.

    If dPSI of the intron retention junction is < 0.05, the LSV IS
        included in the return dictionary. Using 0.05 instead of 0.2
        because, by definition, at this point the LSV is complex.
        Thus, a dPSI of 0.2 is potentially *very* significant.

    Arguments:
    Return_intronic_ids: boolean
    Consider_other_junctions: boolean
        If true, determine
    """
    # names AKA LSV IDs
    check_is_lsv_dict(LSV_dict)
    names = LSV_dict.keys()
    non_intron_names = list()
    intron_names = list()

    for name in names:
        if name == "condition_1_name" or name == "condition_2_name":
            non_intron_names.append(name)
            intron_names.append(name)
            continue
        # If LSV type is intron:
        if LSV_dict[name]["LSV Type"][-1:] == "i":
            # Save these names, too, just in case I want em later
            intron_names.append(name)
        else:
            non_intron_names.append(name)

    for intron in intron_names:
        if intron == "condition_1_name" or intron == "condition_2_name":
            continue
        dPSI_intron = LSV_dict[intron]["E(dPSI) per LSV junction"][-1]
        # If dPSI of the intron junction is not changing significantly, keep it
        #     0.05 was suggested by Yoseph Barash
        if abs(dPSI_intron) < Intron_dPSI_thresh:
            non_intron_names.append(intron)
    if Return_intronic_ids:
        return non_intron_names, intron_names
    return non_intron_names


def get_deseq_diff_expr_genes(deseq_dir,
                              deseq_res_pref_pattern,
                              deseq_fname_pattern="*csv",
                              log2_foldchange=1.2,
                              pval_adj=0.05,
                              deseq_delim="\t",
                              deseq_geneid_colname="ensembl_gene_id",
                              recursive=True):
    """
        Import DESeq results, then identfy which genes are sig diff expr,
            then return a dictionary with the comparison name -> list of Gene IDs.
            Gene ids are sorted by significance

            Sig genes have:
                log 2 fold change           >= Log2FoldChange
                Benjamani Hochberg P-Val    <= BenjHochPV
                Recursive                   Should the deseq files be recursively searched for?
    """
    if deseq_res_pref_pattern[-1] != "_":
        deseq_res_pref_pattern = deseq_res_pref_pattern + "_"
    deseq_fps = find_voila_files.find_files(Path=deseq_dir, Pattern=deseq_fname_pattern, Recursive=recursive)
    if len(deseq_fps) == 0:
        raise RuntimeError(
            "No DESeq files found with pattern: '%s' in Directory:\n%s" % (deseq_fname_pattern, deseq_dir))
    results = dict()
    for deseq_fp in deseq_fps:
        deseq_pa = pa.read_csv(deseq_fp, sep=deseq_delim, header=0)
        log2fc = deseq_pa["log2FoldChange"].abs() >= log2_foldchange
        benhocpv = deseq_pa["padj"] <= pval_adj
        sig_deseq_pa = deseq_pa[log2fc & benhocpv]
        sig_deseq_pa = sig_deseq_pa.sort(["padj"])
        sig_gene_ids = sig_deseq_pa[deseq_geneid_colname].tolist()
        basename = os.path.basename(deseq_fp)
        compname = basename.split(deseq_res_pref_pattern)[1]
        compname = compname.split('_expression_changes.csv')[0]
        results[compname] = sig_gene_ids
    return results


def get_deseq_genes(DESeqDirectory,
                    Prefix,
                    Pattern="*csv",
                    Sep="\t",
                    Gene_ID_ColName="ensembl_gene_id",
                    Recursive=True):
    """
    Import DESeq results, and sort the results by highest to lowest baseMean.
        Return dict[comparison]:[list of ensembl ids sorted by expression level]
    """
    if Prefix[-1] != "_":
        Prefix = Prefix + "_"
    deseq_fps = find_voila_files.find_files(Path=DESeqDirectory, Pattern=Pattern, Recursive=Recursive)
    if len(deseq_fps) == 0:
        raise RuntimeError("No DESeq files found with pattern: '%s' in Directory:\n%s" % (Pattern, DESeqDirectory))
    results = dict()
    for deseq_fp in deseq_fps:
        deseq_pa = pa.read_csv(deseq_fp, sep=Sep, header=0)
        deseq_pa = deseq_pa.sort(["baseMean"], ascending=False)
        gene_ids = deseq_pa[Gene_ID_ColName].tolist()
        basename = os.path.basename(deseq_fp)
        compname = basename.split(Prefix)[1]
        compname = compname.split('_expression_changes.csv')[0]
        results[compname] = gene_ids
    return results


def check_is_quick_import(quick_imp, the_bool=False):
    """
    Check if Quick_import looks like something returned by quick import.
    """
    if not isinstance(quick_imp, dict):
        if not the_bool:
            raise TypeError("Expected a dictionary.")
        else:
            return False
    keys = quick_imp.keys()
    if len(keys) == 0:
        if the_bool:
            return False
        else:
            raise ValueError("This dictionary is empty...")
    try:
        for key in keys:
            check_is_lsv_dict(quick_imp[key])
    except BaseException as err_mess:
        if not the_bool:
            raise ValueError(str(err_mess))
        else:
            return False
    return True


def check_is_lsv_dict(lsv_dict, da_bool=False):
    """
    Check if LSV_dict looks like what is returned by quick_import()
    """
    if not isinstance(lsv_dict, dict):
        if not da_bool:
            raise TypeError("Expected a dictionary.")
        else:
            return da_bool
    lsv_ids = list()
    lsv_ids.extend(lsv_dict.keys())
    try:
        lsv_ids.remove("condition_1_name")
        lsv_ids.remove("condition_2_name")
    except:
        if not da_bool:
            raise ValueError("Expected a LSV dictionary.")
        else:
            return False
    if len(lsv_ids) == 0:
        if da_bool:
            return True  # this happens if the LSV_dict is empty
        else:
            return
    if not isinstance(lsv_dict[lsv_ids[0]], dict):
        if not da_bool:
            raise TypeError("Expected a dictionary of dictionaries.")
        else:
            return False
    if check_is_lsv(lsv_dict, True):
        if not da_bool:
            raise ValueError("Expected a dictionary of LSVs,"
                             " not an LSV itself.")
        else:
            return False

    return True


def check_is_lsv(lsv, Bool=False):
    """
    Check that LSV is, in fact, a LSV.
    """
    if not isinstance(lsv, dict):
        if not Bool:
            raise ValueError("Expected an LSV, instead object is of type: "
                             + str(type(lsv)))
        else:
            return False
    for header in ["Gene Name",
                   "Gene ID",
                   "LSV ID",
                   "E(dPSI) per LSV junction"]:
        if header not in lsv:
            if not Bool:
                raise ValueError("Expected a LSV, but didn't get one.")
            else:
                return False
    keys = list(lsv.keys())
    sub_dict = lsv[keys[0]]
    if isinstance(sub_dict, dict):
        keys = sub_dict.keys()
        if "condition_1_name" in keys:
            if not Bool:
                raise ValueError("Please pick one comparison from the dictionary, "
                                 "not the entire dictionary returned by quick_import()."
                                 " For entire dictionary, use lookup_everywhere() "
                                 "instead of lookup().")
            else:
                return False
    return True


def lsv_dict_subset(dictionary,
                    keys,
                    save_LSV_data=True,
                    new_sub_key=None,
                    new_values=None,
                    remove_empties=True):
    """
    Return subset of Dictionary using Keys.

    Arguments:
        dictionary : LSV dict...
        keys : keys to extract from LSV dict (list or set)
        save_LSV_data: if True, save "dPSI_over_cutoff",
                                        "sig_junctions",
                                        "condition_1_name",
                                        and "condition_2_name"
                                        keys from the dictionary

        new_sub_key/new_values:
            if new_values and new_sub_key, add the new_values to the dictionary. Must be a list
                of the same length as Keys. Assumption: the Dictionary values are
                themselves dictionaries. Thus, Dictionary[key][new_key] = new_value

        remove_empties: if True, remove all key:value pairs where the value is None
    """
    if not isinstance(dictionary, dict):
        raise ValueError("Dictionary needs to be a dict")
    if not isinstance(keys, list):
        if isinstance(keys, set):
            keys = list(keys)
        elif isinstance(keys, str):
            keys = [keys]
        else:
            raise ValueError("Keys needs to be a string, list, or set.")
    if new_values:
        if not new_sub_key:
            raise ValueError("Must provide both a new key and new values")
        if not isinstance(new_values, list):
            raise ValueError("New_values, which is optional, must be a list of provided.")
        if len(new_values) != len(keys):
            raise ValueError(
                "New_values (" + str(len(new_values)) + ") wasn't the same length as provided keys (" + str(
                    len(keys)) + ")")
    if new_sub_key:
        if not new_values:
            raise ValueError("Must provide both new_sub_key AND new_values.")
    # Copy subset of dictionary using found names.
    #     The None is not nec., because I know all keys will be in this
    #     Dict, but for future Caleb/aliens modifying this code, I'm keeping it.
    new_dict = {k: dictionary.get(k, None) for k in keys}
    if remove_empties:
        n_rem = 0
        for key, value in new_dict.items():
            if value == None:
                n_rem += 1
                del new_dict[key]
                # print n_rem, "'None' values removed from LSV dictionary..."

    # If it's an LSV dictionary, and if you want to Save the condition names:
    if save_LSV_data:
        if "condition_1_name" in dictionary:
            # new_dict["dPSI_over_cutoff"]=Dictionary["dPSI_over_cutoff"]
            # new_dict["sig_junctions"]=Dictionary["sig_junctions"]
            new_dict["condition_1_name"] = dictionary["condition_1_name"]
            new_dict["condition_2_name"] = dictionary["condition_2_name"]

    if None in new_dict.values():
        print("Warning: at least 1 key wasn't found in the dictionary."
              "value of None was assigned to such keys.")
    if new_values:
        for key, new_value in zip(keys, new_values):
            new_dict[key][new_sub_key] = new_value
    return new_dict


def lsvs_length(data):
    """
    Given quick import, print the number of LSVs in each comparisons
    """
    check_is_quick_import(data)
    all_lsvs = get_all_unique_lsv_ids(data, verbose=True)
    if len(data.keys()) > 1:
        all_lsvs = list(set(all_lsvs))
        n_all = len(all_lsvs)
        shared_lsvs = get_shared_slv_ids(data, bool=True)
        if shared_lsvs:
            n_shared = len(shared_lsvs)
        else:
            n_shared = 0
        print("%s Shared LSVs between all comparisons." % (n_shared))
        print("Total of %s unique LSVs across all comparisons." % (n_all))


def remove_genes(rm_data,
                 comp_to_lsvid_dict,
                 return_diff_expr_lsvids=False):
    """
        Given a quick import of rm_data and a comp_to_lsvid_dict,
            remove all the LSVs corresponding to each ID.

            return_diff_expr_lsvids : return the LSV IDs that were diff expressed?

            NOTE: THIS MODFIES THE QUICK IMPORT OBJECT.
    """
    check_is_quick_import(rm_data)
    comps_with_deseq_res = comp_to_lsvid_dict.keys()
    for comp in comps_with_deseq_res:
        if comp not in rm_data:
            raise RuntimeError("%s was found in DESeq directory, but not in Majiq results...")
        sig_gene_ids = comp_to_lsvid_dict[comp]
        all_lsv_ids = get_LSV_IDs(rm_data[comp])
        n_s = len(sig_gene_ids)
        n_t = len(all_lsv_ids)
        n_r = 0
        removed_ids = []
        for lsv_id in all_lsv_ids:
            gene = lsv_id.split(":")[0]
            if gene in sig_gene_ids:
                n_r += 1
                removed_ids.append(rm_data[comp].pop(lsv_id))
        print("%s: %s LSVs (%s diff-expr genes) removed leaving %s LSVs" % (comp, n_r, n_s, n_t - n_r))
    data_comps = set(rm_data.keys())
    comps_with_deseq_res = set(comps_with_deseq_res)
    leftover = list(comps_with_deseq_res - data_comps)
    if len(leftover) > 0:
        print("The following comparisons didn't have DESeq results:\n%s" % leftover)
    if return_diff_expr_lsvids:
        return removed_ids


def get_sig_lsv_ids(Data,
                    Cutoff_dPSI=0,
                    Probability_dPSI=0,
                    Sum_for_cutoff=False):
    """
    Given LSV dictionary, return set of unique LSV IDs over cutoff
        Note: recursively handles Quick_import and returns dict of
        comparison_name -> sig_ids set

    Arguments:
        Data: output of import_dpsi
        Cutoff_dPSI: Only return LSV IDs with at least 1 junction abs(dPSI) >= Cutoff_dPSI
        Probability_dPSI: junction must have a dPSI>= this to be considered
        Sum_for_cutoff: If the sum of all +dPSI (that meet probability cutoff)
            is >= Cutoff_dPSI, then include the LSV. Or if abs(sum of dPSI<0) is >= Cutoff_dPSI.
            This is less conservative than default.

    Return:
        Set
    """
    if check_is_quick_import(Data, the_bool=True):
        comparisons = Data.keys()
        results = dict()
        for comparison in comparisons:
            LSV_Dict = Data[comparison]
            results[comparison] = get_sig_lsv_ids(LSV_Dict,
                                                  Cutoff_dPSI,
                                                  Probability_dPSI,
                                                  Sum_for_cutoff)
        return results
    check_is_lsv_dict(Data)

    # names AKA LSV IDs
    names = get_LSV_IDs(Data)
    names_over_cutoff = set()
    if len(names) < 1:
        raise RuntimeError("No LSVs made Cutoff dPSI of %s and Prob of %s" % (Cutoff_dPSI,
                                                                              Probability_dPSI))
    prob_name = get_name_of_prob_key(Data[names[0]])
    for name in names:
        dPSIs = Data[name]["E(dPSI) per LSV junction"]
        prob_dPSIs = Data[name][prob_name]
        if Sum_for_cutoff:
            sum_dPSI_over = 0
            sum_dPSI_under = 0
            meets_prob_cutoff = False
            for dPSI, prob_dPSI in zip(dPSIs, prob_dPSIs):
                if prob_dPSI < Probability_dPSI:
                    continue
                else:
                    meets_prob_cutoff = True
                if dPSI > 0:
                    sum_dPSI_over += dPSI
                elif dPSI < 0:
                    sum_dPSI_under -= dPSI
            if sum_dPSI_over >= Cutoff_dPSI or abs(sum_dPSI_under) >= Cutoff_dPSI:
                if meets_prob_cutoff:
                    names_over_cutoff.add(name)
        else:
            for dPSI, prob_dPSI in zip(dPSIs, prob_dPSIs):
                if (abs(dPSI) >= Cutoff_dPSI) and (prob_dPSI >= Probability_dPSI):
                    names_over_cutoff.add(name)
    return names_over_cutoff


def get_LSV_IDs(LSV_dict):
    """
    Return LSV IDs from dictionary.
    """
    check_is_lsv_dict(LSV_dict)
    lsv_ids = copy.copy(list(LSV_dict.keys()))
    lsv_ids.remove("condition_1_name")
    lsv_ids.remove("condition_2_name")
    return lsv_ids


def get_dpsis(lsv, prob_cutoff=None):
    """
    Given LSV, return list of dPSIs over Prob_Cutoff, if provided.
        If not Prob_cutoff provided, return all dPSIs
    """
    check_is_lsv(lsv)
    if prob_cutoff:
        if prob_cutoff < 0 or prob_cutoff > 1:
            raise ValueError("'Prob_Cutoff' needs to be >=0 or <=1, not: "
                             + str(prob_cutoff))
    all_dPSI = copy.copy(lsv["E(dPSI) per LSV junction"])
    if prob_cutoff:
        probs = get_probs(lsv)
        dPSI = list()
        # Transfer max dPSI from all_dPSI to dPSI until
        #  probability of dPSI is below cutoff.
        while len(probs) > 0 and max(probs) >= prob_cutoff:
            index_of_max_prb = probs.index(max(probs))
            probs.pop(index_of_max_prb)
            dPSI.append(all_dPSI.pop(index_of_max_prb))
        return dPSI
    else:
        return all_dPSI


def get_strand(LSV):
    check_is_lsv(LSV)
    return copy.copy(LSV["strand"])


def get_chr(LSV):
    check_is_lsv(LSV)
    return copy.copy(LSV["chr"])


def get_probs(LSV):
    """
    Given LSV, return exons coords
    """
    check_is_lsv(LSV)
    return copy.copy(LSV[get_name_of_prob_key(LSV)])


def get_juncs(LSV):
    """
    Given LSV, return exons coords
    """
    check_is_lsv(LSV)
    return copy.copy(LSV["Junctions coords"])


def get_exons(LSV):
    """
    Given LSV, return exons coords
    """
    check_is_lsv(LSV)
    return copy.copy(LSV["Exons coords"])


def get_all_dpsis(LSVs, Prob_Cutoff=0):
    """
    Given LSV dictionary or quick_import structure, return list of all dPSIs from all
        LSVs over Prob_Cutoff.
    """
    all_dPSI = list()
    try:
        check_is_lsv_dict(LSVs)
        lsvs = get_LSV_IDs(LSVs)
        # Extract dPSI from each LSV, using cutoff.
        for lsv_id in lsvs:
            all_dPSI.extend(get_dpsis(LSVs[lsv_id], Prob_Cutoff))
    except:
        try:
            check_is_quick_import(LSVs)
            for key in LSVs.keys():
                LSV_dict = LSVs[key]
                all_dPSI.extend(get_all_dpsis(LSV_dict))
        except:
            raise ValueError("Expected LSV dictionary or quick_import return value.")
    return all_dPSI


def get_psis(lsv, cond_1=False, cond_2=False, as_dict=False):
    """
    :param lsv:
    :param cond_1: string, optional
    :param cond_2: string, optional
    :param as_dict: if True, return as dict with cond:psis
    :return: [cond1 psis, cond2 psis]
    """
    check_is_lsv(lsv)
    if not cond_1:
        psi_key_1, psi_key_2 = get_name_of_psi_keys(lsv)
    elif cond_1 and not cond_2 or not cond_1 and cond_2:
        raise RuntimeError("Please specify both conds, not just one...")
    else:
        psi_key_1 = cond_1
        psi_key_2 = cond_2
    cond_1_psi = copy.copy(lsv[psi_key_1])
    cond_2_psi = copy.copy(lsv[psi_key_2])
    if as_dict:
        return {psi_key_1.split(" ")[0]: cond_1_psi, psi_key_2.split(" ")[0]: cond_2_psi}
    return [cond_1_psi, cond_2_psi]


def get_all_unique_lsv_ids(data,
                           verbose=False):
    """
    Given a quick import format, return all unique LSV IDs
        seen across all all LSV DIctionaries.
    """
    check_is_quick_import(data)
    all_lsvs = list()
    comparisons = list(data.keys())
    comparisons.sort()
    for comparison_name in comparisons:
        lsv_dict = data[comparison_name]
        lsvs = get_LSV_IDs(lsv_dict)
        all_lsvs.extend(lsvs)
        if verbose:
            n_lsvs = len(lsvs)
            print("%s LSVs in %s" % (n_lsvs, comparison_name))
    return list(set(all_lsvs))


def get_shared_slv_ids(data, bool=False):
    """
    Given a quick import format, return LSV IDs
        that are shared by all LSV DIctionaries.
    """
    check_is_quick_import(data)
    lsvids = list()
    for comparison in data.keys():
        lsvids.append(set(get_LSV_IDs(data[comparison])))
    union_of_lsv_ids = set.intersection(*lsvids)
    if len(union_of_lsv_ids) == 0:
        if bool:
            return False
        raise ValueError("No LSV IDs common to all comparisons... are you" +
                         " sure they come from the same MAJIQ run? Or perhaps" +
                         " you filtered them somehow before using this function?")
    return union_of_lsv_ids


def get_name_of_psi_keys(lsv):
    """
    Given an LSV, return the LSV key name corresponding to
     the E(PSI) for condition 1 and condition 2.
    """
    count = 0
    keys = []
    for name in lsv.keys():
        if "E(PSI)" in name:
            keys.append(name)
            count += 1
    if count != 2:
        raise RuntimeError("Couldn't find 2 E(PSI) in the LSV...")
    return keys[0], keys[1]


def get_lsvs_quickly(data_qu, lsvids_qu, comparison_qu=False):
    """
    A much faster implementation of get_LSV that isn't nearly
    as careful with displaying errors and such...

    LSV_IDs: may be a single LSV ID or a list of LSV_IDs
    Comparison: Optional, if provided, will use only the
        provided comparison name to get LSVs.
    """
    check_is_quick_import(data_qu)
    if isinstance(lsvids_qu, set):
        lsvids_qu = list(lsvids_qu)
    if check_is_lsv_id(lsvids_qu, Bool=True):
        lsvids_qu = [lsvids_qu]
    else:
        check_is_lsv_id(lsvids_qu[0])
    if len(lsvids_qu) != len(list(set(lsvids_qu))):
        raise RuntimeError("Duplicates in your LSV IDs?")
    if isinstance(comparison_qu, str):
        if comparison_qu not in data_qu:
            raise RuntimeError("Provided comparison isn't valid name..")
        comparisons = [comparison_qu]
    else:
        comparisons = data_qu.keys()
    lsvs = list()
    for lsv_id in lsvids_qu:
        found = False
        for comparison in comparisons:
            if lsv_id in data_qu[comparison]:
                lsvs.append(data_qu[comparison][lsv_id])
                found = True
                break
        if not found:
            raise RuntimeError("%s wasn't found in provided Quick import..." % lsv_id)
    return lsvs


def greatest_dpsi(LSV_dict, Change=1.0, Probability=0.95, Verbose=True):
    """
    Given LSV dictionary, return subset of dictionary where LSVs exhibit
    100% change in exclusion/inclusion of an exon. If nothing is 100%,
    find the biggest dPSI in the LSV dictionary.

    If given a quick_import dict of LSV dicts, this function becomes recursive.

    If a LSV dictionary is emptied, it will not be returned at all.

        Arguments:
            Change: dPSI threshold to start searching for LSVs by.
                Starts at 1.0, then recursively keeps looking for LSVs
                by 5% increments.
            Probability: Set to 0.95 if you want high confidence changing LSVs
                set to 0 if you don't care about the probability.
            Verbose: True/False.

    """
    if check_is_quick_import(LSV_dict, the_bool=True):
        all_changes = dict()
        for LSV_dict_name in LSV_dict.keys():
            if Verbose:
                print("Analyzing " + LSV_dict_name + " ...")
            all_changes[LSV_dict_name] = greatest_dpsi(LSV_dict[LSV_dict_name], Change, Probability, Verbose)
        non_empty_lsv_dicts = remove_empty_lsv_dicts(all_changes, Verbose)
        return non_empty_lsv_dicts

    check_is_lsv_dict(LSV_dict)
    if Change > 1 or Change < 0:
        raise ValueError("Change must be >0 and <=1")
    # Don't sum for the cutoff! (False) <- see function desc
    over_cutoff = get_sig_lsv_ids(LSV_dict, Change, Probability, False)

    # If nothing met cutoff, tell user what highest cuttoff is:
    if len(over_cutoff) == 0:
        if Verbose:
            print("No junctions in the LSVs had a dPSI of " + str(Change) +
                  ", now searching for biggest dPSI...")

        # Keep trying to find LSVs starting at Change %,
        #  and going down bit by bit until a LSV is returned.
        for cutoff in calebs_xrange(Change, 0, -0.01):
            if int(cutoff * 100.0) % int(0.05 * 100.0) == 0:
                if Verbose:
                    print("trying", cutoff, "...")
            # Don't sum for the cutoff! (False) <- see function desc
            over_cutoff = get_sig_lsv_ids(LSV_dict, cutoff, Probability, False)
            if len(over_cutoff) > 0:
                if Verbose:
                    print("Max dPSI identifed as: ", cutoff)
                break

    subset_lsv_dict = lsv_dict_subset(LSV_dict, over_cutoff, True)
    return subset_lsv_dict


def calebs_xrange(start, end, by):
    """
    Does an xrange, but can handle floats instead of integers.

    Warning: floats don't behave quite as you'd expect on computers.
     Be careful how you use this function.
    """
    start = float(start)
    end = float(end)
    by = float(by)
    stuff_after_dec = str(1.0).split(".")[1]
    # Figure out which number has the smallest digits, use that
    # to generate the multiplier (to prevent rounding error)
    for number in start, end, by:
        if len(str(number).split(".")[1]) > len(stuff_after_dec):
            stuff_after_dec = str(number).split(".")[1]
    Multyplier = float(pow(10, len(stuff_after_dec)))
    iteratable = list()
    # I add 'int(By*Multyplier)' to the End here because I like things to be explicitly
    #  defnied. This function starts and ends exactly how you tell it, rather than
    #  Python's stupid way. Blah. Sorry, Caleb was cranky when he wrote this.
    for i in range(int(start * Multyplier), int(end * Multyplier) + int(by * Multyplier), int(by * Multyplier)):
        iteratable.append(float(i) / Multyplier)
    return iteratable


def remove_empty_lsv_dicts(data, print_status=True):
    """
    Given a LSV dictionary or quick_import, return non-empty LSV_dictionaries.
    """
    if check_is_quick_import(data, the_bool=True):
        new_dict = dict()
        for LSV_dict_name in data.keys():
            lsv_dict = data[LSV_dict_name]
            remove_empty = remove_empty_lsv_dicts(lsv_dict, print_status)
            if remove_empty == "empty":
                if print_status:
                    print(LSV_dict_name + " was empty.")
            else:
                new_dict[LSV_dict_name] = lsv_dict
        if len(new_dict) == 0:
            print("Warning! All LSV_dicts were empty...")
        return new_dict
    check_is_lsv_dict(data)
    lsv_ids = get_LSV_IDs(data)
    if len(lsv_ids) > 0:
        return data
    else:
        return "empty"


def get_name_of_prob_key(lsv):
    """
    Given an LSV, return the LSV key name corresponding to
     the probability of dPSI entry.
    """
    for name in list(lsv.keys()):
        if "P(|dPSI|>=" in name and ") per LSV junction" in name:
            return name
        elif "P(|E(dPSI)|>=" in name and ") per LSV junction" in name:
            return name

    raise RuntimeError("Couldn't find probability in the LSV...")


def check_is_lsv_id(Thing, Bool=False):
    """
        Bool: should this function reyturn True/False
            or raise an error? If Bool=False; error.
    """
    if not isinstance(Thing, str):
        return throw_fail(Bool, "Not LSV ID: Expected a string ...")
    if ":" not in Thing:
        return throw_fail(Bool, "Not LSV ID: Expected ':' ...")
    split_thing = Thing.split(":")
    if len(split_thing) != 3:
        return throw_fail(Bool, "Not LSV ID: Expected len of 3 ...")
    if split_thing[2] != "source" and split_thing[2] != "target":
        return throw_fail(Bool, "Not LSV ID: Expected source or target")
    # if split_thing[0][0:4] != "ENSG" and split_thing[0][0:7] != "ENSMUSG":
    #     return throw_fail(Bool, "ENSG* or ENSMUSG*")
    if Bool:
        return True
    else:
        pass


def check_lsv_ids_all_shared(Data, Bool=False):
    check_is_quick_import(Data)
    list_of_ids = list()
    for comparison in Data.keys():
        ids = copy.copy(get_LSV_IDs(Data[comparison]))
        list_of_ids.append(set(ids))
    benchmark = list_of_ids[0]
    for i in range(1, len(list_of_ids)):
        diff = benchmark.difference(list_of_ids[i])
        if len(diff) > 0:
            pdb.set_trace()
            if Bool:
                return False
            raise ValueError("LSV IDs are not all shared between these dictionaries.")
    if Bool:
        return True


def get_cond_1_name(lsv_dict):
    """
    :param lsv_dict: lsv_dict struc ...
    :return: str condition 1 name
    """
    check_is_lsv_dict(lsv_dict)
    return lsv_dict["condition_1_name"]


def get_cond_2_name(lsv_dict):
    """
    :param lsv_dict: lsv_dict struc ...
    :return: str condition 2 name
    """
    check_is_lsv_dict(lsv_dict)
    return lsv_dict["condition_2_name"]


def get_comparison_name(LSV_dict, sep="_"):
    """
    Given LSV dictionary, return condition_1_name[sep]condition_2_name
    """
    check_is_lsv_dict(LSV_dict)
    if not isinstance(sep, str):
        raise ValueError("sep needs to be a string, not a: " + str(type(sep)))
    return LSV_dict["condition_1_name"] + sep + LSV_dict["condition_2_name"]


def get_all_LSV_IDs(Data):
    check_is_quick_import(Data)
    comparisons = Data.keys()
    all_ids = list()
    for comparison in comparisons:
        all_ids.extend(get_LSV_IDs(Data[comparison]))
    return list(set(all_ids))


def throw_fail(Bool=False, Message="Test failed."):
    if Bool:
        return False
    else:
        raise RuntimeError(Message)


def copy_lsv(lsv):
    """
    Given an LSV, return a copy that is a
        DIFFERENT OBJECT IN MEMORY !!!

        Thus, you can modify this LSV without messing up th original.
    """
    check_is_lsv(lsv)
    copied = dict()
    for elem in list(lsv.keys()):
        if isinstance(lsv[elem], list):
            copied[elem] = copy.copy(lsv[elem])
        elif isinstance(lsv[elem], str):
            copied[elem] = lsv[elem]
        elif isinstance(lsv[elem], bool):
            copied[elem] = lsv[elem]
        elif isinstance(lsv[elem], int):
            copied[elem] = lsv[elem]
        elif isinstance(lsv[elem], float):
            copied[elem] = lsv[elem]
        else:
            raise RuntimeError("Not sure what to do with %s" % str(type(lsv[elem])))
    return copied


def impute_lsvs(lsvs,
                imputing_with=0):
    """
    Given a list or single LSV, return the LSV(s) with 0s for
        dPSI, prob(dPSI), and PSI

        imputing_with : impute missing values with what?

        voila link is replaced with "LSV WAS BLANKED"
    """
    dpsi_header = ["Gene Name",
                   "Gene ID",
                   "LSV ID",
                   "E(dPSI) per LSV junction",
                   "TBD",
                   "E(PSI)1",  # This header also has the 1st condition name
                   "E(PSI)2",  # This header also has the 2nd condition name
                   "LSV Type",
                   "A5SS",
                   "A3SS",
                   "ES",
                   "Num. Junctions",
                   "Num. Exons",
                   "De Novo Junctions?",
                   "chr",
                   "strand",
                   "Junctions coords",
                   "Exons coords",
                   "Exons Alternative Start",
                   "Exons Alternative End",
                   "IR coords",
                   "Voila link"]
    was_lonely = False
    if check_is_lsv(lsvs, bool):
        was_lonely = True
        lsvs = [lsvs]
    else:
        check_is_lsv(lsvs[0])
    example_LSV = lsvs[0]
    prob_key = get_name_of_prob_key(example_LSV)
    dpsi_header[4] = prob_key
    LSV_copies = [copy_lsv(x) for x in lsvs]
    if "Voila Link" in lsvs[0]:
        has_voila = True
    else:
        has_voila = False
    for lsv_copy in LSV_copies:
        # make dPSIs 0 (unless imputing_with!=0)
        lsv_copy[dpsi_header[3]] = [imputing_with for x in lsv_copy[dpsi_header[3]]]
        # make probs 0 (unless imputing_with!=0)
        lsv_copy[dpsi_header[4]] = [imputing_with for x in lsv_copy[dpsi_header[4]]]
        # make psis 0 (unless imputing_with!=0)
        psi_key1, psi_key2 = get_name_of_psi_keys(lsv_copy)
        lsv_copy[psi_key1] = [imputing_with for x in lsv_copy[psi_key1]]
        lsv_copy[psi_key2] = [imputing_with for x in lsv_copy[psi_key2]]
        # voila link replace
        if has_voila:
            lsv_copy[dpsi_header[21]] = "LSV WAS BLANKED"
    if was_lonely:
        return LSV_copies[0]
    return LSV_copies


def unimpute_lsv_data(Data, BlankDict):
    """
    Undoes the action of initialize_missing_LSVs()
    """
    for comparison in list(Data.keys()):
        for lsv_id in BlankDict[comparison]:
            Data[comparison].pop(lsv_id)


def change_imputed_values(Data, BlankDict, new_val="NA"):
    """

    :param Data: from quick_import -> impute_missing_lsvs
    :param BlankDict: dict returned if in_place argument used in impute_missing_lsvs
    :param new_val: replace the value (likely 0) with a new value (how about NA??)
    :return: nothing .. alters original dictionary
    """
    for comparison in list(Data.keys()):
        for lsv_id in list(BlankDict[comparison]):
            Data[comparison][lsv_id] = impute_lsvs(Data[comparison][lsv_id],
                                                   imputing_with=new_val)


def impute_missing_lsvs(data,
                        impute_with=0,
                        in_place=True,
                        verbose=True,
                        warnings=True):
    """
    Given a Quick import of Data:
    1) Get list of all LSV IDs from all comparisons
    2) Identify which LSV IDs are 'missing' from each comparison
    3) Add missing LSVs to appropriate comparisons:
        a. Copy of an extant LSV ID
        b. Set all junction E(PSI), E(dPSI), and Prob(E(dPSI)) to 0

    impute_with : (default 0) impute missing values with what ?
    InPlace : if True, then the Data is over-written.
    Verbose : print status statements
    Warnings: should warnings be printed?

    Also generates a BlankDict:
        {comparison -> list of LSV IDs that were blanked}

    Returns:
        if InPlace:
            BlankDict
        else:
            Data with blanked LSVs, BlankDict
    """
    if in_place:
        if warnings:
            print("WARNING: YOUR MAJIQ RESULTS WILL BE OVERWRITTERN SINCE InPlace=True")
    check_is_quick_import(data)
    unique_ids = set(get_all_unique_lsv_ids(data))
    new_dict = dict()
    blanked_dict = dict()
    for comparison in list(data.keys()):
        if verbose:
            print("Filling in the LSV gaps for", comparison, "...")
        this_comps_lsvids = set(get_LSV_IDs(data[comparison]))
        only_in_unique = unique_ids - this_comps_lsvids
        blanked_dict[comparison] = only_in_unique
        if len(only_in_unique) == 0:
            print("%s has all the LSVs already!" % comparison)
            if not in_place:
                if verbose:
                    print("Deep copying...")
                new_dict[comparison] = copy.deepcopy(data[comparison])
            continue
        only_in_unique_lsvs = get_lsvs_quickly(data, only_in_unique)
        if verbose:
            print("Imputing...")
        only_in_unique_lsvs_blanked = impute_lsvs(only_in_unique_lsvs, imputing_with=impute_with)
        if verbose:
            print("Finished imputing...")
        only_in_unique_lsvs_blanked_dict = {x["LSV ID"]: x for x in only_in_unique_lsvs_blanked}
        if not in_place:
            # Copy the dict, new object!
            if verbose:
                print("Deep copying...")
            new_dict[comparison] = copy.deepcopy(data[comparison])
            if verbose:
                print("Updating...")
            new_dict[comparison].update(only_in_unique_lsvs_blanked_dict)
        else:
            data[comparison].update(only_in_unique_lsvs_blanked_dict)
    if in_place:
        return blanked_dict
        # return Data
    return new_dict, blanked_dict
