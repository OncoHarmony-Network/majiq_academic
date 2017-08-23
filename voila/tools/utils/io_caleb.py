import re
import traceback

import numpy as np
import os
import pandas as pa
import pdb
import copy
from voila.tools.utils import find_files
from voila.tools.utils.calebs_xrange import calebs_xrange
from voila.utils.voila_log import voila_log

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'

LOG = voila_log()


def quick_import(input,
                 cutoff_d_psi=0.0,
                 cutoff_prob=0.0,
                 keep_ir=True,
                 cutoff_sum=False,
                 pattern="*deltapsi_deltapsi.tsv",
                 deseq_dir=False,
                 deseq_prefix=False,
                 deseq_pat="*csv",
                 deseq_log2fc=1.2,
                 deseq_pv=0.05,
                 deseq_sep="\t",
                 deseq_id_colname="ensembl_gene_id",
                 just_one=False,
                 stop_at=False,
                 comparisons=None,
                 prefered_type=None,
                 just_file_paths=False):
    """
    Given a directory with '*_quantify_deltapsi' files, import all dPSI
        text files and return dictionary as follows:
        {NAME_OF_COMPARISON : {LSV IDs: {LSV INFO: VALUES}}}

    Arguments:
        input: a directory, a path to a list of voila text files, or a voila text file
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
        comparisons: if provided, only import tsv files with the provided list of comparison names
        prefered_type: If provided, only import "deltapsi" or only "psi" files
        just_file_paths: if True, just return the file paths for found voila txt files

    Assumptions:
        Directory contains files that end in ".deltapsi_quantify_deltapsi.txt"
            or DIrectory points directly at a deltapsi.txt file.

        Expected format of input files:
            ~/path/NAME_OF_COMPARISON.deltapsi_quantify_deltapsi.txt

    Returns all voila dPSI results in a dictionary format.
    """
    if prefered_type:
        if not isinstance(prefered_type, str) or prefered_type not in ["deltapsi", "psi"]:
            raise ValueError("prefered_type must be either 'deltapsi' or 'psi' if specified at all, not '%s'" % prefered_type)
    if not os.path.isdir(input):
        LOG.info("Looks like the user provided a file...")
        if name_looks_like_voila_txt_file(input, pattern=pattern):
            LOG.info("It could be a voila tab file ...")
            basename = os.path.basename(input)
            # Get the file name before .psi_psi or before .deltapsi_deltapsi
            basenames = [get_base_names(basename)]
            voila_txt_files = [input]
        elif os.path.isfile(input):
            LOG.info("It could be a file with file paths to voila tab files ...")
            is_valid_list, voila_txt_files = is_likely_list_of_txtfiles(input)
            if is_valid_list:
                # Get the file names before .psi_psi or before .deltapsi_deltapsi
                basenames = [get_base_names(os.path.basename(x)) for x in voila_txt_files]
            else:
                raise ValueError("%s wasn't  valid list of voila txt files.")
        else:
            raise ValueError(input + " not found.")
    else:
        LOG.info("Searching for %s files ..." % pattern)
        basenames, voila_txt_files = find_files.get_voila_files(input,
                                                                pattern=pattern,
                                                                get_base_names=True)
        if len(basenames) != len(voila_txt_files):
            raise ValueError("Something is probably screwy with the names "
                             "of the voila text files...")
        # Preemptively check if the file paths are all actually pointing at valid voila txt files. Remove the baddies.
        basenames_fix = list()
        voila_txt_files_fix = list()
        check_txtfiles = check_valid_voila_txts(voila_txt_files)
        for bn, fp, is_valid in zip(basenames, voila_txt_files, check_txtfiles):
            if is_valid:
                basenames_fix.append(bn)
                voila_txt_files_fix.append(fp)
        # too lazy to rename things
        basenames = basenames_fix
        voila_txt_files = voila_txt_files_fix
        # Only want to import stuff in list of comparisons...
        if comparisons:
            if isinstance(comparisons, str):
                comparisons = [comparisons]
            elif not isinstance(comparisons, list):
                LOG.error("User provided comparisons is not a list :(")
                exit(1)
            to_remove = list()
            to_keep = list()
            for ii in range(len(voila_txt_files)):
                if not basenames[ii] in comparisons:
                    to_remove.append(ii)
                else:
                    to_keep.append(ii)
            if len(to_keep) < len(comparisons):
                LOG.error("Couldn't find all these comparisons you wanted: %s" % comparisons)
            for index in sorted(to_remove, reverse=True):
                del basenames[index]
                del voila_txt_files[index]
    if len(set(basenames)) != len(set(voila_txt_files)):
        # this means more than one tsv with the same name is in this directory substructure
        # possibly because user ran different thresholds (e.g. 10% and 20%) but kept the same names
        # no worries, I think I have a fix... maybe... this could come back and bite me in the @$$
        dupes = [x for n, x in enumerate(basenames) if x in basenames[:n]]
        for dup in dupes:
            i = 0
            while dup in basenames:
                new_name = "%s_duplicate%s" % (dup, i)
                basenames[basenames.index(dup)] = new_name
                i += 1

    if len(voila_txt_files) == 0:
        raise RuntimeError("Didn't find any voila txt files...")

    LOG.info("Found " + str(len(voila_txt_files)) +
             " dPSI text files ...")

    imported_files = dict()
    funky_ids = list()

    if just_one:
        voila_txt_files = [voila_txt_files[0]]
    valid_fps = list()
    for f, comparison_name in zip(voila_txt_files, basenames):
        if prefered_type:
            expected = prefered_type
        else:
            expected = "deltapsi" if "deltapsi_deltapsi" in os.path.basename(f) else "psi"
        imported_file = import_voila_txt(f,
                                         cutoff_d_psi,
                                         cutoff_prob,
                                         stop_at=stop_at,
                                         expected_type=expected,
                                         just_checking_validity=just_file_paths)

        # imported_file will be False if import_dpsi thinks it is not a valid tsv...
        if not imported_file:
            continue
        if just_file_paths:
            valid_fps.append(f)
        else:
            imported_files[comparison_name] = imported_file
    if just_file_paths:
        return valid_fps
    if cutoff_d_psi != 0 or cutoff_prob != 0 or not keep_ir:
        imported_files = subset_significant(imported_files,
                                            cutoff_dpsi=cutoff_d_psi,
                                            cutoff_prob=cutoff_prob,
                                            keep_introns=keep_ir,
                                            cutoff_sum=cutoff_sum)
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


def check_valid_voila_txts(file_path_list):
    """
    For each file in the list, check if it is a valid voila txt file, return list of bools
    :param file_path_list:
    :return: list of bools
    """
    passed = list()
    for fp in file_path_list:
        expected = "deltapsi" if "deltapsi_deltapsi" in os.path.basename(fp) else "psi"
        test_imp = import_voila_txt(fp,
                                    stop_at="Gene Name",  # my trick to not actually import the file ;)
                                    expected_type=expected,
                                    just_checking_validity=True)
        if test_imp:
            passed.append(True)
        else:
            passed.append(False)
    return passed


def import_voila_txt(fp,
                     cutoff_d_psi=0,
                     cutoff_prob=0,
                     stop_at=False,
                     expected_type="deltapsi_deltapsi",
                     just_checking_validity=False):
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
        stop_at: if provided, stop reading voila file when you reach this LSV ID
        expected_type: If provided, only import "deltapsi" or only "psi" files
        just_checking_validity: Bool. If True, don't print "Importing..."


    """
    if not isinstance(stop_at, list) and not isinstance(stop_at, str):
        raise ValueError("%s needs to be a list or str, not %s" % (stop_at, type(stop_at)))
    if not isinstance(fp, str):
        raise TypeError("Expected file path to be string, instead it was %s" % type(fp))
    if not have_permission(fp):
        LOG.info("Uhhh this is embarrassing, it looks like you don't have permissions to view"
                 " this file: %s" % fp)
        return False
    if check_if_file_binary(fp):
        LOG.info("%s matched search pattern, but it is binary, so the file is not"
                 " a voila tsv ... skipping it" % fp)
        return False
    lsv_dictionary = dict()
    has_voila = True
    file_headers = list()
    with open(fp, "r") as handle:
        line_i = 0
        found_stop_at = False if isinstance(stop_at, str) else [False for x in stop_at]
        can_stop = False
        for line in handle:
            # if "ENSG00000078061" in line:
            #     pdb.set_trace()
            if isinstance(stop_at, str):
                if stop_at in line:
                    found_stop_at = True
                    can_stop = True
                else:
                    found_stop_at = False
                    if line_i > 0 and not can_stop:
                        line_i += 1
                        continue
                    elif line_i > 0 and can_stop:
                        break
            elif isinstance(stop_at, list):
                gene_in_line_bools = [gene in line for gene in stop_at]
                if True in gene_in_line_bools:
                    gene_ii = gene_in_line_bools.index(True)
                    found_stop_at[gene_ii] = True
                else:
                    if line_i > 0 and not can_stop:
                        line_i += 1
                        continue
                    elif line_i > 0 and can_stop:
                        break
            line_split = line.rstrip("\r\n").split("\t")
            if line_i == 0:
                if expected_type == "deltapsi":
                    if not has_valid_voila_dpsi_tsv_header(line):
                        LOG.info("%s matched search pattern, but this file doesn't appear"
                                 " to be a valid voila deltapsi_deltapsi output... skipping it..." % fp)
                        return False
                else:
                    if not has_valid_voila_psi_tsv_header(line):
                        LOG.info("%s matched search pattern, but this file doesn't appear"
                                 " to be a valid psi_psi tsv output... skipping it..." % fp)
                        return False
                # Fix pound sign silliness
                line_split[0] = line_split[0].replace("#", "")
                file_headers.extend(line_split)
                if expected_type == "deltapsi":
                    condition_1_name = line_split[5].split(" ")[0]
                    condition_2_name = line_split[6].split(" ")[0]
                    if not just_checking_validity:
                        threshold = get_threshold(line_split)
                        LOG.info("Importing %s vs %s deltapsi data (P>=%s) ..." % (condition_1_name,
                                                                                   condition_2_name,
                                                                                   threshold))
                # Else its a psi file
                else:
                    sample_id = get_base_names(fp)
                    if not just_checking_validity:
                        LOG.info("Importing %s psi_psi data ..." % (sample_id))
                if "Voila Link" in file_headers:
                    has_voila = True
                    pre_voila_1_0_0 = True
                else:
                    pre_voila_1_0_0 = False
                    has_voila = False

                line_i += 1
                continue

            if expected_type == "deltapsi":
                the_data = get_dpsi_data(line_split,
                                         pre_voila_1_0_0,
                                         file_headers,
                                         has_voila)
            else:
                the_data = get_psi_data(line_split,
                                        pre_voila_1_0_0,
                                        file_headers,
                                        has_voila)

            line_i += 1
            if can_stop:
                if isinstance(stop_at, str):
                    if not found_stop_at:
                        break
                else:
                    if False not in found_stop_at:
                        break
            # Add the line's data to the dict
            lsv_dictionary.update(the_data)
            if just_checking_validity:
                return True
            if isinstance(stop_at, list):
                if False not in found_stop_at:
                    can_stop = True
                else:
                    can_stop = False
    # add the meta_info for the experiment
    lsv_dictionary["meta_info"] = dict()
    lsv_dictionary["meta_info"]["abs_path"] = os.path.abspath(fp)
    if expected_type == "deltapsi":
        lsv_dictionary["meta_info"]["condition_1_name"] = condition_1_name
        lsv_dictionary["meta_info"]["condition_2_name"] = condition_2_name
        if not just_checking_validity:
            lsv_dictionary["meta_info"]["prob_thresh"] = threshold
    else:
        lsv_dictionary["meta_info"]["sample_id"] = sample_id
    return lsv_dictionary


def get_threshold(line_split):
    """
    From split dpsi voila text file header, return the probabily threshold as a float
    :param line_split: the first line of the file
    :return: float
    """
    try:
        the_thresh = float(line_split[4].split("P(|dPSI|>=")[1].split(") per LSV junction")[0])
    except:
        err_message = ""
        if "P(|dPSI|>=" not in line_split[4]:
            err_message += "...Uh oh, 'P(|dPSI|>=' wasn't in the prob column header..."
        if ") per LSV junction" not in line_split[4]:
            err_message += "\n...Uh oh, ') per LSV junction' wasn't in the prob column header..."
        if len(err_message) > 0:
            LOG.error(err_message)
        raise
    return the_thresh


def get_dpsi_data(line_split,
                  pre_voila_1_0_0,
                  file_headers,
                  has_voila):
    """

    :param line_split:
    :return:
    """
    lsv_dictionary = dict()
    gene_name = str(line_split[0])

    gene_id = str(line_split[1])

    d_psi_floated = [float(x) for x in line_split[3].split(';')]

    prob_d_psis_floated = [float(x) for x in line_split[4].split(';')]

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

    if pre_voila_1_0_0:
        de_novo_junct = int(line_split[13])
    else:  # Else it is boolean
        de_novo_junct = str(line_split[13])

    chrom = str(line_split[14])

    strand = str(line_split[15])

    junct_coord = str(line_split[16]).split(";")

    exon_coord = str(line_split[17]).split(";")

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

    lsv_dictionary[LSV_ID]["Reference_Type"] = ref_type
    return lsv_dictionary


def get_psi_data(line_split,
                 pre_voila_1_0_0,
                 file_headers,
                 has_voila):
    """

    :param split_line:
    :return:
    """
    lsv_dictionary = dict()
    gene_name = str(line_split[0])

    gene_id = str(line_split[1])

    psi_floated = [float(x) for x in line_split[3].split(';')]

    var_psi_float = [float(x) for x in line_split[4].split(';')]

    LSV_ID = str(line_split[2])
    if "target" in LSV_ID:
        ref_type = "target"
    elif "source" in LSV_ID:
        ref_type = "source"

    lsv_type = str(line_split[5])

    a5ss = bool(line_split[6])

    a3ss = bool(line_split[7])

    es = bool(line_split[8])

    n_junctions = int(line_split[9])

    n_exons = int(line_split[10])

    if pre_voila_1_0_0:
        de_novo_junct = int(line_split[11])
    else:  # Else it is boolean
        de_novo_junct = str(line_split[11])

    chrom = str(line_split[12])

    strand = str(line_split[13])

    junct_coord = str(line_split[14]).split(";")

    exon_coord = str(line_split[15]).split(";")

    exon_alt_start = str(line_split[16])

    exon_alt_end = str(line_split[17])

    ir_coords = str(line_split[18])

    if has_voila:
        voila_link = str(line_split[19])

    lsv_dictionary[LSV_ID] = dict({
        file_headers[0]: gene_name,
        file_headers[1]: gene_id,
        file_headers[2]: LSV_ID,
        file_headers[3]: psi_floated,
        file_headers[4]: var_psi_float,
        file_headers[5]: lsv_type,
        file_headers[6]: a5ss,
        file_headers[7]: a3ss,
        file_headers[8]: es,
        file_headers[9]: n_junctions,
        file_headers[10]: n_exons,
        file_headers[11]: de_novo_junct,
        file_headers[12]: chrom,
        file_headers[13]: strand,
        file_headers[14]: junct_coord,
        file_headers[15]: exon_coord,
        file_headers[16]: exon_alt_start,
        file_headers[17]: exon_alt_end,
        file_headers[18]: ir_coords})
    if has_voila:
        lsv_dictionary["Voila Link"] = voila_link

    lsv_dictionary[LSV_ID]["Reference_Type"] = ref_type
    return lsv_dictionary


def is_likely_list_of_txtfiles(the_file):
    """
    Check if the_file is a line-by-line list of voila text files
    :param the_file:
    :return: True, [list, of, files] OR False, [list of files found until oe isn't a txt file]
    """
    poss_files = list()
    if not os.path.isfile(the_file):
        raise ValueError("Sorry, this doesn't appear to be a valid filepath:\n%s" % the_file)
    if not have_permission(the_file):
        raise ValueError("Sorry, you don't have permission to open:\n%s" % the_file)
    if is_voila_txt_file(the_file):
        return False, poss_files
    if check_if_file_binary(the_file):
        rtype = "rb"
    else:
        rtype = "r"
    with open(the_file, rtype) as handle:
        for line in handle:
            poss_file = line.rstrip("\n\r")
            if len(poss_file) == 0:
                continue
            if not is_voila_txt_file(poss_file):
                LOG.info("This doesn't look like a voila txt file:\n%s" % poss_file)
                return False, poss_files
            else:
                poss_files.append(poss_file)
    if len(poss_files) == 0:
        LOG.info("No voila text files found in:\n%s" % the_file)
        return False, poss_files
    return True, poss_files


def is_voila_txt_file(the_file):
    """
    Check if the file is a voila txt file
    :param the_file:
    :return: Bool
    """
    if not os.path.isfile(the_file):
        LOG.error("Sorry, first line of supposed list doesn't appear to be a valid filepath:\n%s" % the_file)
        exit(1)
    if not have_permission(the_file):
        LOG.error("Sorry, you don't have permission to open:\n%s" % the_file)
        exit(1)
    if check_if_file_binary(the_file):
        LOG.info("This is binary, so it's definitely not a voila text file:\n%s" % the_file)
        return False
    looks_like_voila = False
    line_i = 0
    with open(the_file, "r") as handle:
        for line in handle:
            if line_i == 0:
                if has_valid_voila_dpsi_tsv_header(line):
                    looks_like_voila = True
            break
    return looks_like_voila


def have_permission(file_path, f_type="r"):
    """
    Check if you have user permission to read file
    :param f_type: "r" or "rb"
    :param file_path:
    :return: Bool
    """
    try:
        i = 0
        with open(file_path, f_type) as handle:
            for line in handle:
                i += 1
                if i > 10:
                    break
    except UnicodeDecodeError:
        return have_permission(file_path=file_path, f_type="rb")
    except PermissionError:
        return False
    return True


def check_if_file_binary(file_path):
    """
    https://stackoverflow.com/a/7392391/7378802
    :param file_path:
    :return: Bool
    """
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f})
    the_bytes = open(file_path, 'rb').read(1024)
    if bool(the_bytes.translate(None, textchars)):
        return True
    try:
        i = 0
        with open(file_path, "r") as handle:
            for line in handle:
                i += 1
                if i > 10:
                    break
    except UnicodeDecodeError:
        return True
    return False


def has_valid_voila_dpsi_tsv_header(header_line):
    """
    Make sure the supposed tsv file is actually a tsv file
    :param header_line: str, from file's first line
    :return: True or False
    """
    is_valid = True
    expect_to_be_in = ["E(dPSI) per LSV junction",
                       "#Gene Name",
                       "Gene ID",
                       "LSV ID",
                       "E(PSI)",
                       "LSV Type",
                       "Exons coords",
                       "Junctions coords",
                       "strand",
                       "chr"]
    for expected in expect_to_be_in:
        if expected not in header_line:
            return False
    return is_valid


def has_valid_voila_psi_tsv_header(header_line):
    """
    Make sure the supposed tsv file is actually a tsv file
    :param header_line: str, from file's first line
    :return: True or False
    """
    is_valid = True
    expect_to_be_in = ["E(PSI) per LSV junction",
                       "Var(E(PSI)) per LSV junction",
                       "#Gene Name",
                       "Gene ID",
                       "LSV ID",
                       "E(PSI)",
                       "LSV Type",
                       "Exons coords",
                       "Junctions coords",
                       "strand",
                       "chr"]
    for expected in expect_to_be_in:
        if expected not in header_line:
            return False
    return is_valid


def comp_without_dup(comp_name):
    """
    If there were duplicate tsv files with the same comparison name,
        I add duplicate_# to the end. This will return the name without
        the duplicate_#, or just the name if it doesn't have duplicate_
    :param comp_name:
    :return:
    """
    fixed_name = comp_name
    if "_duplicate" in comp_name:
        fixed_name = comp_name[0:comp_name.index("_duplicate")]
    return fixed_name


def subset_significant(data,
                       cutoff_dpsi=0.2,
                       cutoff_prob=0.95,
                       keep_introns=True,
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
        over_cutoff_ids = get_sig_lsv_ids(data,
                                          cutoff_dpsi,
                                          cutoff_prob,
                                          cutoff_sum)
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


def get_base_names(file):
    if os.path.exists(file):
        basename = os.path.basename(file)
    else:
        basename = file
    split_file_name = basename.split(".")
    comparison_name = split_file_name[0]
    return comparison_name


def name_looks_like_voila_txt_file(the_file, pattern="*deltapsi_deltapsi.tsv"):
    if not os.path.exists(the_file):
        raise ValueError("Supposed deltapsi txt file doesn't exist.")
    newpat = pattern.replace("*", ".*")
    if re.search(newpat, the_file):
        if the_file[:-4].endswith("deltapsi_deltapsi") or the_file[:-4].endswith("psi_psi"):
            return True
        else:
            return False
    else:
        return False


def recursive_dirs(Directory):
    """
    Given a directory, return all sub directories, recursively.
    """
    is_dir = list()
    for thing in os.listdir(Directory):
        thing = os.path.join(Directory, thing)
        if os.path.isdir(thing):
            is_dir.append(thing)
            is_dir.extend(recursive_dirs(thing))
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
    names = list(LSV_dict.keys())
    non_intron_names = list()
    intron_names = list()

    for name in names:
        if name == "meta_info":
            # non_intron_names.append(name)
            # intron_names.append(name)
            continue
        # If LSV type is intron:
        if LSV_dict[name]["LSV Type"][-1:] == "i":
            # Save these names, too, just in case I want em later
            intron_names.append(name)
        else:
            non_intron_names.append(name)

    for intron in intron_names:
        # if intron == "condition_1_name" or intron == "condition_2_name":
        #     continue
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
    deseq_fps = find_files.find_files(Path=deseq_dir, Pattern=deseq_fname_pattern, Recursive=recursive)
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
    deseq_fps = find_files.find_files(Path=DESeqDirectory, Pattern=Pattern, Recursive=Recursive)
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


def check_is_ignant(data, users_prob):
    """
    You should only use P(Expected dPSI) to ID significantly changing LSVs at the
        threshold voila was run with. So, only 0 or exactly the threshold. Raise
        error if this isn't the case for any of the data
    :param data:
    :param users_prob:
    """
    comps = get_comparisons(data, sort=True)
    lsv_dictlist = [data[comp] for comp in comps]
    all_threshes = [get_prob_threshold(thed) for thed in lsv_dictlist]
    for comp, tsv_thresh in zip(comps, all_threshes):
        if users_prob != tsv_thresh and users_prob != 0:
            raise ValueError("WARNING !!! You are using an ill-advised threshold (%s) "
                             " ... %s was run with voila thresh of %s!" % (users_prob,
                                                                           comp,
                                                                           tsv_thresh))


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
        lsv_ids.remove("meta_info")
        # lsv_ids.remove("condition_1_name")
        # lsv_ids.remove("condition_2_name")
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


def psi_or_deltapsi(the_thing):
    if check_is_lsv_dict(the_thing,
                         da_bool=True):
        if "sample_id" in the_thing["meta_info"]:
            return "psi"
        elif "condition_1_name" in the_thing["meta_info"]:
            return "deltapsi"
        else:
            LOG.error("No sample_id or condition_1_name... are you sure this is an LSV dict??")
            exit(1)
    if check_is_lsv(the_thing,
                    Bool=True):
        if "E(PSI) per LSV junction" in the_thing:
            return "psi"
        elif "E(dPSI) per LSV junction" in the_thing:
            return "deltapsi"
        else:
            LOG.error("Not psi or deltapsi... are you sure this is an LSV??")
            exit(1)


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
                   "Exons coords"]:
        if header not in lsv:
            if not Bool:
                raise ValueError("Expected a LSV, but didn't get one.")
            else:
                return False
    keys = list(lsv.keys())
    sub_dict = lsv[keys[0]]
    if isinstance(sub_dict, dict):
        keys = list(sub_dict.keys())
        # if "condition_1_name" in keys:
        if "meta_info" in keys:
            if not Bool:
                raise ValueError("Please pick one comparison from the dictionary, "
                                 "not the entire dictionary returned by quick_import()."
                                 " For entire dictionary, use lookup_everywhere() "
                                 "instead of lookup().")
            else:
                return False
    return True


def quick_import_subset(data,
                        lsv_ids,
                        in_place=False):
    """
    Assuming the data is imputed or shares all lsv ids
    :param data: quick import (maybe imputed)
    :param lsv_ids:
    :param in_place: if True, overwerite data
    :return:
    """
    comparisons = get_comparisons(data, sort=True)
    res = dict()
    for comp in comparisons:
        this_subset = lsv_dict_subset(data[comp], lsv_ids)
        if in_place:
            data[comp] = this_subset
        else:
            res[comp] = this_subset
    if not in_place:
        return res


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
        # if "condition_1_name" in dictionary:
        if "meta_info" in dictionary:
            # new_dict["dPSI_over_cutoff"]=Dictionary["dPSI_over_cutoff"]
            # new_dict["sig_junctions"]=Dictionary["sig_junctions"]
            # new_dict["condition_1_name"] = dictionary["condition_1_name"]
            # new_dict["condition_2_name"] = dictionary["condition_2_name"]
            new_dict["meta_info"] = dictionary["meta_info"]
    if None in new_dict.values():
        LOG.info("Warning: at least 1 key wasn't found in the dictionary."
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
        shared_lsvs = get_shared_lsv_ids(data, bool=True)
        if shared_lsvs:
            n_shared = len(shared_lsvs)
        else:
            n_shared = 0
        LOG.info("%s Shared LSVs between all comparisons." % (n_shared))
        LOG.info("Total of %s unique LSVs across all comparisons." % (n_all))


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
        all_lsv_ids = get_lsv_ids(rm_data[comp])
        n_s = len(sig_gene_ids)
        n_t = len(all_lsv_ids)
        n_r = 0
        removed_ids = []
        for lsv_id in all_lsv_ids:
            gene = lsv_id.split(":")[0]
            if gene in sig_gene_ids:
                n_r += 1
                removed_ids.append(rm_data[comp].pop(lsv_id))
        LOG.info("%s: %s LSVs (%s diff-expr genes) removed leaving %s LSVs" % (comp, n_r, n_s, n_t - n_r))
    data_comps = set(rm_data.keys())
    comps_with_deseq_res = set(comps_with_deseq_res)
    leftover = list(comps_with_deseq_res - data_comps)
    if len(leftover) > 0:
        LOG.info("The following comparisons didn't have DESeq results:\n%s" % leftover)
    if return_diff_expr_lsvids:
        return removed_ids


def get_comparisons(data, sort=True):
    """

    :param data: quick import
    :param sort: return names sorted? True or False
    :return: list of comparison names
    """
    check_is_quick_import(data)
    comps = list(data.keys())
    if sort:
        comps.sort()
    return comps


def get_sig_lsv_ids(data,
                    cutoff_d_psi=0.0,
                    prob_d_psi=0.0,
                    sum_for_cutoff=False,
                    collapse=False):
    """
    Given LSV dictionary, return set of unique LSV IDs over cutoff

        Note: recursively handles Quick_import and returns dict of
        comparison_name -> sig_ids set, unless collapse=True

    Arguments:
        data: output of import_dpsi
        cutoff_d_psi: Only return LSV IDs with at least 1 junction abs(dPSI) >= Cutoff_dPSI
        prob_d_psi: junction must have a dPSI>= this to be considered
        sum_for_cutoff: If the sum of all +dPSI (that meet probability cutoff)
            is >= Cutoff_dPSI, then include the LSV. Or if abs(sum of dPSI<0) is >= Cutoff_dPSI.
            This is less conservative than default.
        Collapse: if data is quick import, and collapse=True, collapse all sets into
            one big set to return

    Return:
        set or dict of sets
    """
    if check_is_quick_import(data, the_bool=True):
        comparisons = data.keys()
        results = dict()
        for comparison in comparisons:
            LSV_Dict = data[comparison]
            results[comparison] = get_sig_lsv_ids(LSV_Dict,
                                                  cutoff_d_psi,
                                                  prob_d_psi,
                                                  sum_for_cutoff)
        if collapse:
            collapsed_res = set()
            for this_set in results.values():
                collapsed_res = collapsed_res | this_set
            results = collapsed_res
        return results
    check_is_lsv_dict(data)

    # names AKA LSV IDs
    names = get_lsv_ids(data)
    names_over_cutoff = set()
    if len(names) < 1:
        raise RuntimeError("No LSVs made Cutoff dPSI of %s and Prob of %s" % (cutoff_d_psi,
                                                                              prob_d_psi))
    prob_name = get_name_of_prob_key(data[names[0]])
    for name in names:
        dPSIs = data[name]["E(dPSI) per LSV junction"]
        prob_dPSIs = data[name][prob_name]
        if sum_for_cutoff:
            sum_dPSI_over = 0
            sum_dPSI_under = 0
            meets_prob_cutoff = False
            for dPSI, prob_dPSI in zip(dPSIs, prob_dPSIs):
                if prob_dPSI < prob_d_psi:
                    continue
                else:
                    meets_prob_cutoff = True
                if dPSI > 0:
                    sum_dPSI_over += dPSI
                elif dPSI < 0:
                    sum_dPSI_under -= dPSI
            if sum_dPSI_over >= cutoff_d_psi or abs(sum_dPSI_under) >= cutoff_d_psi:
                if meets_prob_cutoff:
                    names_over_cutoff.add(name)
        else:
            for dPSI, prob_dPSI in zip(dPSIs, prob_dPSIs):
                if (abs(dPSI) >= cutoff_d_psi) and (prob_dPSI >= prob_d_psi):
                    names_over_cutoff.add(name)
    return names_over_cutoff


def get_lsv_ids(lsv_dict):
    """
    Return LSV IDs from dictionary.
    """
    check_is_lsv_dict(lsv_dict)
    lsv_ids = copy.copy(list(lsv_dict.keys()))
    # lsv_ids.remove("condition_1_name")
    # lsv_ids.remove("condition_2_name")
    lsv_ids.remove("meta_info")
    return lsv_ids


def get_dpsis(lsv,
              prob_cutoff=None,
              as_np_array=False):
    """
    Given LSV, return list of dPSIs over Prob_Cutoff, if provided.
        If not Prob_cutoff provided, return all dPSIs

        as_np_array: bool
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
        return np.array(dPSI) if as_np_array else dPSI
    else:
        return np.array(all_dPSI) if as_np_array else all_dPSI


def get_psis(lsv, cond_1=False, cond_2=False, as_dict=False, as_np_array=False):
    """
    :param lsv:
    :param cond_1: string, optional
    :param cond_2: string, optional
    :param as_dict: if True, return as dict with cond:psis
    :param as_np_array: if True, return as np arrays instead of lists
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
    if as_np_array:
        cond_1_psi = np.array(cond_1_psi)
        cond_2_psi = np.array(cond_2_psi)
    if as_dict:
        return {psi_key_1.split(" ")[0]: cond_1_psi, psi_key_2.split(" ")[0]: cond_2_psi}
    return [cond_1_psi, cond_2_psi]


def genename_from_id(lsvdict, lsvid):
    """
    Given an LSV dictionaryor quick import  and an lsv id, return the gene name
    :param lsvdict:
    :param lsvid:
    :return: str
    """
    if check_is_lsv_dict(lsvdict, da_bool=True):
        return get_gene_name(lsvdict[lsvid])
    check_is_quick_import(lsvdict)
    for comp in lsvdict:
        for thislsvid in lsvdict[comp]:
            if lsvid in thislsvid:
                return get_gene_name(lsvdict[comp][thislsvid])
    raise ValueError("%s not found." % lsvid)


def get_gene_name(lsv):
    check_is_lsv(lsv)
    return lsv["Gene Name"]


def get_strand(lsv):
    check_is_lsv(lsv)
    return copy.copy(lsv["strand"])


def get_chr(lsv):
    check_is_lsv(lsv)
    return copy.copy(lsv["chr"])


def get_probs(lsv,
              as_np_array=False):
    """
    Given LSV, return P(E(dPSI))
    """
    check_is_lsv(lsv)
    res = copy.copy(lsv[get_name_of_prob_key(lsv)])
    if as_np_array:
        res = np.array(res)
    return res


def get_juncs(lsv):
    """
    Given LSV, return exons coords
    """
    check_is_lsv(lsv)
    return copy.copy(lsv["Junctions coords"])


def get_exons(lsv):
    """
    Given LSV, return exons coords
    """
    check_is_lsv(lsv)
    return copy.copy(lsv["Exons coords"])


def list_d_psis(lsv_dict,
                as_np_array=False):
    """
    Return dictionary of LSV IDs pointing at list of dPSIs:
    (this col
    not in dict) list:
    Junction:   | dPSI |
            0   |   #  |
            1   |   #  |
            2   |   #  |
            ... |   ...|
    """
    check_is_lsv_dict(lsv_dict)
    lsv_to_psi_dict = dict()
    lsvs = get_lsv_ids(lsv_dict)
    # Extract dPSI from each LSV, using cutoff.
    for lsv_id in lsvs:
        lsv_to_psi_dict[lsv_id] = get_dpsis(lsv_dict[lsv_id],
                                            as_np_array=as_np_array)
    return lsv_to_psi_dict


def list_probs(lsv_dict,
               as_np_array=False):
    """
    Return dictionary of LSV IDs pointing at list of P(dPSIs):
    (this col
    not in dict) list (or array):
    Junction:   | Prob |
            0   |   #  |
            1   |   #  |
            2   |   #  |
            ... |   ...|
    """
    check_is_lsv_dict(lsv_dict)
    lsv_to_prob_dict = dict()
    lsvs = get_lsv_ids(lsv_dict)
    # Extract dPSI from each LSV, using cutoff.
    for lsv_id in lsvs:
        lsv_to_prob_dict[lsv_id] = get_probs(lsv_dict[lsv_id],
                                             as_np_array=as_np_array)
    return lsv_to_prob_dict


def list_psi(lsv_dict,
             as_np_array=False):
    """
    Return dictionary of conditions pointing at dictionary of
        LSV IDs pointing at list of PSIs:

    condition_* points at:
    (this col
    not in dict) list:
    Junction:   | PSI  |
            0   |   #  |
            1   |   #  |
            2   |   #  |
            ... |   ...|
    """
    check_is_lsv_dict(lsv_dict)
    # cond_1_name = lsv_dict["condition_1_name"]
    # cond_2_name = lsv_dict["condition_2_name"]
    cond_1_name = get_cond_1_name(lsv_dict)
    cond_2_name = get_cond_2_name(lsv_dict)
    condtion_dict = dict()
    condtion_dict[cond_1_name] = dict()
    condtion_dict[cond_2_name] = dict()
    lsvs = get_lsv_ids(lsv_dict)
    # Extract dPSI from each LSV, using cutoff.
    for lsv_id in lsvs:
        PSIs = get_psis(lsv_dict[lsv_id], as_np_array=as_np_array)
        cond_1_psi = PSIs[0]
        cond_2_psi = PSIs[1]
        condtion_dict[cond_1_name][lsv_id] = cond_1_psi
        condtion_dict[cond_2_name][lsv_id] = cond_2_psi
    return condtion_dict


def get_all_dpsis(lsvs,
                  prob_cutoff=0,
                  as_np_arrays=False):
    """
    Given LSV dictionary or quick_import structure, return list of all dPSIs from all
        LSVs over Prob_Cutoff.
    """
    all_dpsi = list()
    try:
        check_is_lsv_dict(lsvs)
        lsvs = get_lsv_ids(lsvs)
        # Extract dPSI from each LSV, using cutoff.
        for lsv_id in lsvs:
            all_dpsi.extend(get_dpsis(lsvs[lsv_id], prob_cutoff, as_np_arrays))
    except:
        try:
            check_is_quick_import(lsvs)
            for key in lsvs.keys():
                lsv_dict = lsvs[key]
                all_dpsi.extend(get_all_dpsis(lsv_dict, prob_cutoff, as_np_arrays))
        except:
            raise ValueError("Expected LSV dictionary or quick_import return value.")
    return all_dpsi


def get_num_d_psi(data,
                  return_comparisons=False,
                  use_binary_index_info=None):
    """
    Given dictionary of LSV dictionaries return numpy arrays
        of all dPSIs from all comparisons for each LSV. Such that
        each LSV ID in the dictionary is pointing at an array that
        looks like this:

                |A vs B|A vs C  | all other comparisons [IN SORTED ORDER]
    Junction:   | dPSI | dPSI   | ...
            0   |   #  |  #     | ...
            1   |   #  |  #     | ...
            2   |   #  |  #     | ...
            ... |   ...|  ...   | ...

    Arguments:
        data: quick imp
        return_comparisons: Boolean. If True, also return a list
            of the comparison names, in the same order as the columns
            of the numpy array.
        use_binary_index_info: None, "closer", or "further"
            If the LSV dictionaries were returned by get_binary_LSVs(),
            then they will have a "binary_indices" key that points at
            the binary indices. Use this info to return only the dPSI
            value that corresponds to the closer or further junction
            from the reference exon.

    NOTES:
        - only those LSV IDs that are shared by all LSV dictionaries
        in Data are evaluated. LSV IDs that are unique to or missing from
        any of the LSV Dictionaries are not returned by this function.
        - columns are sorted by name


    Returns the numpy array and a list of LSV_IDs
    """
    binary_index = "cow"  # stupid pep
    if use_binary_index_info:
        if not use_binary_index_info == "closer" and not use_binary_index_info == "further":
            raise ValueError("Use_binary_index_info needs to be 'closer' or 'further' if provided")
    if use_binary_index_info == "closer":
        binary_index = 0
    if use_binary_index_info == "further":
        binary_index = 1
    check_is_quick_import(data)
    comparison_dict = dict()
    for comparison in data.keys():
        d_psis = list_d_psis(data[comparison])
        comparison_dict[comparison] = d_psis
    union_of_lsv_ids = get_shared_lsv_ids(data)
    lsv_dict = dict()
    comparisons = list(comparison_dict.keys())
    comparisons.sort()  # alphabatized!!!
    for lsv_id in list(union_of_lsv_ids):
        list_of_dpsis = list()
        for comparison in comparisons:
            d_psis = comparison_dict[comparison][lsv_id]
            if use_binary_index_info:
                binary_i = data[comparison][lsv_id]["binary_indices"][binary_index]
                d_psis = d_psis[binary_i]
            list_of_dpsis.append(d_psis)
        lsv_dict[lsv_id] = np.array(list_of_dpsis).T
    if return_comparisons:
        return lsv_dict, comparisons
    return lsv_dict


def get_num_prob(data,
                 return_comparisons=False,
                 use_binary_index_info=False):
    """
    Given dictionary of LSV dictionaries return numpy arrays
        of all P(dPSIs) from all comparisons for each LSV. Such that
        each LSV ID in the dictionary is pointing at an array that
        looks like this:

                |A vs B|A vs C  | all other comparisons [ SORTED!!!]...
    Junction:   | P(dPSI) | P(dPSI)   | ...
            0   |   #     |     #     | ...
            1   |   #     |     #     | ...
            2   |   #     |     #     | ...
            ... |   ...   |     ...   | ...
    Arguments:
        data: Quick Import structure
        return_comparisons: Boolean. If True, also return a list
            of the comparison names, in the same order as the columns
            of the numpy array.
        use_binary_index_info: None, "closer", or "further"
            If the LSV dictionaries were returned by get_binary_LSVs(),
            then they will have a "binary_indices" key that points at
            the binary indices. Use this info to return only the P(dPSI)
            value that corresponds to the closer or further junction
            from the reference exon.

    NOTES:
        - only those LSV IDs that are shared by all LSV dictionaries
        in Data are evaluated. LSV IDs that are unique to or missing from
        any of the LSV Dictionaries are not returned by this function.
        - columns are sorted by name


    Returns the numpy array and a list of LSV_IDs
    """
    # stupid pep rules
    binary_index = "cow"
    if use_binary_index_info:
        if not use_binary_index_info == "closer" and not use_binary_index_info == "further":
            raise ValueError("Use_binary_index_info needs to be 'closer' or 'further' if provided")
        if use_binary_index_info == "closer":
            binary_index = 0
        elif use_binary_index_info == "further":
            binary_index = 1
        else:
            raise RuntimeError("I don't know what to do here")
    check_is_quick_import(data)
    comparison_dict = dict()
    for nup_comparison in data.keys():
        nu_probs = list_probs(data[nup_comparison])
        comparison_dict[nup_comparison] = nu_probs
    union_of_lsv_ids = get_shared_lsv_ids(data)
    lsv_dict = dict()
    comparisons = list(comparison_dict.keys())
    comparisons.sort()  # alphabatized
    for lsv_id in list(union_of_lsv_ids):
        list_of_probs = list()
        for nup_comparison in comparisons:
            probs = comparison_dict[nup_comparison][lsv_id]
            if use_binary_index_info:
                binary_i = data[nup_comparison][lsv_id]["binary_indices"][binary_index]
                probs = probs[binary_i]
            list_of_probs.append(probs)
        lsv_dict[lsv_id] = np.array(list_of_probs).T
    if return_comparisons:
        return lsv_dict, comparisons
    return lsv_dict


def get_num_psi(data,
                return_comparisons=False,
                use_binary_index_info=None):
    """
    Given dictionary of LSV dictionaries return numpy array
        of all PSIs from all conditions for each LSV. Such that
        each LSV ID in the dictionary is pointing at an array that
        looks like this:

                |   A  |  B     | all other conditions [SORTED]...
    Junction:   |  PSI |  PSI   | ...
            0   |   #  |  #     | ...
            1   |   #  |  #     | ...
            2   |   #  |  #     | ...
            ... |   ...|  ...   | ...

    Arguments:
        return_comparisons: if True, also return a list of conditions
            found in the same order as the numpy array columns.
        use_binary_index_info: None, "closer", or "further"
            If the LSV dictionaries were returned by get_binary_LSVs(),
            then they will have a "binary_indices" key that points at
            the binary indices. Use this info to return only the PSI
            value that corresponds to the closer or further junction
            from the reference exon.

    NOTE: only those LSV IDs that are shared by all LSV dictionaries
        in Data are evaluated. LSV IDs that are unique to or missing from
        any of the LSV Dictionaries are not returned by this function.

    Returns the numpy array and a list of LSV_IDs
    """
    if use_binary_index_info:
        if not use_binary_index_info == "closer" and not use_binary_index_info == "further":
            raise ValueError("Use_binary_index_info needs to be 'closer' or 'further' if provided")
    if use_binary_index_info == "closer":
        binary_index = 0
    if use_binary_index_info == "further":
        binary_index = 1
    comparison = "thing to make PyCharm happy"
    binary_index = "thing to make PyCharm happy"
    check_is_quick_import(data)
    condition_dict = dict()
    lsv_id_lists = list()
    for comparison in data.keys():
        lsv_dict = data[comparison]
        # cond_1_name = lsv_dict["condition_1_name"]
        # cond_2_name = lsv_dict["condition_2_name"]
        cond_1_name = get_cond_1_name(lsv_dict)
        cond_2_name = get_cond_2_name(lsv_dict)
        PSIs = list_psi(lsv_dict)
        cond_1_PSIs = PSIs[cond_1_name]
        cond_2_PSIs = PSIs[cond_2_name]
        condition_dict[cond_1_name] = cond_1_PSIs
        condition_dict[cond_2_name] = cond_2_PSIs
    union_of_lsv_ids = get_shared_lsv_ids(data)
    lsv_dict = dict()
    conditions = list(condition_dict.keys())
    conditions.sort()
    for lsv_id in list(union_of_lsv_ids):
        # if lsv_id == "ENSG00000173744:228395807-228395926:target":
        #     pdb.set_trace()
        list_of_PSIs = list()
        for condition in conditions:
            PSI = condition_dict[condition][lsv_id]
            if use_binary_index_info:
                binary_i = data[comparison][lsv_id]["binary_indices"][binary_index]
                PSI = PSI[binary_i]
            list_of_PSIs.append(PSI)
        lsv_dict[lsv_id] = np.array(list_of_PSIs).T
    if return_comparisons:
        return lsv_dict, conditions
    return lsv_dict


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
        lsvs = get_lsv_ids(lsv_dict)
        all_lsvs.extend(lsvs)
        if verbose:
            n_lsvs = len(lsvs)
            LOG.info("%s LSVs in %s" % (n_lsvs, comparison_name))
    return list(set(all_lsvs))


def get_shared_lsv_ids(data, bool=False):
    """
    Given a quick import format, return LSV IDs
        that are shared by all LSV DIctionaries.
    """
    check_is_quick_import(data)
    lsvids = list()
    for comparison in data.keys():
        lsvids.append(set(get_lsv_ids(data[comparison])))
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
                LOG.info("Analyzing " + LSV_dict_name + " ...")
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
            LOG.info("No junctions in the LSVs had a dPSI of " + str(Change) +
                     ", now searching for biggest dPSI...")

        # Keep trying to find LSVs starting at Change %,
        #  and going down bit by bit until a LSV is returned.
        for cutoff in calebs_xrange(Change, 0, -0.01):
            if int(cutoff * 100.0) % int(0.05 * 100.0) == 0:
                if Verbose:
                    LOG.info("trying", cutoff, "...")
            # Don't sum for the cutoff! (False) <- see function desc
            over_cutoff = get_sig_lsv_ids(LSV_dict, cutoff, Probability, False)
            if len(over_cutoff) > 0:
                if Verbose:
                    LOG.info("Max dPSI identifed as: ", cutoff)
                break

    subset_lsv_dict = lsv_dict_subset(LSV_dict, over_cutoff, True)
    return subset_lsv_dict


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
                    LOG.info(LSV_dict_name + " was empty.")
            else:
                new_dict[LSV_dict_name] = lsv_dict
        if len(new_dict) == 0:
            LOG.info("Warning! All LSV_dicts were empty...")
        return new_dict
    check_is_lsv_dict(data)
    lsv_ids = get_lsv_ids(data)
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
        ids = copy.copy(get_lsv_ids(Data[comparison]))
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
    return lsv_dict["meta_info"]["condition_1_name"]


def get_cond_2_name(lsv_dict):
    """
    :param lsv_dict: lsv_dict struc ...
    :return: str condition 2 name
    """
    check_is_lsv_dict(lsv_dict)
    return lsv_dict["meta_info"]["condition_2_name"]


def get_sample_id(lsv_dict):
    check_is_lsv_dict(lsv_dict)
    return lsv_dict["meta_info"]["sample_id"]


def get_prob_threshold(lsv_dict):
    check_is_lsv_dict(lsv_dict)
    return lsv_dict["meta_info"]["prob_thresh"]


def gen_comparison_name(LSV_dict, sep="_"):
    """
    Given LSV dictionary, return condition_1_name[sep]condition_2_name
    """
    check_is_lsv_dict(LSV_dict)
    if not isinstance(sep, str):
        raise ValueError("sep needs to be a string, not a: " + str(type(sep)))
    return get_cond_1_name(LSV_dict) + sep + get_cond_2_name(LSV_dict)


def get_abs_path(lsv_dict):
    """
    :param lsv_dict: lsv dict..
    :return: str abs file path to txt file
    """
    check_is_lsv_dict(lsv_dict)
    return lsv_dict["meta_info"]["abs_path"]


def get_base_name(lsv_dict, sep="_"):
    """
    Given LSV dictionary, return condition_1_name[sep]condition_2_name
    """
    check_is_lsv_dict(lsv_dict)
    if not isinstance(sep, str):
        raise ValueError("sep needs to be a string, not a: " + str(type(sep)))
    type = psi_or_deltapsi(lsv_dict)
    if type == "deltapsi":
        return get_cond_1_name(lsv_dict) + sep + get_cond_2_name(lsv_dict)
    else:
        return get_sample_id(lsv_dict)


def get_all_lsv_ids(data):
    if not check_is_quick_import(data,
                                 the_bool=True):
        if check_is_lsv_dict(data,
                             da_bool=True):
            return get_lsv_ids(data)
        else:
            raise RuntimeError("Expected a LSV dictionary or Quick Import...")
    comparisons = data.keys()
    all_ids = list()
    for comparison in comparisons:
        all_ids.extend(get_lsv_ids(data[comparison]))
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


def comparisons_quantifiable(comparisons, blank_dict, lsv_id):
    """
    Determine which dPSI comparisons were able to quantify the LSV.
        Comparisons must be SORTED.

    :param comparisons: list of comparisons
    :param blank_dict: returned from impute_missing_lsvs()
    :param lsv_id: check which comparisons could quantify this LSV ID
    :return: numpy array of Bools corresponding to the comparisons able to quantify the LSV
    """
    bools = list()
    for comp in comparisons:
        this_bool = lsv_id in blank_dict[comp]
        bools.append(not this_bool)
    return np.array(bools)


def impute_lsvs(lsvs,
                imputing_with=0):
    """
    Given a list of or a single LSV, return the LSV(s) with 0s for
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
            LOG.info("WARNING: YOUR MAJIQ RESULTS WILL BE OVERWRITTERN SINCE InPlace=True")
    check_is_quick_import(data)
    unique_ids = set(get_all_unique_lsv_ids(data, verbose=False))
    new_dict = dict()
    blanked_dict = dict()
    for comparison in list(data.keys()):
        if verbose:
            LOG.info("Filling in the LSV gaps for %s ... " % comparison)
        this_comps_lsvids = set(get_lsv_ids(data[comparison]))
        only_in_unique = unique_ids - this_comps_lsvids
        blanked_dict[comparison] = only_in_unique
        if len(only_in_unique) == 0:
            LOG.info("%s has all the LSVs already!" % comparison)
            if not in_place:
                if verbose:
                    LOG.info("Deep copying...")
                new_dict[comparison] = copy.deepcopy(data[comparison])
            else:
                new_dict[comparison] = data[comparison]
            continue
        only_in_unique_lsvs = get_lsvs_quickly(data, only_in_unique)
        if verbose:
            LOG.info("Imputing...")
        only_in_unique_lsvs_blanked = impute_lsvs(only_in_unique_lsvs, imputing_with=impute_with)
        if verbose:
            LOG.info("Finished imputing...")
        only_in_unique_lsvs_blanked_dict = {x["LSV ID"]: x for x in only_in_unique_lsvs_blanked}
        if not in_place:
            # Copy the dict, new object!
            if verbose:
                LOG.info("Deep copying...")
            new_dict[comparison] = copy.deepcopy(data[comparison])
            if verbose:
                LOG.info("Updating...")
            new_dict[comparison].update(only_in_unique_lsvs_blanked_dict)
        else:
            data[comparison].update(only_in_unique_lsvs_blanked_dict)
    if in_place:
        return blanked_dict
        # return Data
    return new_dict, blanked_dict


def import_dpsi_pandas(tsv_file, columns=None):
    """

    :param tsv_file: path to tsv file
    :param columns: if provided as a list, only import columns at [indices provided]
    :return: pandas dataframe
    """
    # first, get the header of the tsv file (there is a # sign we need to remove)
    with open(tsv_file, "r") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            header = line
            header.replace("#", "")
            header = header.split("\t")
            # if columns:
            #     header = [header[ii] for ii in columns]
            break
    if columns:
        if isinstance(columns, int):
            columns = [columns]
        elif not isinstance(columns, list):
            LOG.error("Expected columns to be in int or [list of ints], "
                      "instead it was: %s and looked like: %s" % type(columns), columns)
            exit(1)
    else:  # else columns was empty, so import everything
        columns = range(len(header))
    pa_dataframe = pa.read_csv(tsv_file,
                               sep="\t",
                               header=None,
                               skiprows=1,
                               names=header,
                               usecols=columns)
    return pa_dataframe


def file_to_list(filepath):
    """
    Given a file path, return list of each line rstriped
    :param filepath: str to file..
    :return: list
    """
    filelist = list()
    with open(filepath, "r") as handle:
        for line in handle:
            filelist.append(line.rstrip("\n\r"))
    return filelist
