import pdb

from voila.tools import Tool
from voila.tools.utils import io_caleb
from voila.utils.voila_log import voila_log


# Caleb Matthew Radens
# radlinsky@gmail.com

__author__ = 'cradens'


LOG = voila_log()


class ThisisLookup(Tool):
    help = 'Given a directory and Gene Name, or Gene ID, or LSV ID, prettily-print all LSVs. Alternatively, ' \
           'just return all the voila text files that match your pattern and/or comparisons.'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        mutually_excl_grp = parser.add_mutually_exclusive_group(required=True)
        mutually_excl_grp.add_argument('--lookup_val',
                            type=str,
                            help='Gene Names, Gene IDs, or LSV IDs to lookup. Comma separated, or a file'
                                 'with each line as a lookup value.')
        help_mes = 'Flag: only return file paths for the voila txt files found..'
        mutually_excl_grp.add_argument('--just_file_paths',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = 'Flag: use this if you just want to lookup result from one (random) text file'
        parser.add_argument('--just_one',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = 'Flag: Don\'t just print essential stats, print all the columns from LSV file.'
        parser.add_argument('--dont_abbreviate',
                            default=False,
                            action="store_true",
                            help=help_mes)
        help_mes = 'Which comparisons or samples to lookup ID in? Single space or comma separated please.'
        parser.add_argument('--names',
                            '--comparisons',
                            type=str,
                            help=help_mes)
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.0)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            type=float,
                            help=help_mes,
                            default=0.0)
        return parser

    def run(self, args):
        # parse the comparisons argument
        if args.names:
            if "," in args.names or " " in args.names:
                args.names.replace(" ", ",")
                to_lookup = args.names.split(",")
            else:
                to_lookup = [args.names]
            dont_remove_dups = False
        else:
            to_lookup = None
            dont_remove_dups=True
        if args.lookup_val:
            the_lookup_vals = args.lookup_val.replace(" ","").split(",")
        imported = io_caleb.quick_import(input=args.directory,
                                         cutoff_d_psi=args.dpsi_thresh,
                                         cutoff_prob=args.prob_dpsi_thresh,
                                         pattern=args.pattern,
                                         keep_ir=True,
                                         just_one=args.just_one,
                                         stop_at=the_lookup_vals,
                                         comparisons=to_lookup,
                                         just_file_paths=args.just_file_paths)
        io_caleb.check_is_ignant(imported, args.dpsi_thresh)
        if args.just_file_paths:
            print("\n".join(imported))
            return
        to_lookup = list(imported.keys())

        # coding for readability here...
        # default is True..
        abbreviated_bool = True
        if args.dont_abbreviate:
            abbreviated_bool = False
        lookup_everywhere(dictionary_lookup=imported,
                          name=the_lookup_vals,
                          just_one=args.just_one,
                          abbreviated=abbreviated_bool,
                          comparisons_lookup=to_lookup,
                          dont_rem_dup=dont_remove_dups)


def lookup_everywhere(dictionary_lookup,
                      name,
                      save_lsv_structure_lookup=True,
                      just_one=False,
                      print_bool=True,
                      abbreviated=True,
                      comparisons_lookup=False,
                      dont_rem_dup=False):
    """
    Given dictionary of LSV dictionaries, return or print LSV dictionaries
        within the gene if gene id or ENSG ID provided. Return the LSVs if LSV IDs
        provided.

        Args:
            dictionary_lookup: dict of LSV_dicitonaries (which is returned by quick_import)
            name: May be a gene name OR Ensembl Gene ID or LSV ID
            save_lsv_structure_lookup: if True, saves "condition_1/2_name"
            just_one: should just one LSV be printed?
            print_bool: should results be returned or printed?
            abbreviated: only print/return the useful data?
            comparisons_lookup: if provided, only look in these comparisons
            dont_remove_dup : When True, if comparisons_lookup has *duplicate# in it, then keep that
    """
    io_caleb.check_is_quick_import(dictionary_lookup)
    found_dicts = dict()
    found_data = False
    for lsv_dict_name in dictionary_lookup.keys():
        if comparisons_lookup:
            if not dont_rem_dup:
                if io_caleb.comp_without_dup(lsv_dict_name) not in \
                        [io_caleb.comp_without_dup(x) for x in comparisons_lookup]:
                    continue
        if isinstance(name, list):
            for na in name:
                found_data = lookup(dictionary_lookup[lsv_dict_name],
                            name=na,
                            printable=print_bool,
                            save_lsv_structure=save_lsv_structure_lookup,
                            not_found_error=False,
                            abbreviated=abbreviated)
                if found_data == "gene_not_found" or found_data == "lsv_id_not_found":
                    if print_bool:
                        print(na + " not found in " + lsv_dict_name + "\n")
                elif not found_data:
                    print("Nothing found...")
                    return
                else:
                    if lsv_dict_name in found_dicts:
                        if print_bool:
                            found_dicts[lsv_dict_name] += found_data
                        else:
                            found_dicts[lsv_dict_name].update(found_data)
                    else:
                        found_dicts[lsv_dict_name] = found_data
        else:
            found_data = lookup(dictionary_lookup[lsv_dict_name],
                                name=name,
                                printable=print_bool,
                                save_lsv_structure=save_lsv_structure_lookup,
                                not_found_error=False,
                                abbreviated=abbreviated)
            if found_data == "gene_not_found" or found_data == "lsv_id_not_found":
                if print_bool:
                    print(name + " not found in " + lsv_dict_name + "\n")
            elif not found_data:
                print("Nothing found...")
                return
            else:
                found_dicts[lsv_dict_name] = found_data
        if just_one:
            break
    if print_bool:
        print_lookup_everywhere(found_dicts)
    else:
        return found_dicts


def lookup(lsv_dictionary,
           name,
           save_lsv_structure=True,
           printable="print",
           not_found_error=True,
           abbreviated=True):
    """
        Given LSV dictionary and a Gene name OR Ensembl ID OR LSV ID, look up LSVs.

        Args:
            lsv_dictionary: lsv dict ...
            name: May be a gene name OR Ensembl Gene ID or LSV ID
            save_lsv_structure: if True, saves "condition_1/2_name"
            printable: True/False or 'print'
                True/False:
                Should the results be returned as a LSV dictionary (False)
                or as a printable string (True)?
                or simply printed? ('print')
            not_found_error: Leave this default...
            abbreviated: only print important info

        Returns subset of the dictionary pointing with LSVs in the gene.
    """
    if not isinstance(printable, bool) and printable != "print":
        raise ValueError("printable needs to be True/False or 'print'")
    if type(name) is not str:
        raise ValueError("Name needs to be a string.")
    if ":" in name:
        if "source" in name or "target" in name:
            lsv_id = True
        else:
            raise RuntimeError("%s is weird. Is a full LSV ID?" % name)
    else:
        lsv_id = False
    io_caleb.check_is_lsv_dict(lsv_dictionary)

    if lsv_id:
        possible_ids = io_caleb.get_lsv_ids(lsv_dictionary)
        if name in possible_ids:
            new_dict = io_caleb.lsv_dict_subset(lsv_dictionary, name, save_lsv_structure)
        else:
            return "lsv_id_not_found"
    else:  # gene name or Ensembl ID !
        if name[0:4] == "ENSG" or name[0:4] == "ENSM":
            ensembl_id = True
        else:
            ensembl_id = False
        matched_ids = list()
        all_ids = io_caleb.get_lsv_ids(lsv_dictionary)
        for lsvid in all_ids:
            if lsvid == "meta_info":
                continue
            if lsv_dictionary[lsvid]["Gene Name"].lower() == name.lower():
                matched_ids.append(lsvid)
                continue
            if ensembl_id:
                if lsv_dictionary[lsvid]["Gene ID"].lower() == name.lower():
                    matched_ids.append(lsvid)

        if len(matched_ids) == 0:
            if not_found_error:
                raise ValueError(name
                                 + " wasn't found in the provided LSV_dictionary.")
            else:
                return "gene_not_found"

        new_dict = io_caleb.lsv_dict_subset(lsv_dictionary, matched_ids, save_lsv_structure)
        if save_lsv_structure:
            correction = 1
        else:
            correction = 0
        if len(new_dict) - correction != len(matched_ids):
            raise BaseException("Something weird happened when extracting"
                                + " lsvs from the LSV dictionary...")
    if isinstance(printable, bool) and printable:
        return print_lsv(new_dict, print_bool=False, abbreviated=abbreviated)
    elif printable == "print":
        print_lsv(new_dict, print_bool=True, abbreviated=abbreviated)
        return
    return new_dict


def get_lsv(data, lsv_id, comparison=False):
    """
    Given LSV ID, and quick import, return any LSV from the
     data.

     Args:
         data : quick import
         lsv_id : lsv id...
        comparison: if comparison name (Str) given, use that to get the LSV.
    """
    io_caleb.check_is_quick_import(data)
    lsv_dicts = lookup_everywhere(data,
                                         lsv_id,
                                         save_lsv_structure_lookup=False,
                                         print_bool=False)
    if comparison:
        comparison = comparison
    else:
        comparison = lsv_dicts.keys()[0]
    lsv = lsv_dicts[comparison].values()[0]
    return lsv


def get_voila_link(vo_data,
                   vo_lsvid,
                   vo_comp=0,
                   vo_print=True,
                   fix_link=["data", "Volumes/data"]):
    """
    Given quick import of Data and LSV ID (or list of IDs),
        get the Voila links. Auto converts links to
        Samba-ready link.

        Data: quick import of data
        LSV_ID: a single LSV ID or a list of such
        Comp: int. Which comparison to use from Data to get link?
        Print: should the results be printed? IF False, returns them instead
        Fix_link: The voila link is broken. fix it by replacing:
        'data' with 'Volumes/data' in the link.
    """
    if isinstance(vo_lsvid, set):
        vo_lsvid = list(vo_lsvid)
    if not io_caleb.check_is_lsv_id(vo_lsvid, Bool=True):
        if isinstance(vo_lsvid, list):
            all_links = list()
            all_ids = list()
            for lsvid in vo_lsvid:
                link = get_voila_link(vo_data=vo_data,
                                      vo_comp=vo_comp,
                                      vo_lsvid=lsvid,
                                      fix_link=fix_link,
                                      vo_print=False)
                all_links.append(link)
                all_ids.append(lsvid)
            if vo_print:
                for lsv, link in zip(all_ids, all_links):
                    print("===")
                    print(lsv)
                    print(link)
                return
            else:
                return all_links, all_ids
    io_caleb.check_is_quick_import(vo_data)
    lsv_dict = vo_data[vo_data.keys()[vo_comp]]
    if vo_lsvid not in lsv_dict:
        raise RuntimeError("%s not found in %s" % (vo_lsvid, vo_data.keys()[vo_comp]))
    link = lsv_dict[vo_lsvid]["Voila link"]
    link = link.replace(fix_link[0], fix_link[1])
    if vo_print:
        print(vo_lsvid)
        print(link)
    else:
        return link


def print_lookup_everywhere(printable_stuff):
    """
        Printable_stuff returned by lookup_everywhere()
    """
    return_text = ""
    for lsv_dict_name in list(printable_stuff.keys()):
        lsv_dict = printable_stuff[lsv_dict_name]
        return_text += "==##=##== " + lsv_dict_name + " ==##==##==\n"
        return_text += lsv_dict
        return_text += "\n"
    print(return_text)


def print_lsv(lsv, print_bool=True, abbreviated=True):
    """ Given a LSV or dictionary of LSVs, return pretty-printable list
            of data or print the LSV, if Print=True
        Args:
            abbreviated: if True, only print a subset of the LSV that I think is useful.
            print_bool: Should dictionary be printed, or just returned?
            lsv: quick import, lsv dict, or lsv
    """
    printable = ""
    if io_caleb.check_is_lsv(lsv, True):
        helped_print = help_print_lsv(lsv, abbreviated)
        printable += helped_print
    elif io_caleb.check_is_lsv_dict(lsv, True):
        lsv_ids = io_caleb.get_lsv_ids(lsv)
        comparison_name = io_caleb.get_base_name(lsv)
        printable += "\n=====" + comparison_name + "=====\n"
        for lsv_id_p in lsv_ids:
            printable += print_lsv(lsv[lsv_id_p], print_bool=False, abbreviated=abbreviated)
            printable += "\n"
    elif io_caleb.check_is_quick_import(lsv, True):
        comparison_names = lsv.keys()
        for comparison_p in comparison_names:
            printable += print_lsv(lsv[comparison_p], print_bool=False, abbreviated=abbreviated)
    else:
        pdb.set_trace()
        raise ValueError("Can't print that, sorry.")

    if print_bool:
        print(printable)
    else:
        return printable


def help_print_lsv(lsv, abbreviated=True):
    """Given LSV, return a tab-delim string of all the data.

        This function is exclusively for print_lsv.
    """
    type_lsv = io_caleb.psi_or_deltapsi(lsv)
    if type_lsv == "deltapsi":
        the_header = ["Gene Name",
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
                       "De Novo Junctions",
                       "chr",
                       "strand",
                       "Junctions coords",
                       "Exons coords",
                       "Exons Alternative Start",
                       "Exons Alternative End",
                       "IR coords"]
        the_header[4] = io_caleb.get_name_of_prob_key(lsv)
        the_header[5], the_header[6] = io_caleb.get_name_of_psi_keys(lsv)
        headers_to_keep = [0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 16, 17]
    else:  # it is a psi file
        the_header = ["Gene Name",
                       "Gene ID",
                       "LSV ID",
                       "E(PSI) per LSV junction",
                       "Var(E(PSI)) per LSV junction",
                       "LSV Type",
                       "A5SS",
                       "A3SS",
                       "ES",
                       "Num. Junctions",
                       "Num. Exons",
                       "De Novo Junctions",
                       "chr",
                       "strand",
                       "Junctions coords",
                       "Exons coords",
                       "Exons Alternative Start",
                       "Exons Alternative End",
                       "IR coords"]
        headers_to_keep = [0, 1, 2, 3, 4, 5, 12, 13, 14, 15]
    stringed_lsv = ""
    if abbreviated:
        for item_i in headers_to_keep:
            item = the_header[item_i]
            stringed_lsv += item + "\t" + str(lsv[item]) + "\n"
    else:
        header_i_copy = list(range(0, len(the_header), 1))
        headers_to_keep = [0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 16, 17]
        all_headers = list()
        all_headers.extend(lsv.keys())
        for header_i in headers_to_keep:
            header_i_copy.remove(header_i)
            all_headers.remove(the_header[header_i])
            stringed_lsv += the_header[header_i] + "\t" + str(lsv[the_header[header_i]]) + "\n"
        for remainders in header_i_copy:
            stringed_lsv += the_header[remainders] + "\t" + str(lsv[the_header[remainders]]) + "\n"
            all_headers.remove(the_header[remainders])
        for custom_extra in all_headers:
            stringed_lsv += custom_extra + "\t" + str(lsv[custom_extra]) + "\n"
    return stringed_lsv
