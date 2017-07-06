from voila.tools import Tool
from voila.tools import io_voila_caleb
import os
import subprocess
import platform
import fnmatch
import argparse
import pdb

class ThisisLookup(Tool):
    help = 'Given a directory and Gene Name, or Gene ID, or LSV ID, pretty-print all LSVs'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        parser.add_argument('lookup_val',
                            type=str,
                            help='Gene Name, Gene ID, or LSV ID to lookup')
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
        help_mes = 'Which comparison(s) too lookup ID in? Possible formats:\n' \
                   '    --comparisons name_of_single_ID\n' \
                   '    --comparisons name1 name2 name3\n' \
                   '    --comparisons name1,name2,name3'
        parser.add_argument('--comparisons',
                            type=str,
                            help=help_mes)

        return parser

    def run(self, args):
        imported = io_voila_caleb.quick_import(dir=args.directory,
                                               cutoff_d_psi=0,
                                               cutoff_prob=0,
                                               pattern=args.pattern,
                                               keep_ir=True,
                                               just_one=args.just_one,
                                               stop_at=args.lookup_val)
        if args.comparisons:
            if "," in args.comparisons or " " in args.comparisons:
                args.comparisons.replace(" ", ",")
                to_lookup = args.comparisons.split(",")
            else:
                to_lookup = list(args.comparisons)
        else:
            to_lookup = list(imported.keys())
        avail_comps = list(imported.keys())
        for tolook in to_lookup:
            if tolook not in avail_comps:
                raise ValueError("%s is not a valid comparison name..." % tolook)
        # coding for readability here...
        # default is True..
        abbreviated_bool = True
        if args.dont_abbreviate:
            abbreviated_bool = False
        lookup_everywhere(dictionary_lookup=imported,
                          name=args.lookup_val,
                          just_one=args.just_one,
                          abbreviated=abbreviated_bool,
                          comparisons_lookup=to_lookup)


def lookup_everywhere(dictionary_lookup,
                      name,
                      save_lsv_structure_lookup=True,
                      just_one=False,
                      print_bool=True,
                      abbreviated=True,
                      comparisons_lookup=False):
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
    """
    io_voila_caleb.check_is_quick_import(dictionary_lookup)
    # if not isinstance(Dictionary,dict) or not isinstance(Dictionary[Dictionary.keys()[0]],dict):
    #     raise ValueError("Please provide a dictionary of dictionaries.")
    # if Dictionary[Dictionary.keys()[0]].has_key("Gene ID"):
    #     raise ValueError("Please provide a dictionary of LSV dictionaries,"
    #                         " not a single LSV dictionary. Or, use lookup().")
    # if not Dictionary[Dictionary.keys()[0]].has_key("condition_1_name"):
    #     raise ValueError("Please provide a dictionary of LSV dictionaries.")

    found_dicts = dict()
    found_data = False
    for LSV_dict_name in dictionary_lookup.keys():
        if comparisons_lookup:
            if LSV_dict_name not in comparisons_lookup:
                continue
        found_data = lookup(dictionary_lookup[LSV_dict_name],
                            name=name,
                            printable=print_bool,
                            save_lsv_structure=save_lsv_structure_lookup,
                            not_found_error=False,
                            abbreviated=abbreviated)
        if found_data == "gene_not_found" or found_data == "lsv_id_not_found":
            if print_bool:
                print(name + " not found in " + LSV_dict_name + "\n")
        elif not found_data:
            print("Nothing found...")
            return
        else:
            found_dicts[LSV_dict_name] = found_data
        if just_one:
            break
    if print_bool:
        _print_lookup_everywhere(found_dicts)
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
    io_voila_caleb.check_is_lsv_dict(lsv_dictionary)

    if lsv_id:
        possible_ids = io_voila_caleb.get_LSV_IDs(lsv_dictionary)
        if name in possible_ids:
            new_dict = io_voila_caleb.lsv_dict_subset(lsv_dictionary, name, save_lsv_structure)
        else:
            return "lsv_id_not_found"
    else:  # gene name or Ensembl ID !
        if name[0:4] == "ENSG" or name[0:4] == "ENSM":
            ensembl_id = True
        else:
            ensembl_id = False
        matched_ids = list()
        all_ids = io_voila_caleb.get_LSV_IDs(lsv_dictionary)
        for lsvid in all_ids:
            if lsvid == "condition_1_name" or lsvid == "condition_2_name":
                continue
            if lsv_dictionary[lsvid]["Gene Name"] == name:
                matched_ids.append(lsvid)
                continue
            if ensembl_id:
                if lsv_dictionary[lsvid]["Gene ID"] == name:
                    matched_ids.append(lsvid)

        if len(matched_ids) == 0:
            if not_found_error:
                raise ValueError(name
                                 + " wasn't found in the provided LSV_dictionary.")
            else:
                return "gene_not_found"

        new_dict = io_voila_caleb.lsv_dict_subset(lsv_dictionary, matched_ids, save_lsv_structure)
        if save_lsv_structure:
            correction = 2
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
    io_voila_caleb.check_is_quick_import(data)
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
    if not io_voila_caleb.check_is_lsv_id(vo_lsvid, Bool=True):
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
    io_voila_caleb.check_is_quick_import(vo_data)
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


def _print_lookup_everywhere(Printable_stuff):
    """
        Printable_stuff returned by lookup_everywhere()
    """
    return_text = ""
    for LSV_dict_name in Printable_stuff.keys():
        LSV_dict = Printable_stuff[LSV_dict_name]
        return_text += "==##=##== " + LSV_dict_name + " ==##==##==\n"
        return_text += LSV_dict
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
    if io_voila_caleb.check_is_lsv(lsv, True):
        helped_print = _help_print_lsv(lsv, abbreviated)
        printable += helped_print
    elif io_voila_caleb.check_is_lsv_dict(lsv, True):
        lsv_ids = io_voila_caleb.get_LSV_IDs(lsv)
        comparison_name = io_voila_caleb.get_comparison_name(lsv)
        printable += "\n=====" + comparison_name + "=====\n"
        for lsv_id_p in lsv_ids:
            printable += print_lsv(lsv[lsv_id_p], print_bool=False, abbreviated=abbreviated)
            printable += "\n"
    elif io_voila_caleb.check_is_quick_import(lsv, True):
        comparison_names = lsv.keys()
        for comparison_p in comparison_names:
            printable += print_lsv(lsv[comparison_p], print_bool=False, abbreviated=abbreviated)
    else:
        raise ValueError("Can't print that, sorry.")

    if print_bool:
        print(printable)
    else:
        return printable


def _help_print_lsv(lsv, abbreviated=True):
    """Given LSV, return a tab-delim string of all the data.

        This function is exclusively for print_lsv.

        Depends on keys defined in DPSI_HEADER.
    """
    io_voila_caleb.check_is_lsv(lsv)
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
                   "De Novo Junctions",
                   "chr",
                   "strand",
                   "Junctions coords",
                   "Exons coords",
                   "Exons Alternative Start",
                   "Exons Alternative End",
                   "IR coords"]
    dpsi_header[4] = io_voila_caleb.get_name_of_prob_key(lsv)
    dpsi_header[5], dpsi_header[6] = io_voila_caleb.get_name_of_psi_keys(lsv)
    stringed_lsv = ""
    if abbreviated:
        headers_to_keep = [0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 16, 17]
        for item_i in headers_to_keep:
            item = dpsi_header[item_i]
            stringed_lsv += item + "\t" + str(lsv[item]) + "\n"
    else:
        header_i_copy = list(range(0, len(dpsi_header), 1))
        headers_to_keep = [0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 16, 17]
        all_headers = list()
        all_headers.extend(lsv.keys())
        for header_i in headers_to_keep:
            header_i_copy.remove(header_i)
            all_headers.remove(dpsi_header[header_i])
            stringed_lsv += dpsi_header[header_i] + "\t" + str(lsv[dpsi_header[header_i]]) + "\n"
        for remainders in header_i_copy:
            stringed_lsv += dpsi_header[remainders] + "\t" + str(lsv[dpsi_header[remainders]]) + "\n"
            all_headers.remove(dpsi_header[remainders])
        for custom_extra in all_headers:
            stringed_lsv += custom_extra + "\t" + str(lsv[custom_extra]) + "\n"

    return stringed_lsv