from voila.tools import Tool
import os
import subprocess
import platform
import fnmatch
import pdb


class ThisisFindVOilaTexts(Tool):
    help = 'Given a directory, return all voila files inside it, recursively'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        help_mes = "Optional pattern matching to identify the voila text files\n" \
                   "Default for voila txt file: *tsv\n" \
                   "Default for deltapsi_voila: *.deltapsi.voila\n" \
                   "Default for deltapsi_prior: *.priomatrix.pkl"
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = 'Optional flag: return comparison names ?'
        parser.add_argument('--return-names',
                            action='store_true',
                            help=help_mes)
        help_mes = "Voila file type:"
        parser.add_argument('--file_type',
                            choices={"tsv", "deltapsi_voila", "prior_matrix"},
                            default="tsv",
                            type=str,
                            help=help_mes)

        return parser

    def run(self, args):
        print(find_voila_files(directory=args.directory,
                               pattern=args.pattern,
                               file_type=args.file_type,
                               get_comp_names=args.return_names))


def find_voila_files(directory, pattern, file_type, get_comp_names=False):
    """
    Find and return all "tsv", "deltapsi_voila", or "prior_matrix" files.
    :param directory:
    :param pattern:
    :param file_type: "tsv", "deltapsi_voila",or "prior_matrix"
    :param get_comp_names: True/False ... only change this to True if tsv and if you want
    :return:
    """
    if file_type not in ["tsv", "deltapsi_voila", "prior_matrix"]:
        raise ValueError("%s is not a valid file_type" % file_type)
    if file_type == "tsv":
        files = find_voila_files(directory=directory,
                                 pattern=pattern,
                                 get_comp_names=get_comp_names)
        return files
    if pattern == "*tsv":
        if file_type == "deltapsi_voila":
            pattern_match = "*.deltapsi.voila"
            files = find_voila_files(directory=directory,
                                     pattern=pattern_match,
                                     get_comp_names=False)
        else:  # file_type == "prior_matrix":
            pattern_match = "*.priomatrix.pkl"
            files = find_voila_files(directory=directory,
                                     pattern=pattern_match,
                                     get_comp_names=False)
        return files
    else:
        return find_voila_files(directory=directory,
                                pattern=pattern,
                                get_comp_names=False)


def get_voila_files(directory, pattern, get_comp_names=False):
    """
    Recursive search directory for files matching pattern. Speed-optimized for searching
        a directory where majiq was run.
    :param directory: path where txt files are
    :param pattern: grep-capable pattern match for file
    :param get_comp_names: True/False ... only change this to True if tsv and if you want
    :return: list of full file paths
    """
    if not (os.path.isdir(directory)):
        raise ValueError("Directory doesn't seem to exist. ~_~")
    dpsi_comparison_name = list()
    dpsi_files = list()
    # Find all voila dPSI text files in directory
    for root, subdirs, files in os.walk(directory):
        if os.path.basename(root) == "summaries" or os.path.basename(root) == "static":
            continue  # this folder contains all the html files.. skip over this dir
        elif os.path.basename(root) == "doc":
            if os.path.basename(os.path.dirname(root)) == "static":
                continue  # lots of docs in this dir, don't need to look through em..
        elif os.path.basename(root) == "lsvs":
            if os.path.basename(os.path.dirname(root)) == "doc":
                continue  # lots of .gtfs in this dir, don't need to look through em..
        # ID files that contain the Pattern
        found = _deltpsitextfinder(root, pattern=pattern)
        if found:
            for found_file in found:
                if get_comp_names:
                    dpsi_comparison_name.append(_get_deltapsi_txt_file_comparison(found_file))
                dpsi_files.append(os.path.join(root, found_file))
        else:
            continue
        break
    if get_comp_names:
        # dpsi_comparison_name = [x.decode('ascii') for x in dpsi_comparison_name]
        return dpsi_comparison_name, dpsi_files
    return dpsi_files


def _get_deltapsi_txt_file_comparison(filename):
    if os.path.exists(filename):
        basename = os.path.basename(filename)
    else:
        basename = filename
    split_file_name = basename.split(".")
    comparison_name = split_file_name[0]
    return comparison_name


def _deltpsitextfinder(directory, pattern):
    files = find_files(directory, pattern)
    if len(files) > 0:
        return files
    else:
        return False


def find_files(path, pattern, recursive=True):
    """
    Given a path, return all files that match pattern.

        If Mac or Linux, this uses subprocess.Popen
        else, uses os.walk(), which is much slower
    """
    if path:
        if not os.path.exists(path):
            raise IOError(path + " doesn't exist.")
    else:
        path = os.getcwd()
    if not pattern or not isinstance(path, str):
        raise ValueError("Please provide a valid str Pattern to match for file searching.")
    if not isinstance(recursive, bool):
        raise ValueError("Recursive needs to be True or False")
    file_list = []
    if check_is_unix():
        if recursive:
            recursive = ""
        else:
            recursive = " -maxdepth 1"
        find_cmd_line = 'find ' + path + recursive + ' -name "' + pattern + '"'
        out = subprocess.Popen(find_cmd_line, shell=True, stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Get standard out and error
        (stdout, stderr) = out.communicate()
        # need to decode bytes .. tested on Mac.. I bet this is going to break someday
        stdout = stdout.decode('ascii')
        # Save found files to list
        file_list = stdout.split()
    else:
        # Walk through directory
        for dName, sdName, fList in os.walk(path):
            for fileName in fList:
                if fnmatch.fnmatch(fileName, pattern):  # Match search string
                    file_list.append(os.path.join(dName, fileName))
            if not recursive:
                return file_list
    return file_list


def check_is_unix():
    system = platform.system()
    if "Darwin" in system:
        return True
    elif "Linux" in system:
        return True
    else:
        return False
