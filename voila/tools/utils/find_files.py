import os
import subprocess
import platform
import fnmatch
import pdb

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


# TODO use find better

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


def get_voila_files(directory, pattern, get_base_names=False):
    """
    Recursive search directory for files matching pattern. Speed-optimized for searching
        a directory where majiq was run.
    :param directory: path where txt files are
    :param pattern: grep-capable pattern match for file
    :param get_base_names: True/False ... only change this to True if you really want to
    :return: list of full file paths
    """
    from voila.tools.utils import io_caleb
    if not (os.path.isdir(directory)):
        raise ValueError("Directory doesn't seem to exist. ~_~")
    base_names = list()
    voila_txt_files = list()
    # Ensure ./blah becomes /abs/path/blah
    directory = os.path.abspath(directory)
    # Find all voila dPSI text files in directory
    for root, subdirs, files in os.walk(directory):
        exclude = ["summaries",
                   "static",
                   "doc",
                   "lsvs",
                   "tmp"]
        # don't os.walk down these directories..
        subdirs[:] = [d for d in subdirs if d not in exclude]
        found = voila_txt_finder(root, pattern=pattern)
        if found:
            for found_file in found:
                if get_base_names:
                    base_names.append(io_caleb.get_base_names(found_file))
                voila_txt_files.append(os.path.join(root, found_file))
        else:
            continue
        break
    if get_base_names:
        return base_names, voila_txt_files
    return voila_txt_files


def voila_txt_finder(directory, pattern):
    files = find_files(directory, pattern)
    if len(files) > 0:
        return files
    else:
        return False


def find_files(path,
               pattern,
               recursive=True,
               must_be_files=True):
    """
    Given a path, return all files that match pattern.

    If Mac or Linux, this uses subprocess.Popen
    else, uses os.walk(), which is much slower
    :param path:
    :param pattern: grep style patterm matching
    :param recursive: search inside directories recursively?
    :param must_be_files: after finding matches, make sure they are files?
        If this is False, this function will also return directories that match
        the pattern...
    :return: [list of file paths]
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
        # Save found files to list (each result has a newline character
        file_list = stdout.split(os.linesep)
        # remove empy elements (probably last thing in the stdout is '\n' which is '' now)
        while '' in file_list:
            file_list.remove('')
    else:
        # Walk through directory
        for dName, sdName, fList in os.walk(path):
            for fileName in fList:
                if fnmatch.fnmatch(fileName, pattern):  # Match search string
                    file_list.append(os.path.join(dName, fileName))
            if not recursive:
                return file_list
    if must_be_files:
        to_remove = list()
        to_keep = list()
        for ii in range(len(file_list)):
            if not os.path.isfile(file_list[ii]):
                to_remove.append(ii)
            else:
                to_keep.append(ii)
        file_list = [file_list[i] for i in to_keep]
    return file_list


def check_is_unix():
    system = platform.system()
    if "Darwin" in system:
        return True
    elif "Linux" in system:
        return True
    else:
        return False