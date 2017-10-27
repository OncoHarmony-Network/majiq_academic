"""
Script used for the identification and discovery of potential cryptic exons in RNA-seq datasets
the program takes as input a list of paths to .tsv files from the Voila output and creates a
dictionary of ctrl LSVs that can be compared to future samples to search for cryptic exons. The
script was developed as an add-on to voila, voila tools.

Developed by Nicholas Page @ the Barash (BioCiphers) Lab, UPenn
Any problems or questions please email: npage98@gmail.com
"""

import os.path
import pickle

from voila.tools import Tool


class CtrlPsiDictionaryTool(Tool):
    help = 'Makes a dictionary of Ctrl PSI values for each junction provided to the program from a .tsv file. These ' \
           'values will be compared to dPSI values from new experiments to help identify potential cryptic exons.'

    def arguments(self):
        parser = self.get_parser()

        parser.add_argument('--input',
                            metavar='INPUT',
                            required=True,
                            help='File path to the input file list containing the paths to all of the .tsv '
                                 'files that will be used to build the dPSI dictionary. The default option '
                                 'assumes that the .tsv file contains a header. [Default: None]')

        parser.add_argument('--output',
                            metavar='OUTDIR',
                            default=os.getcwd(),
                            help='File path to the directory where the dictionary file should be '
                                 'created and saved. [Default: ./]')

        parser.add_argument('--has_header',
                            default='True',
                            help='Specify whether or not the .tsv files listed in the input filelist '
                                 'contain headers. [Default: True] otherwise, False.')

        parser.add_argument('--get_both_ePSI_values',
                            default='False',
                            help='Adds the ePSI values from both index 5 and 6 for each line in a .tsv '
                                 'file. Use this option if the .tsv file is from a comparison between '
                                 'two control (non disease state) tissues. [Default: False] otherwise, True.')
        return parser

    def run(self, args):
        makeDictionary(args)


def getMaxPSI(ctrl_PSIs):
    """
    Gets the maximum Ctrl PSI value for each junction in an LSV
    :param ctrl_PSIs: The list of Ctrl PSI values for each junction found in multiple datasets
                        that is already saved within the dictionary.
    :return: max_PSI -> The maximum Ctrl PSI value found at a certain junction.
    """
    ctrl_PSI_list = [float(x) for x in ctrl_PSIs.split(';')]
    max_PSI = max(ctrl_PSI_list)
    return max_PSI


def makeDictionary(args):
    """
    Makes a dictionary of the maximum Ctrl PSI values for each LSV provided by the user via a list of .tsv files.
    The dictionary is saved as a text file in the output directory specified by the user
    :param args: a list containing the input file list and the output directory provided by the
                   user that will be used to create and store the Ctrl PSI dictionary.
    :return: None
    """
    dict = {}

    filelist = open(args.input, 'r')

    # Adds each junction to the dictionary from every file in the file list
    for line in filelist:

        file = open(line.rstrip('\n'), 'r')

        # Removes the header from the input .tsv files if specified by the user.
        if args.has_header == 'True':
            file.readline()

        # Adds each junction to the dictionary from every LSV (line) in the file list
        for line in file:

            # Obtains the gene ID, the list of Ctrl PSI values and the list of junction
            # coordinates for each LSV
            LSV = line.split('\t')
            geneID = LSV[1]
            ctrl_PSIs = LSV[5].split(';')
            junction_coords = LSV[16].split(';')

            # Gets both ePSI values from the .tsv file if specified by the user. This option should be
            # used when the .tsv file contains the comparison between two control tissues.
            if args.get_both_ePSI_values == 'True':

                tmp_PSIs = LSV[6].split(';')
                tmp_list = []

                # Combines the ePSI values from junctions with the same coordinates into a list.
                for ctrl_PSI, tmp_PSI in zip(ctrl_PSIs, tmp_PSIs):
                    tmp_list.append(ctrl_PSI + ';' + tmp_PSI)

                ctrl_PSIs = tmp_list

            # Appends new junctions to an existing gene in the dictionary
            if geneID in dict.keys():

                # Appends previously unseen junctions to the dictionary based on the gene in which
                # the junction in found
                for junction_coord, ctrl_PSI in zip(junction_coords, ctrl_PSIs):

                    if junction_coord in dict[geneID].keys():

                        dict[geneID][junction_coord] = dict[geneID][junction_coord] + ';' + ctrl_PSI

                    else:

                        dict[geneID].update({junction_coord: ctrl_PSI})

            # Creates a new gene in the dictionary for previously unseen genes
            else:

                dictJunc = {}

                for junction_coord, ctrl_PSI in zip(junction_coords, ctrl_PSIs):
                    dictJunc.update({junction_coord: ctrl_PSI})

                dict.update({geneID: dictJunc})

    # Replaces the list of Ctrl PSI values for each junction in a dictionary
    # with the maximum Ctrl PSI value in that list
    for geneID in dict:
        for junc in dict[geneID]:
            dict[geneID][junc] = getMaxPSI(dict[geneID][junc])

    # Dumps the Ctrl PSI dictionary into a text file in the directory provided by the user
    output = open(args.output + 'Ctrl_PSI_Dictionary', 'w')
    pickle.dump(dict, output)
