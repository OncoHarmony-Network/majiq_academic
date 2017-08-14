"""
Script used for the identification and discovery of potential cryptic junctions in RNA-seq datasets
the program takes as input a list of paths to .tsv files from the Voila output and outputs a
list of putative cryptic junctions in the standard voila output format.

Developed by Nicholas Page @ the Barash (BioCiphers) Lab, UPenn
Any problems or questions please email: npage98@gmail.com
"""
import os.path
import pickle

from voila.tools import Tool


class PutativeCryticJunctionsTool(Tool):
    help = 'This script takes a file list as input that contains the ' \
           'paths to the .tsv files provided by the user as input ' \
           'and finds potential cryptic exons in the ' \
           'dataset by searching for low PSI values in the control ' \
           'condition and removing annotated junctions for each LSV. ' \
           'The script outputs a list of potential cryptic exons in ' \
           'the standard .tsv format.  Note: Output from both ' \
           'ctrl_psi_dictionary and junction_dictionary are required to run this tool.'

    def arguments(self):
        parser = self.get_parser()

        parser.add_argument('--input',
                            metavar='INPUT',
                            default=os.getcwd(),
                            help='File path to the filelist containing the paths to the .tsv '
                                 'files to be used as input. [Default: ./]')

        parser.add_argument('--output',
                            metavar='OUTDIR',
                            default=os.getcwd(),
                            help='File path to where the output .tsv file should be saved. [Default: ./]')

        parser.add_argument('--ctrl_PSI_threshold',
                            default=0.05,
                            help='The PSI value below which and exon is considered not present in the '
                                 'control sample. [Default: 0.05]')

        parser.add_argument('--dPSI_threshold',
                            metavar='dPSI_THRESHOLD',
                            default=0.2,
                            help='The dPSI value which an LSV must exceed to be considered. [Default: 0.2]')

        parser.add_argument('--dPSI_confidence_threshold',
                            metavar='dPSI_CONFIDENCE_THRESHOLD',
                            default=0.95,
                            help='The confidence threshold for dPSI that an LSV must exceed to be considered. '
                                 '[Default: 0.95]')

        parser.add_argument('--ctrl_PSI_dictionary',
                            help='The filepath to a dictionary text file containing the control PSI values for '
                                 'junctions found in healthy tissues. [Default: None] REQUIRED.')

        parser.add_argument('--junction_dictionary',
                            help='The filepath to a s dictionary text file containing a list of all previously '
                                 'annotated splicing junctions based off a transcript database. [Default: None] '
                                 'REQUIRED.')

        return parser

    def run(self, args):
        putative_cryptic_junctions(args)


def maybe_cryptic(args, line, threshold):
    """
    Determines if the given LSV is cryptic based on the thresholds provided by the user.
    :param args:
    :param line: .tsv output representing an individual LSV
    :param threshold: An array of the threshold values used to determine if an LSV is cryptic
                        This contains the ctrl_PSI_threshold which is used to determine if the
                        LSV is present in the control sample (i.e. it should be a very small value),
                        the dPSI_threshold is used as a cutoff for the amount of change in an LSV
                        to be considered [Default: 0.2], the dPSI_confidence_threshold is the
                        probability the dPSI value reaches the threshold that was set while
                        running voila (i.e. this should be a very large value [Default: 0.95].
    :return: True if the LSV is cryptic and False otherwise
    """
    LSVdictFile = open(args.ctrl_PSI_dictionary, 'r')
    LSVdict = pickle.load(LSVdictFile)

    coordDictFile = open(args.junction_dictionary, 'r')
    coordDict = pickle.load(coordDictFile)

    # Extracts the threshold values from the list passed to the method
    ctrl_PSI_threshold = float(threshold[0])
    dPSI_threshold = float(threshold[1])
    dPSI_confidence_threshold = float(threshold[2])

    # Obtains the dPSI value, ctrl PSI value, confidence value, and de novo boolean for each junction
    # in an LSV and stores each of them in a list.
    LSV = line.split('\t')
    geneID = LSV[1]
    dPSI_list = [abs(float(x)) for x in LSV[3].split(";")]
    dPSI_confidence_list = [abs(float(x)) for x in LSV[4].split(";")]
    ctrl_PSI_list = [abs(float(x)) for x in LSV[5].split(";")]
    coords = [x for x in LSV[16].split(';')]

    # Iterates through each junction in the LSV provided to the maybeCryptic method.
    for dPSI, dPSI_confidence, ctrl_PSI, coord in \
            zip(dPSI_list, dPSI_confidence_list, ctrl_PSI_list, coords):

        # Checks if a junction may be cryptic based on the parameters specified by the user. It uses a try-except statement
        # because an exception is thrown if a coordinate that is checked for in the LSVdictionary is not there.
        try:

            # Checks if a junction may be cryptic based on the parameters specified by the user.
            if (dPSI >= dPSI_threshold) and (dPSI_confidence >= dPSI_confidence_threshold) \
                    and (ctrl_PSI <= ctrl_PSI_threshold) and (LSVdict[geneID][coord] <= ctrl_PSI_threshold) \
                    and coord not in coordDict[geneID]:
                return True

        # Checks if a junction may be cryptic if an exception is thrown when searching for the junction coordinates
        # in the LSVdictionary.
        except:

            if (dPSI >= dPSI_threshold) and (dPSI_confidence >= dPSI_confidence_threshold) \
                    and (ctrl_PSI <= ctrl_PSI_threshold) and coord not in coordDict[geneID]:
                return True

    return False


def putative_cryptic_junctions(args):
    """
    Handles the arguments provided by the user and runs the isCryptic methods
    :param args:
    :return: Outputs a .tsv file with all of the potential cryptic exons
    """
    # Creates a list of all of the paramaters to pass to the isCryptic method
    thresholds = [args.ctrl_PSI_threshold, args.dPSI_threshold, args.dPSI_confidence_threshold]

    filelist = open(args.input, 'r')

    # Iterates through each file in the filelist and searches for potential cryptic junctions.
    for line in filelist:

        file = open(line.rstrip('\n'), 'r')

        # Creates a new .tsv file in the output directory named cryptic_PREVIOUS_FILE_NAME
        output = open(args.output + '/cryptic_' + line.split('/')[-1], 'w')
        output.write(file.readline())

        # Checks if each line represnting an LSV in the .tsv file is cryptic
        for line in file:
            if maybe_cryptic(args, line, thresholds) == True:
                output.write(line)

        output.close()
