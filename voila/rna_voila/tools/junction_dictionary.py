"""
Script used for the identification and discovery of potential cryptic junctions in RNA-seq datasets
the program takes as input a .gff3 annotation database file and converts it into a dictionary
of all the annotated junction coordinates in the databases where they can be accessed by the key
which is the geneID of the gene that the junction belongs to.

Developed by Nicholas Page @ the Barash (BioCiphers) Lab, UPenn
Any problems or questions please email: npage98@gmail.com
"""
import os.path
import pickle

import gffutils

from rna_voila.tools import Tool


class JunctionDictionaryTool(Tool):
    help = 'Makes a dictionary of all of the junctions in an .gff3 annotation database indexed by the gene that they ' \
           'belong to. This annotation database will be used to check if a certain junction may belong to a cryptic ' \
           'junction.'

    def arguments(self):
        parser = self.get_parser()

        parser.add_argument('--input',
                            metavar='INPUT',
                            required=True,
                            help='File path to the input .gff3 file annotation database. [Default: None]')

        parser.add_argument('--output',
                            metavar='OUTDIR',
                            default=os.getcwd(),
                            help='File path to the directory where the dictionary file should be '
                                 'created and saved. [Default: ./]')

        return parser

    def run(self, args):
        junction_dictionary(args)


def junction_dictionary(args):
    """
    Creates a dictionary of all of the junctions in an .gff3 annotation database where their key is the geneID of the
    gene the junction belongs to.
    :param args: A list containg the path to the input .gff3 file that will be used to build
                   the junction database and the directory to which the output file containing the
                   dictionary should be written to.
    :return:
    """

    junction_dict = {}
    id_dict = {}

    annotation_file = open(args.input, 'r')

    # Iterates through every line in the .gff3 annotation database
    for line in annotation_file:

        # Skips over the lines in the .gff3 files that contain comments.
        if '#' in line: continue

        # Creates a gffutils feature object from each line in the .gff3 file.
        feature = gffutils.feature.feature_from_line(line)

        # Extracts information from each entry in the annotation database concerning the descriptive
        # annotation of each entry, the coordinates of that entry, and the parent structures of
        # each annotated entry.
        annotation = line.split('\t')[2]
        coord1 = line.split('\t')[3]
        coord2 = line.split('\t')[4]
        coords = '-'.join([coord1, coord2])

        # Adds each new annotated gene to a dictionary containing all of the geneIDs that is
        # accessible by the key corresponding to the gene name.
        if annotation == 'gene':
            id = feature.attributes['ID'][0]
            name = feature.attributes['Name'][0]

            # Creates a 'Nonetype' entry in the junction dictionary as a placeholder for each new gene detected.
            id_dict.update({name: id})
            junction_dict.update({id: None})

        # Adds each new annotated mRNA or transcript to the junction dictionary, accessible by a key
        # corresponding to the transcript ID.
        if annotation == 'mRNA' or annotation == 'transcript':

            # Converts the gene name in the attributes list to the gene ID.
            geneName = feature.attributes['parent_name'][0]
            geneID = id_dict.get(geneName)

            transcriptID = feature.attributes['ID'][0]

            # Adds the transcript to the junction dictionary if no previous transcripts
            # have been encountered from that gene.
            if junction_dict[geneID] == None:

                junction_dict[geneID] = {transcriptID: None}

            # Updates the junction dictionary with each additional transcript
            # found for each gene.
            else:

                junction_dict[geneID].update({transcriptID: None})

        # Adds the coordinates of each annotated junction to the junction dictionary.
        if annotation == 'exon':

            # Converts the gene name in the attributes list to the gene ID.
            geneName = feature.attributes['gene_name'][0]
            geneID = id_dict.get(geneName)

            transcriptID = feature.attributes['Parent'][0]

            # Adds an exon to a transcript dictionary for the given gene if an exon
            # for that gene has not been previously encountered.
            if junction_dict[geneID][transcriptID] == None:

                junction_dict[geneID][transcriptID] = [coords]

            # Adds an exon to the transcript dictionary for a given gene if an exon
            # for that gene has been previously encountered.
            else:

                junction_dict[geneID][transcriptID].append(coords)

    # Replaces the nested dictionaries for each gene with the list of all possible
    # annotated junctions for that gene.
    for geneID_key in junction_dict.keys():

        all_junctions = []

        # Iterates through each transcript that exists for a given gene.
        for transcriptID_key in junction_dict[geneID_key].keys():

            exons = junction_dict[geneID_key][transcriptID_key]
            all_coords = []

            # Used because it is possible to have a gene with no transcipt
            try:

                for exon in exons:
                    coord1 = exon.split('-')[0]
                    all_coords.append(coord1)
                    coord2 = exon.split('-')[1]
                    all_coords.append(coord2)

            except:

                pass

            # Sorts all of the exon coordinates for a given gene.
            all_coords = [int(i) for i in all_coords]
            all_coords.sort()
            all_coords = [str(i) for i in all_coords]

            # Converts the exon coordinate list provided into pairs of coordinates
            # representing the junctions between those exons.
            if len(all_coords) > 2:

                coord_num = len(all_coords)

                i = 1
                while i < coord_num - 2:
                    junc = []
                    junc.append(all_coords[i])
                    junc.append(all_coords[i + 1])
                    all_junctions.append('-'.join(junc))
                    i += 2

        # Stores all possible junctions for each gene in the dictionary.
        junction_dict[geneID_key] = all_junctions

    # Dumps the junction database as a text file in the directory provided by the user.
    output = open(args.output + '/Junction_Dictionary', 'w')
    pickle.dump(junction_dict, output)
