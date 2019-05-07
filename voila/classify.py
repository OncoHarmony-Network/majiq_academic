from voila.config import ClassifyConfig
from voila import constants
from voila.voila_log import voila_log
from voila.exceptions import VoilaException, UnknownAnalysisType, UnsupportedAnalysisType
from voila.api import SpliceGraph, Matrix
from voila.classifier.as_types import Graph
from voila.classifier.tsv_writer import TsvWriter

import os

class Classify:
    def __init__(self):
        """
        Factory class used to create the Classification for the specific analysis type.
        """
        config = ClassifyConfig()
        analysis_type = config.analysis_type

        voila_log().info(analysis_type + ' CLASSIFY')

        if analysis_type == constants.ANALYSIS_PSI:
            run_classifier()
        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            raise UnsupportedAnalysisType(analysis_type)
        elif analysis_type == constants.ANALYSIS_HETEROGEN:
            raise UnsupportedAnalysisType(analysis_type)
        else:
            raise UnknownAnalysisType(analysis_type)



def run_classifier():

    config = ClassifyConfig()

    if not config.gene_ids:
        with SpliceGraph(config.splice_graph_file) as sg:
           gene_ids = list(g['id'] for g in sg.genes())
    else:
        gene_ids = config.gene_ids


    if not os.path.exists(config.directory):
        os.makedirs(config.directory)

    voila_log().info("Classifying %d gene(s)" % len(gene_ids))
    voila_log().info("Writing TSVs to %s" % os.path.abspath(config.directory))


    TsvWriter.delete_tsvs()
    for gene_id in gene_ids:
        graph = Graph(gene_id)

        writer = TsvWriter(graph, gene_id)
        writer.cassette()

        writer.alt3prime()
        writer.alt5prime()
        # writer.alt3and5prime()

        writer.mutually_exclusive()
        writer.intron_retention()
        writer.multi_exon_skipping()
        writer.summary()

        writer.alternate_first_exon()
        writer.alternate_last_exon()
