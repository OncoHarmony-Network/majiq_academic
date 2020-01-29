from voila.config import ClassifyConfig
from voila import constants
from voila.voila_log import voila_log
from voila.exceptions import VoilaException, UnknownAnalysisType, UnsupportedAnalysisType
from voila.api import SpliceGraph, Matrix
from voila.classifier.as_types import Graph
from voila.classifier.tsv_writer import TsvWriter
from voila.classifier.training_writer import TrainingWriter
from math import ceil
import time
import os, sys
from multiprocessing import Manager, Pool
import glob
import traceback
from itertools import islice


class Classify:
    def __init__(self):
        """
        Factory class used to create the Classification for the specific analysis type.
        """
        config = ClassifyConfig()
        analysis_type = config.analysis_type

        voila_log().info(analysis_type + ' CLASSIFY')

        run_classifier()



def print_filter_summary():
    config = ClassifyConfig()
    jst = 64
    voila_log().info("╔═══════════════════════════════════════════════════════════════╗")

    if config.decomplexify_reads_threshold or config.decomplexify_deltapsi_threshold or config.decomplexify_psi_threshold:
        voila_log().info(("╠╡ Before Classification:").ljust(jst)[:jst] + "║")
        if config.decomplexify_psi_threshold:
            voila_log().info(("║     Dropping junctions with max(PSI) < %0.3f" % config.decomplexify_psi_threshold).ljust(jst)[:jst] + "║")
        if config.decomplexify_deltapsi_threshold:
            voila_log().info(("║     Dropping junctions with max(dPSI) < %0.3f" % config.decomplexify_deltapsi_threshold).ljust(jst)[:jst] + "║")
        if config.decomplexify_reads_threshold:
            voila_log().info(("║     Dropping junctions with max(number_of_reads) < %d" % config.decomplexify_reads_threshold).ljust(jst)[:jst] + "║")
    else:
        voila_log().info(("╠╡ Not dropping any junctions").ljust(jst)[:jst] + "║")

    if config.changing or config.non_changing:
        voila_log().info(("╠╡ After Classification:").ljust(jst)[:jst] + "║")
        if config.changing:
            voila_log().info(("║     Dropping modules if none(max(E(|dPSI|)) >= %.3f)" % config.non_changing_threshold).ljust(jst)[
                             :jst] + "║")
            voila_log().info(("║     Dropping modules if none(max(P(|dPSI|>=%.3f)) >= %.3f)" %
                             (config.changing_threshold, config.probability_changing_threshold)).ljust(jst)[:jst] + "║")
        elif config.non_changing:
            voila_log().info(("║     Dropping modules if any(max(E(|dPSI|)) <= %.3f)" % config.non_changing_threshold).ljust(jst)[
                             :jst] + "║")
            voila_log().info(("║     Dropping modules if any(max(P(|dPSI|<=%.3f)) >= %.3f)" %
                             (config.non_changing_threshold, config.probability_non_changing_threshold)).ljust(jst)[:jst] + "║")

    voila_log().info("╚═══════════════════════════════════════════════════════════════╝")

def classify_gene(args):

    gene_id, experiment_names, q = args
    config = ClassifyConfig()

    try:
        graph = Graph(gene_id, experiment_names)


        if not config.output_training_data:

            writer = TsvWriter(graph, gene_id)


            if config.putative_multi_gene_regions:
                writer.p_multi_gene_region()

            else:
                if 'events' in config.enabled_outputs:
                    writer.cassette()

                    writer.alt3prime()
                    writer.alt5prime()
                    writer.alt3and5prime()

                    writer.p_alt3prime()
                    writer.p_alt5prime()

                    writer.mutually_exclusive()
                    writer.alternative_intron()

                    writer.alternate_first_exon()
                    writer.alternate_last_exon()
                    writer.p_alternate_first_exon()
                    writer.p_alternate_last_exon()

                    writer.multi_exon_spanning()
                    writer.tandem_cassette()
                    writer.exitron()
                    writer.orphan_junction()

                if 'summary' in config.enabled_outputs:
                    writer.summary()

                if 'heatmap' in config.enabled_outputs:
                    writer.heatmap()
                if 'junctions' in config.enabled_outputs:
                    writer.junctions()
                if 'mpe' in config.enabled_outputs:
                    writer.mpe()

                if ClassifyConfig().keep_constitutive and 'events' in config.enabled_outputs:
                    writer.constitutive()

        else:
            writer = TrainingWriter(graph, gene_id)

            if 'matrices' in config.enabled_outputs:
                writer.adjacency_matrix()

            if 'exons' in config.enabled_outputs:
                writer.exons_tsv()

            if 'junctions' in config.enabled_outputs:
                writer.junctions_tsv()

            if 'paths' in config.enabled_outputs:
                writer.paths_tsv()


    except KeyboardInterrupt:
        raise

    except:
        if config.debug:
            print(traceback.format_exc())
        voila_log().warning("Some error processing gene %s , turn on --debug for more info" % gene_id)

    if q:
        q.put(None)

def run_classifier():

    config = ClassifyConfig()

    experiment_names = set()
    for voila_file in config.voila_files:
        with Matrix(voila_file) as m:
            for grp in m.experiment_names:
                for exp in grp:
                    if exp:
                        experiment_names.add(exp)


    if not config.gene_ids:
        with SpliceGraph(config.splice_graph_file) as sg:
            if config.debug_num_genes:
                gene_ids = list(g['id'] for g in islice(sg.genes(), config.debug_num_genes))
            else:
                gene_ids = list(g['id'] for g in sg.genes())
    else:
        gene_ids = config.gene_ids

    #gene_ids = gene_ids[:20]

    if not os.path.exists(config.directory):
        os.makedirs(config.directory)

    if os.listdir(config.directory):
        if config.overwrite:
            log_file_relpath = config.logger.split('/')[-1] if config.logger else None
            for _file in os.listdir(config.directory):
                if _file != log_file_relpath:
                    os.remove(os.path.join(config.directory, _file))
        else:
            voila_log().critical(
                "Files already present in %s; not running classifier "
                "(--overwrite to delete files in this directory and run anyway)" % config.directory)
            sys.exit(1)

    voila_log().info("Classifying %d gene(s)" % len(gene_ids))
    voila_log().info("Quantifications based on %d input file(s)" % len(config.voila_files))

    print_filter_summary()

    voila_log().info("Writing TSVs to %s" % os.path.abspath(config.directory))

    #total_genes = len(gene_ids)
    if not config.output_training_data:
        writerClass = TsvWriter
    else:
        writerClass = TrainingWriter
        TrainingWriter.delete_tsvs()

    writerClass.delete_tsvs()

    work_size = len(gene_ids)

    if config.debug:
        try:
            for i, gene_id in enumerate(gene_ids):
                # t1 = time.time()
                classify_gene((gene_id, experiment_names, None,))
                # t2 = time.time()
                # print(t2-t1)
                print('Processing Genes and Modules [%d/%d]\r' % (i, work_size), end="")
        except KeyboardInterrupt:
            print('                                                  \r', end="")


        print('                                                  \r', end="")

    else:
        manager = Manager()
        q = manager.Queue()

        p = Pool(config.nproc)



        # voila_index = p.map(self._heterogen_pool_add_index, zip(lsv_ids, range(work_size), repeat(work_size)))
        classifier_pool = p.map_async(classify_gene, ((x, experiment_names, q) for x in gene_ids),)

        # monitor loop
        while True:

            if classifier_pool.ready():
                break
            else:
                size = q.qsize()
                print('Processing Genes and Modules [%d/%d]\r' % (size, work_size), end="")
                time.sleep(2)

        print('                                                  \r', end="")
        res = classifier_pool.get()

    voila_log().info("Concatenating Results")


    writer = writerClass(None, None)
    writer.start_all_headers()

    for _tsv in writerClass.tsv_names():
        read_files = glob.glob(os.path.join(config.directory, _tsv) + ".*")
        with open(os.path.join(config.directory, _tsv), "rb") as outfile:
            headers = outfile.read()
        with open(os.path.join(config.directory, _tsv), "wb") as outfile:
            outfile.write(headers)
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
                os.remove(f)

    if config.output_training_data:
        # in this case we can't know each specific prefix ahead of time, we need to search for them
        # (group files starting with X name, then '.hdf5.' then some number
        groups = {}
        read_files = glob.glob(config.directory + '/*.hdf5.*')
        for _file in read_files:

            group, ext, proc = _file.split('.')
            dest = "%s.%s" % (group, ext)
            if not dest in groups:
                groups[dest] = []
            groups[dest].append(_file)

        for group in groups:
            writer.combine_hdf5s(group, groups[group])



    voila_log().info("Classification Complete")


# import multiprocessing
#
#
# def Writer(dest_filename, some_queue, some_stop_token):
#     with open(dest_filename, 'w') as dest_file:
#         while True:
#             line = some_queue.get()
#             if line == some_stop_token:
#                 return
#             dest_file.write(line)
#
#
# def the_job(some_queue):
#     for item in something:
#         result = process(item)
#         some_queue.put(result)
#
#
# if __name__ == "__main__":
#     queue = multiprocessing.Queue()
#
#     STOP_TOKEN = "STOP!!!"
#
#     writer_process = multiprocessing.Process(target=Writer, args=("output.txt", queue, STOP_TOKEN))
#     writer_process.start()
#
#     # Dispatch all the jobs
#
#     # Make sure the jobs are finished
#
#     queue.put(STOP_TOKEN)
#     writer_process.join()
#     # There, your file was written.