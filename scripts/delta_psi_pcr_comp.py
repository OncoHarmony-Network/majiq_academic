# from __future__ import division
from collections import defaultdict
import pickle
import numpy
from scipy.stats import pearsonr
import argparse
import ast
import os
import matplotlib.pyplot as pyplot
from scripts.utils import _save_or_show

__author__ = 'abarrera'


def plot_rtpcr_majiq_boxplot(rt_pcr, majiq, plotpath):
    #figure out how many groups of events exist

    fig = pyplot.plot(figsize=[6, 6], dpi=300)
    pyplot.title("Delta PSI comparison: RT-PCR Vs MAJIQ (N=%d)" % len(rt_pcr))

    print majiq, rt_pcr
    fit = numpy.polyfit(majiq, rt_pcr, 1)
    fit_fn = numpy.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    pyplot.plot(majiq, fit_fn(majiq), '--k')

    diagonal = numpy.linspace(-1, 1, num=len(rt_pcr))
    pyplot.plot(diagonal, diagonal, '--', color="#cccccc")
    pyplot.ylim([-1, 1])
    pyplot.plot(majiq, rt_pcr, '.', color='r')

    pyplot.xlabel('MAJIQ')
    pyplot.ylabel('RT-PCR')
    _save_or_show(plotpath, "delta_psi_comparison_rtpcr_majiq")


def barchart_expected(expected_psis, plotpath, mfile):
    pyplot.figure()
    pyplot.hist(expected_psis, bins=40)
    name = os.path.basename(mfile)
    _save_or_show(plotpath, name +'_expected_dist')


def collapse_matrix(matrix):
    "Collapse the diagonals probabilities in 1-D and return them"
    collapsed = []
    matrix_corner = matrix.shape[0]
    for i in xrange(-matrix_corner, matrix_corner):
        collapsed.append(numpy.diagonal(matrix, offset=i).sum())
    return numpy.array(collapsed)


def delta_expected_psi(bins):
    step = 2.0 / bins.size
    projection_prod = bins * numpy.arange(-1+step/2, 1, step)
    # if numpy.any(numpy.apply_along_axis(numpy.isnan, 0, numpy.array(projection_prod))):
    print numpy.sum(projection_prod)

    return numpy.sum(projection_prod)


def avg_expected_delta_psi(bins_list):
    return numpy.mean([delta_expected_psi(collapse_matrix(bins)) for bins in bins_list])


def coverage_from_file(file, result_dict, cov_suffix):
    with open(file) as majiq_file:
        majiq_builder = pickle.load(majiq_file)
        for lsv in majiq_builder[1]:
            for i, junc in enumerate(lsv.junction_list):
                result_dict[lsv.id+"#"+str(i)].append(junc.nnz)
    # Save results
    pickle.dump(result_dict, open(file+cov_suffix, 'w'))
    print "Coverage saved in: %s" % file+cov_suffix


def load_coverage(cov_file_list, cov_suffix):
    coverage_list = []
    try:
        for file in cov_file_list:
            coverage_list.append(pickle.load(open(file+cov_suffix)))
    except IOError:
        print "[INFO]::Coverage not found for %s..." % file
        return []

    return coverage_list

def get_coverage(majiq_files_or_dir, save_coverage=True):
    rest_file_list = []
    stim_file_list = []
    COV_SUFFIX = ".coverage"
    coverage_lists = []

    if os.path.isdir(majiq_files_or_dir):
        for file in os.listdir(majiq_files_or_dir):
            if file.endswith(".majiq") and "toJuan.resting" in file:
                rest_file_list.append(majiq_files_or_dir+'/'+file)
            if file.endswith(".majiq") and "toJuan.stim" in file:
                stim_file_list.append(majiq_files_or_dir+'/'+file)
    else:
        for file in majiq_files_or_dir:
            if "toJuan.resting" in file:
                rest_file_list.append(file)
            if "toJuan.stim" in file:
                stim_file_list.append(file)

    coverage_lists.append(load_coverage(rest_file_list, COV_SUFFIX))
    coverage_lists.append(load_coverage(stim_file_list, COV_SUFFIX))

    if len(coverage_lists[0]) and len(coverage_lists[1]):
        return coverage_lists

    from multiprocessing import Pool
    pool = Pool()
    rest_cov_dict_list = []
    stim_cov_dict_list = []
    jobs_list = []

    for rest_file in rest_file_list:
        rest_cov_dict = defaultdict(list)
        jobs_list.append(pool.apply_async(coverage_from_file, args=[rest_file, rest_cov_dict, COV_SUFFIX]))

    for stim_file in stim_file_list:
        stim_cov_dict = defaultdict(list)
        jobs_list.append(pool.apply_async(coverage_from_file, args=[stim_file, stim_cov_dict, COV_SUFFIX]))

    pool.close()
    pool.join()


    coverage_lists.append(load_coverage(rest_file_list, COV_SUFFIX))
    coverage_lists.append(load_coverage(stim_file_list, COV_SUFFIX))
    return coverage_lists


def main():
    parser = argparse.ArgumentParser(description="Compare Delta PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--majiq-deltas", required=True, dest='majiq_deltas', nargs='+', help='MAJIQ Delta PSI predictions for resting and stimuli RNA-Seq data.')
    parser.add_argument("--names-map-file", required=True, dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')
    parser.add_argument("--builder-files", required=True, dest='majiq_builder_files', help='File containing MAJIQ builder output file.')
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    ## Analyze coverage, saving it for later
    coverage_list = get_coverage(args.majiq_builder_files)

    # Read name mapping file
    names_pcr2majiq_dict = {}  # key RT-PCR, value MAJIQ
    names_majiq2pcr_dict = {}
    names_junc_majiq = defaultdict(list)
    gene_names_counts = defaultdict(lambda: 0)
    with open(args.names_map_file) as names_map:
        for name_map in names_map:
            # MAJIQ
            mapped_name = ast.literal_eval(name_map)
            [majiq_name, junc_num] = mapped_name[1].split('#')
            names_junc_majiq[majiq_name].append(int(junc_num))

            # RT-PCR
            if ':' in mapped_name[0]:
                names_pcr2majiq_dict[mapped_name[0].split(':')[1]] = mapped_name[1]
                names_majiq2pcr_dict[mapped_name[1]] = mapped_name[0].split(':')[1]
                gene_names_counts[mapped_name[0].split(':')[1]] += 1
            else:
                names_pcr2majiq_dict[mapped_name[0]] = mapped_name[1]
                names_majiq2pcr_dict[mapped_name[1]] = mapped_name[0]
                gene_names_counts[mapped_name[0]] += 1

    print "Number of events names for RT-PCR in names-map-file: %d" % len(names_pcr2majiq_dict.keys())
    print "Number of events names for MAJIQ in names-map-file: %d" % len(names_junc_majiq.keys())

    # Parse RT-PCR results
    pcr_rest_stim_delta = defaultdict(list)
    with open(args.pcr, 'r') as pcr_file:
        for i, pcr_line in enumerate(pcr_file):
            if i<1: continue  # headers
            pcr_fields = pcr_line.rstrip().split()

            if pcr_fields[0] not in names_pcr2majiq_dict.keys(): continue  # If the event is not in the list of events selected, skip it
            try:
                pcr_rest_stim_delta[pcr_fields[0]].append((float(pcr_fields[1]) - float(pcr_fields[2]))/100)
            except IndexError:
                print "%s value not found, assigned 0..." % pcr_fields[0]
                pcr_rest_stim_delta[pcr_fields[0]] = [0]
            except ValueError, e:
                print e.message
                print "Event PSI could not be converted into float, wrong parsing??"
                pcr_rest_stim_delta[pcr_fields[0]] = [0, 0]

    print "Number of events found in RT-PCR file: %d" % len(pcr_rest_stim_delta)

    # Process MAJIQ files for resting RNA-Seq data
    majiq_delta_dict = defaultdict(list)

    if os.path.isdir(args.majiq_deltas[0]):
        majiq_delta_files = []
        for file in os.listdir(args.majiq_deltas[0]):
            if file.endswith("deltamatrix.pickle"):
                majiq_delta_files.append(args.majiq_deltas[0]+'/'+file)
    else:
        majiq_delta_files = args.majiq_deltas

    for mfile in majiq_delta_files:
        # expected_delta_psi = []
        with open(mfile) as mfile_open:
            mpickle = pickle.load(mfile_open)
            for i, lsv_info in enumerate(mpickle[1]):
                if len(mpickle[0][i])>2:
                    continue
                # expected_delta_psi.append(delta_expected_psi(mpickle[0][i][0]))
                if lsv_info[1] not in names_junc_majiq.keys(): continue # If the event is not in the list of events selected, skip it

                # Find all junctions from that LSV that are included
                for j, lsv_way in enumerate(mpickle[0][i]):
                    if j in names_junc_majiq[lsv_info[1]]:
                        majiq_delta_dict[lsv_info[1]+"#"+str(j)].append(mpickle[0][i][j])  # add bins from delta psis
        # barchart_expected(expected_psis, args.plotpath, mfile)
        print "Number of events found in %s: %d" % (mfile, len(majiq_delta_dict.keys()))


    ## Intersect names from RT-PCR and MAJIQ
    common_names_set = set([names_pcr2majiq_dict[k] for k in pcr_rest_stim_delta.keys()]).intersection(set(majiq_delta_dict.keys()))

    rt_pcr = []
    majiq = []

    flipped_thres = 1.0
    flipped_lsv_dict = defaultdict(str)  # List of strings with flipped LSVs info
    for common_name in common_names_set:
        for name_majiq, name_pcr in names_majiq2pcr_dict.iteritems():
            if names_majiq2pcr_dict[common_name] == name_pcr:

                name = name_majiq
                if not len(majiq_delta_dict[name]):
                    print "%s - %s: Not found in majiq_dict" % (name_majiq, name_pcr)
                    continue
                # For Majiq, compute mean over Expected PSIs
                majiq_delta = avg_expected_delta_psi(majiq_delta_dict[name])
                min_delta = 10
                rtpcr_delta = pcr_rest_stim_delta[names_majiq2pcr_dict[name]][0]
                for rtpcr_delta_entry in pcr_rest_stim_delta[names_majiq2pcr_dict[name]]:
                    if abs(rtpcr_delta_entry - majiq_delta) < min_delta:
                        rtpcr_delta = rtpcr_delta_entry
                # check if event has expected PSIs suspicious of being flipped
                if abs(majiq_delta - rtpcr_delta) > flipped_thres:
                    flipped_lsv_dict[name] = "%s\t%s\t%f\t%f\t%d\t%d\n" % (names_majiq2pcr_dict[name], name, rtpcr_delta, majiq_delta, int(gene_names_counts[names_majiq2pcr_dict[name]]<2), int(len(names_junc_majiq[str(name).split('#')[0]])<2) )
                    continue

                print "Found: %s - %s" % (name, names_majiq2pcr_dict[name])
                rt_pcr.append(rtpcr_delta)
                majiq.append(majiq_delta)

    plot_rtpcr_majiq_boxplot(rt_pcr, majiq, args.plotpath)

    print "%d elements flipped" % len(flipped_lsv_dict.keys())
    # Save presumably flipped events
    flipped_lsv_names = [str(name_str).split("#")[0] for name_str in flipped_lsv_dict.keys()]
    with open('flipped_events_delta_threshold_%.2f.txt' % flipped_thres, 'w') as flipped_file:
        flipped_file.write("name_rtpcr\tname_majiq\trtpcr_delta\tmajiq_delta\tunique_rtpcr_name?\tlsv_1_way?\ttarget_junction_coord\n")
        # with open(args.majiq_builder_files) as majiq_builder_files:
        #     majiq_builder = pickle.load(majiq_builder_files)
        #     for lsv in majiq_builder[1]:
        #         if lsv.id in flipped_lsv_names:
        #             for flip_key in flipped_lsv_dict:
        #                 if lsv.id in flip_key:
        #                     flipped_lsv_dict[flip_key] = flipped_lsv_dict[flip_key] + "\t%s\n" % str(lsv.junction_id[int(str(flip_key).split('#')[1])])
        for flipped_lsv in flipped_lsv_dict:
            flipped_file.write(flipped_lsv_dict[flipped_lsv])


if __name__ == '__main__':
    main()