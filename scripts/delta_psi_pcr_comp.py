# from __future__ import division
from collections import defaultdict
import pickle
import numpy as np
import argparse
import ast
import os
import matplotlib.pyplot as ppl
from scripts.utils import _save_or_show

__author__ = 'abarrera'


def get_brewer_color(n):
    BP = [
       [228,26,28],
       [55,126,184],
       [77,175,74],
       [152,78,163],
       [255,127,0],
       [255,255,51],
       [166,86,40],
       [247,129,191],
       [153,153,153],

       [28,126,128],
       [155,226,29],
       [177,275,19],
       [252,178,8],
       [55,227,100],
       [55,55,151],
       [266,186,140],
       [47,229,36],
       [253,253,253]
    ]
    return "#"+"".join(format(n, 'x') for n in BP[n])

def plot_rtpcr_majiq_boxplot(rt_pcr, majiq, coverage, plotpath):
    #figure out how many groups of events exist

    LOW_COV = 15
    MED_COV = 40
    coverage = np.array(coverage)
    low_cov_mask = (coverage <= LOW_COV)
    med_cov_mask = (coverage > LOW_COV) & (coverage <= MED_COV)
    high_cov_mask = (coverage > MED_COV)

    fig, axx = ppl.subplots(2, 2, sharex=True, sharey=True, figsize=[12, 12], dpi=300)
    ppl.suptitle("Delta PSI comparison: RT-PCR Vs MAJIQ (N=%d)" % len(rt_pcr))
    diagonal = np.linspace(-1, 1, num=len(rt_pcr))

    # All
    axx[0][0].scatter(np.array(majiq)[low_cov_mask], np.array(rt_pcr)[low_cov_mask], c=get_brewer_color(0), alpha=.5, s=50)
    axx[0][0].scatter(np.array(majiq)[med_cov_mask], np.array(rt_pcr)[med_cov_mask], c=get_brewer_color(1), alpha=.5, s=50)
    axx[0][0].scatter(np.array(majiq)[high_cov_mask], np.array(rt_pcr)[high_cov_mask], c=get_brewer_color(2), alpha=.5, s=50)
    axx[0][0].set_ylim([-1, 1])
    axx[0][0].set_xlim([-1, 1])
    axx[0][0].set_xlabel('MAJIQ')
    axx[0][0].set_ylabel('RT-PCR')
    axx[0][0].set_title('All')
    fit = np.polyfit(majiq, rt_pcr, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    axx[0][0].plot(majiq, fit_fn(majiq), '--k')
    axx[0][0].plot(diagonal, diagonal, '--', color="#cccccc")

    majiq_arrays = [np.array(majiq)[low_cov_mask], np.array(majiq)[med_cov_mask], np.array(majiq)[high_cov_mask]]
    rtpcr_arrays = [np.array(rt_pcr)[low_cov_mask], np.array(rt_pcr)[med_cov_mask], np.array(rt_pcr)[high_cov_mask]]
    titles = [
        "Cov <= %d (N=%d)" % (LOW_COV, np.count_nonzero(low_cov_mask)),
        "%d < Cov <= %d (N=%d)" % (LOW_COV, MED_COV, np.count_nonzero(med_cov_mask)),
        "Cov > %d (N=%d)" % (MED_COV, np.count_nonzero(high_cov_mask)),
    ]

    for n in xrange(3):
        if n==0:
            x=0
            y=1
        else:
            x=1
            y=n-1
        axx[x][y].scatter(majiq_arrays[n], rtpcr_arrays[n], c=get_brewer_color(n), alpha=.5, s=50)
        axx[x][y].set_ylim([-1, 1])
        axx[x][y].set_xlim([-1, 1])
        axx[x][y].set_xlabel('MAJIQ')
        axx[x][y].set_ylabel('RT-PCR')
        axx[x][y].set_title(titles[n])
        fit = np.polyfit(majiq_arrays[n], rtpcr_arrays[n], 1)
        fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
        axx[x][y].plot(majiq_arrays[n], fit_fn(majiq_arrays[n]), '--k')
        axx[x][y].plot(diagonal, diagonal, '--', color="#cccccc")

    print majiq, rt_pcr
    # fit = np.polyfit(majiq, rt_pcr, 1)
    # fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    # ppl.plot(majiq, fit_fn(majiq), '--k')

    # diagonal = np.linspace(-1, 1, num=len(rt_pcr))
    # ppl.plot(diagonal, diagonal, '--', color="#cccccc")
    # ppl.ylim([-1, 1])
    # pyplot.plot(majiq, rt_pcr, '.', color='r')
    #
    # ppl.xlabel('MAJIQ')
    # ppl.ylabel('RT-PCR')
    _save_or_show(plotpath, "delta_psi_comparison_rtpcr_majiq")


def barchart_expected(expected_psis, plotpath, mfile):
    ppl.figure()
    ppl.hist(expected_psis, bins=40)
    name = os.path.basename(mfile)
    _save_or_show(plotpath, name +'_expected_dist')


def collapse_matrix(matrix):
    "Collapse the diagonals probabilities in 1-D and return them"
    collapsed = []
    matrix_corner = matrix.shape[0]
    for i in xrange(-matrix_corner, matrix_corner):
        collapsed.append(np.diagonal(matrix, offset=i).sum())
    return np.array(collapsed)


def delta_expected_psi(bins, file_name):
    step = 2.0 / bins.size
    projection_prod = bins * np.arange(-1+step/2, 1, step)
    print "%s: %f" % (file_name, np.sum(projection_prod))

    return np.sum(projection_prod)


def avg_expected_delta_psi(bins_list):
    return np.mean([delta_expected_psi(collapse_matrix(bins[0]), bins[1]) for bins in bins_list])


def coverage_from_file(file, cov_suffix):
    result_dict = defaultdict(list)  # First the num. reads, second the positions
    with open(file) as majiq_file:
        majiq_builder = pickle.load(majiq_file)
        for lsv in majiq_builder[1]:
            for i, junc in enumerate(lsv.junction_list):
                result_dict[lsv.id+"#"+str(i)].append([np.sum(junc.data[0]), junc.nnz])
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
    jobs_list = []

    for rest_file in rest_file_list:
        jobs_list.append(pool.apply_async(coverage_from_file, args=[rest_file, COV_SUFFIX]))

    for stim_file in stim_file_list:
        jobs_list.append(pool.apply_async(coverage_from_file, args=[stim_file, COV_SUFFIX]))

    pool.close()
    pool.join()

    coverage_lists.append(load_coverage(rest_file_list, COV_SUFFIX))
    coverage_lists.append(load_coverage(stim_file_list, COV_SUFFIX))
    return coverage_lists


def get_min_coverage(coverage_list, name):
    return min(np.mean([c[name][0] for c in coverage_list[0]]), np.mean([c[name][0] for c in coverage_list[1]]))


def main():
    parser = argparse.ArgumentParser(description="Compare Delta PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--majiq-deltas", required=True, dest='majiq_deltas', nargs='+', help='MAJIQ Delta PSI predictions for resting and stimuli RNA-Seq data.')
    parser.add_argument("--names-map-file", required=True, dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')
    parser.add_argument("--builder-files", required=True, dest='majiq_builder_files', help='File containing MAJIQ builder output file.')
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    ## Analyze coverage, saving it for later
    coverage_list = get_coverage(args.majiq_builder_files)  # List of lists of dicts

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
                        majiq_delta_dict[lsv_info[1]+"#"+str(j)].append([mpickle[0][i][j], mfile])  # add bins from delta psis
        # barchart_expected(expected_psis, args.plotpath, mfile)
        print "Number of events found in %s: %d" % (mfile, len(majiq_delta_dict.keys()))


    ## Intersect names from RT-PCR and MAJIQ
    common_names_set = set([names_pcr2majiq_dict[k] for k in pcr_rest_stim_delta.keys()]).intersection(set(majiq_delta_dict.keys()))

    rt_pcr = []
    majiq = []
    coverage = []

    flipped_thres = 1.0
    flipped_lsv_dict = defaultdict(str)  # List of strings with flipped LSVs info
    for common_name in common_names_set:
        for name_majiq, name_pcr in names_majiq2pcr_dict.iteritems():
            if names_majiq2pcr_dict[common_name] == name_pcr:

                name = name_majiq
                if not len(majiq_delta_dict[name]):
                    print "%s - %s: Not found in majiq_dict" % (name_majiq, name_pcr)
                    continue
                print "Found: %s - %s" % (name, names_majiq2pcr_dict[name])

                # For Majiq, compute mean over Expected PSIs
                majiq_delta = avg_expected_delta_psi(majiq_delta_dict[name])
                min_delta = 10
                rtpcr_delta = pcr_rest_stim_delta[names_majiq2pcr_dict[name]][0]
                for rtpcr_delta_entry in pcr_rest_stim_delta[names_majiq2pcr_dict[name]]:
                    if abs(rtpcr_delta_entry - majiq_delta) < min_delta:
                        rtpcr_delta = rtpcr_delta_entry
                # check if event has expected PSIs suspicious of being flipped
                if abs(majiq_delta - rtpcr_delta) > .5:
                    print "%s has an abs. expected difference in delta psi of %.2f. Coverage: %d " % (name, abs(majiq_delta - rtpcr_delta), get_min_coverage(coverage_list, name))
                if abs(majiq_delta - rtpcr_delta) > flipped_thres:
                    flipped_lsv_dict[name] = "%s\t%s\t%f\t%f\t%d\t%d\n" % (names_majiq2pcr_dict[name], name, rtpcr_delta, majiq_delta, int(gene_names_counts[names_majiq2pcr_dict[name]]<2), int(len(names_junc_majiq[str(name).split('#')[0]])<2) )
                    continue


                rt_pcr.append(rtpcr_delta)
                majiq.append(majiq_delta)
                coverage.append(get_min_coverage(coverage_list, name))

    plot_rtpcr_majiq_boxplot(rt_pcr, majiq, coverage, args.plotpath)

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