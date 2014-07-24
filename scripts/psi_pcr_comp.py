# from __future__ import division
from collections import defaultdict
import pickle
import numpy
from scipy.stats import pearsonr
import argparse
import ast
import os
import matplotlib.pyplot as pyplot

__author__ = 'abarrera'


def _save_or_show(plotpath, name):
    if plotpath:
        if os.path.isdir(plotpath):
            plot_base_path, plot_name = plotpath, name
        else:
            plot_base_path, plot_name = os.path.split(plotpath)
            if not os.path.exists(plot_base_path):
                os.makedirs(plot_base_path)
            if not plot_name:
                plot_name = name
        pyplot.savefig("%s/%s.png"%(plot_base_path, plot_name), width=300, height=300, dpi=100)
        print "Saved in:\n%s/%s" % (plot_base_path, plot_name)

        pyplot.clf()
    else:
        pyplot.show()

def plot_rtpcr_majiq_boxplot(rt_pcr, majiq, plotpath):
    #figure out how many groups of events exist

    fig, (ax1, ax2) = pyplot.subplots(1, 2, sharex=True, sharey=True, figsize=[12, 6], dpi=300)
    fig.suptitle("PSI comparison: RT-PCR Vs MAJIQ (N=%d)" % len(rt_pcr[0]))

    diagonal = numpy.linspace(0, 1, num=len(rt_pcr[0]))
    print majiq[0], rt_pcr[0]
    fit = numpy.polyfit(majiq[0], rt_pcr[0], 1)
    fit_fn = numpy.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    ax1.plot(majiq[0], fit_fn(majiq[0]), '--k')
    ax1.plot(diagonal, diagonal, '--', color="#cccccc")
    ax1.plot(majiq[0], rt_pcr[0], '.', color='r')

    ax1.set_xlabel('MAJIQ')
    ax1.set_ylabel('RT-PCR')
    ax1.set_title('Resting')

    fit = numpy.polyfit(majiq[1], rt_pcr[1], 1)
    fit_fn = numpy.poly1d(fit)
    ax2.plot(diagonal, diagonal, '--', color="#cccccc")
    ax2.plot(majiq[1], rt_pcr[1], '.', color='b')
    ax2.plot(majiq[1], fit_fn(majiq[1]), '--k')

    ax2.set_xlabel('MAJIQ')
    ax2.set_title('Stimuli')

    _save_or_show(plotpath, "psi_comparison_rtpcr_majiq")


def barchart_expected(expected_psis, plotpath, mfile):
     f = pyplot.figure()
     pyplot.hist(expected_psis, bins=40)
     name = os.path.basename(mfile)
     _save_or_show(plotpath, name +'_expected_dist')


def expected_psi(bins):
    bins = numpy.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * numpy.arange(step / 2, 1, step)
    return numpy.sum(projection_prod)


def avg_expected_psi(bins_list):
    return numpy.mean([expected_psi(bins) for bins in bins_list])


def expected_avg_psi(bins_list):
    return expected_psi(numpy.mean(bins_list, axis=0))


def expected_delta_psi(bins_list, rtpcr_psi):
    step = 1.0 / 40
    bins_index = numpy.arange(step / 2, 1, step)
    return numpy.mean([numpy.sum([bins[i]*abs(psi - rtpcr_psi) for i, psi in enumerate(bins_index)]) for bins in bins_list])


def main():
    parser = argparse.ArgumentParser(description="Compare PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--majiq-rest", required=True, dest='majiq_rest', nargs='+', help='MAJIQ PSI predictions for resting RNA-Seq data.')
    parser.add_argument("--majiq-stim", required=True, dest='majiq_stim', nargs='+', help='MAJIQ PSI predictions for stimulated RNA-Seq data.')
    parser.add_argument("--names-map-file", required=True, dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')
    parser.add_argument("--builder-file", required=True, dest='majiq_builder_file', help='File containing MAJIQ builder output file.')
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    # Read name mapping file
    names_pcr2majiq_dict = {}  # key RT-PCR, value MAJIQ
    names_majiq2pcr_dict = {}
    names_junc_majiq = defaultdict(list)
    gene_names_counts = defaultdict(lambda : 0)
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

    print "Number of events in names-map-file RT-PCR: %d" % len(names_pcr2majiq_dict.keys())
    print "Number of events in names-map-file MAJIQ: %d" % len(names_junc_majiq.keys())

    # Parse RT-PCR results
    pcr_rest_stim = defaultdict(list)
    with open(args.pcr, 'r') as pcr_file:
        for i, pcr_line in enumerate(pcr_file):
            if i<1: continue  # headers
            pcr_fields = pcr_line.rstrip().split()

            if pcr_fields[0] not in names_pcr2majiq_dict.keys(): continue  # If the event is not in the list of events selected, skip it
            try:
                pcr_rest_stim[pcr_fields[0]] = [float(pcr_fields[1])/100, float(pcr_fields[2])/100]
            except IndexError:
                print "%s value not found, assigned 0..." % pcr_fields[0]
                pcr_rest_stim[pcr_fields[0]] = [0, 0]
            except ValueError, e:
                print e.message
                print "Event PSI could not be converted into float, wrong parsing??"
                pcr_rest_stim[pcr_fields[0]] = [0, 0]

    print "Number of events found in RT-PCR file: %d" % len(pcr_rest_stim)

    # Process MAJIQ files for resting RNA-Seq data
    majiq_rest_dict = defaultdict(list)
    for mfile in args.majiq_rest:
        expected_psis = []
        with open(mfile) as mfile_open:
            mpickle = pickle.load(mfile_open)
            for i, lsv_info in enumerate(mpickle[1]):
                if len(mpickle[0][i])<3:
                    expected_psis.append(expected_psi(mpickle[0][i][0]))
                if lsv_info[1] not in names_junc_majiq.keys(): continue # If the event is not in the list of events selected, skip it

                # Find all junctions from that LSV that are included
                for j, lsv_way in enumerate(mpickle[0][i]):
                    if j in names_junc_majiq[lsv_info[1]]:
                        majiq_rest_dict[lsv_info[1]+"#"+str(j)].append(mpickle[0][i][j])
        # barchart_expected(expected_psis, args.plotpath, mfile)
    print "Number of events found in MAJIQ rest RNA-Seq: %d" % len(majiq_rest_dict.keys())

    # Process MAJIQ files for stimuli RNA-Seq data
    majiq_stim_dict = defaultdict(list)
    for mfile in args.majiq_stim:
        expected_psis = []
        with open(mfile) as mfile_open:
            mpickle = pickle.load(mfile_open)
            for i, lsv_info in enumerate(mpickle[1]):
                if len(mpickle[0][i])<3:
                    expected_psis.append(expected_psi(mpickle[0][i][0]))
                if lsv_info[1] not in names_junc_majiq.keys(): continue # If the event is not in the list of events selected, skip it

                # Find all junctions from that LSV that are included
                for j, lsv_way in enumerate(mpickle[0][i]):
                    if j in names_junc_majiq[lsv_info[1]]:
                        majiq_stim_dict[lsv_info[1]+"#"+str(j)].append(mpickle[0][i][j])
                # majiq_stim_dict[lsv_info[1]+"#"+names_junc_majiq[lsv_info[1]]].append(get_mean_step(mpickle[0][i][int(names_junc_majiq[lsv_info[1]])]))
        # barchart_expected(expected_psis, args.plotpath, mfile)
    print "Number of events found in MAJIQ stim RNA-Seq: %d" % len(majiq_stim_dict.keys())

    ## Intersect names from RT-PCR and MAJIQ
    common_names_set = set([names_pcr2majiq_dict[k] for k in pcr_rest_stim.keys()]).intersection(set(majiq_rest_dict.keys())).intersection(set(majiq_stim_dict.keys()))

    rt_pcr = [[], []]
    majiq = [[], []]

    flipped_thres = .35
    flipped_lsv_dict = defaultdict(str)  # List of strings with flipped LSVs info
    for common_name in common_names_set:
        for name_majiq, name_pcr in names_majiq2pcr_dict.iteritems():
            if names_majiq2pcr_dict[common_name] == name_pcr:
                name = name_majiq

                # For Majiq, compute mean over Expected PSIs
                majiq_rest_stat = avg_expected_psi(majiq_rest_dict[name])
                majiq_stim_stat = avg_expected_psi(majiq_stim_dict[name])

                # check if event has expected PSIs suspicious of being flipped
                if abs(majiq_rest_stat - pcr_rest_stim[names_majiq2pcr_dict[name]][0]) > flipped_thres or abs(majiq_stim_stat - pcr_rest_stim[names_majiq2pcr_dict[name]][1]) > flipped_thres:
                    flipped_lsv_dict[name] = "%s\t%s\t%f\t%f\t%f\t%f\t%d\t%d\n" % (names_majiq2pcr_dict[name], name, pcr_rest_stim[names_majiq2pcr_dict[name]][0], majiq_rest_stat, pcr_rest_stim[names_majiq2pcr_dict[name]][1], majiq_stim_stat, int(gene_names_counts[names_majiq2pcr_dict[name]]<2), int(len(names_junc_majiq[str(name).split('#')[0]])<2) )
                    continue

                print "%s - %s" % (name, names_majiq2pcr_dict[name])
                print "---- RT-PCR ----"
                rt_pcr[0].append(pcr_rest_stim[names_majiq2pcr_dict[name]][0])
                rt_pcr[1].append(pcr_rest_stim[names_majiq2pcr_dict[name]][1])
                print "[Resting]:\t%f" % (pcr_rest_stim[names_majiq2pcr_dict[name]][0])
                print "[Stimuli]:\t%f" % (pcr_rest_stim[names_majiq2pcr_dict[name]][1])

                print "---- MAJIQ ----"
                print "[Resting]: Mean of expected:\t%f" % (float(majiq_rest_stat))
                print "[Resting]: Expected delta psi statistic: %f" % (float(expected_delta_psi(majiq_rest_dict[name], pcr_rest_stim[names_majiq2pcr_dict[name]][0])))
                print "[Stimuli]: Mean of expected:\t%f" % (float(majiq_stim_stat))
                print "[Stimuli]: Expected delta psi statistic: %f" % (float(expected_delta_psi(majiq_stim_dict[name], pcr_rest_stim[names_majiq2pcr_dict[name]][1])))
                # print "[Resting]: Mean of %s:\t%f" % (majiq_rest_dict[name], float(majiq_rest_stat))
                # print "[Stimuli]: Mean of %s:\t%f" % (majiq_stim_dict[name], float(majiq_stim_stat))
                majiq[0].append(majiq_rest_stat)
                majiq[1].append(majiq_stim_stat)
    # print repr(rt_pcr), repr(majiq)

    plot_rtpcr_majiq_boxplot(rt_pcr, majiq, args.plotpath)

    # Save presumably flipped events
    flipped_lsv_names = [str(name_str).split("#")[0] for name_str in flipped_lsv_dict.keys()]
    with open('flipped_events_threshold_%.2f.txt' % flipped_thres, 'w') as flipped_file:
        flipped_file.write("name_rtpcr\tname_majiq\trtpcr_rest%f\tmajiq_rest\trtpcr_stim\tmajiq_stim\tunique_rtpcr_name?\tlsv_1_way?\ttarget_junction_coord\n")
        with open(args.majiq_builder_file) as majiq_builder_file:
            majiq_builder = pickle.load(majiq_builder_file)
            for lsv in majiq_builder[1]:
                if lsv.id in flipped_lsv_names:
                    for flip_key in flipped_lsv_dict:
                        if lsv.id in flip_key:
                            flipped_lsv_dict[flip_key] += "\t%s" % str(lsv.junction_id[str(flip_key).split('#')[1]])
        for flipped_lsv in flipped_lsv_dict:
            flipped_file.write(flipped_lsv_dict[flipped_lsv])


if __name__ == '__main__':
    main()