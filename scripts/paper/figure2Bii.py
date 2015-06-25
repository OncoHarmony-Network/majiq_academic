__author__ = 'abarrera'
import matplotlib as mpl
mpl.use('Agg')
from scripts.utils import save_or_show

from collections import defaultdict
import argparse
import cPickle as pkl
import numpy as np
import os
from matplotlib import pyplot as plt

def load_coverage(cov_file_list, cov_suffix):
    coverage_list = []
    try:
        for file in cov_file_list:
            coverage_list.append(pkl.load(open(file+cov_suffix)))
    except IOError:
        print "[INFO]::Coverage not found for %s..." % file
        return []

    return coverage_list

def coverage_from_file(file, cov_suffix):
    if os.path.isfile(file+cov_suffix):
        return pkl.load(open(file+cov_suffix))

    result_dict = defaultdict(list)  # First the num. reads, second the positions
    with open(file) as majiq_file:
        majiq_builder = pkl.load(majiq_file)
        for lsv in majiq_builder[1]:
            num_reads = 0
            num_pos = 0
            for i, junc in enumerate(lsv.junction_list):
               num_reads += np.sum(junc.data[0])
               num_pos += junc.nnz
            result_dict[lsv.id].append([num_reads, num_pos])
    # Save results
    pkl.dump(result_dict, open(file+cov_suffix, 'w'))
    print "Coverage saved in: %s" % file+cov_suffix
    return result_dict


def load_clean_reads(file):
    with open(file) as clean_f:
        lsv_l = pkl.load(clean_f)
        return dict(lsv_l)


def cov_combined(coverages):
    common_dict = defaultdict(list)
    for cov_dict in coverages:
        for k, v in cov_dict.iteritems():
            common_dict[k].append(float(sum(v)))
    for k in common_dict:
        common_dict[k] = np.mean(np.array(common_dict[k]))
    # print repr(common_dict)
    return common_dict


def plot_bins(read_bins, names, rtpcr, plotpath, V=0.1, extension='png'):

    V *= 100
    def autolabel(rects, rtpcrs):
        for ii, rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., 0.9*height, '%d/%d'% (np.count_nonzero(np.array(rtpcrs[ii])>V), len(rtpcrs[ii])),
                    ha='center', va='bottom', color='white', size='large')

    repro_lsvs = [(sum(r)*1.0/max(1.0, len(r))*100) for r in read_bins]

    plotname="Reproducibility divided by coverage (N=%d)." % (np.sum([len(i) for i in read_bins]))
    plotname+="\n%.0f%% of LSVs validated (V=%.1f)" % (100*(1.*sum([np.count_nonzero(np.array(pcr_cov)>V) for pcr_cov in rtpcr])/sum([len(pcr_cov) for pcr_cov in rtpcr])), V/100)
    n_groups = len(names)

    fig, ax = plt.subplots()

    bar_width = .9
    index = np.arange(n_groups)+.05


    opacity = 0.5

    rects1 = plt.bar(index, repro_lsvs, bar_width,
             alpha=opacity,
             color='b',
             label='MAJIQ lsvs')

    autolabel(rects1, rtpcr)

    plt.xlabel('# Positions with reads')
    plt.ylabel('% Reproduced')
    plt.ylim([0, 100])
    plt.title(plotname)
    # plt.grid()
    # plt.title('Delta in expected PSI between %s.' % (replica_names_joint))
    plt.xticks(index + bar_width/2.0, names)
    # plt.legend(loc=2)

    plt.tight_layout()
    save_or_show(plotpath,
                  plotname.lower().replace('\n', ' - ').replace(' ', '_').replace('(', '').replace(')', '').replace('=', ''),
                  exten=extension)


def get_rtpcr_std(lcer, lliv):
    deltas=[]
    for c in lcer:
        for l in lliv:
            deltas.append(c-l)
    return np.std(deltas)


def main():

    EXTENSION_TYPES = ['png', 'pdf']

    """Distribution of reproduced LSVs by coverage"""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('events', type=str, help='Events files from reproducibility ranking.')
    parser.add_argument('--builder-files', required=True, dest='builder_f', nargs='+', help='MAJIQ builder files used to compute lsv coverage.')
    parser.add_argument('--rtpcr-files', dest='rtpcr_f', help='RT-PCR results for MAJIQ')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--V', default=.1, type=float, action='store', help='Output file path')
    parser.add_argument('--extension', default=EXTENSION_TYPES[1], choices=EXTENSION_TYPES, help='Extension of the created figure (%s).' % ', '.join(EXTENSION_TYPES))
    args = parser.parse_args()

    lsvs = pkl.load(open(args.events))
    coverages = []
    for majiq_f in args.builder_f:
        # coverages.append(coverage_from_file(majiq_f, '.coverage'))
        coverages.append(load_clean_reads(majiq_f))
    common_cov = cov_combined(coverages)

    BINNAMES = ['0-15', '15-20', '20-40', '40-100', '100-xx']
    ranges = [15, 20, 40, 100, 900000]

    read_bins = [[],[],[],[],[]]
    rtpcr_bins = [[],[],[],[],[]]

    rtpcr = defaultdict()
    if args.rtpcr_f:
        with open(args.rtpcr_f) as rtpcr_f:
            for l in rtpcr_f:
                if l.startswith("#"): continue
                fields = l.split()
                try:
                    rtpcr[fields[0].split("#")[0]]=[
                        float(fields[3])-float(fields[4]),
                        # float(fields[5])-float(fields[6]),
                        fields[1],
                        get_rtpcr_std([float(ss) for ss in fields[-6:-4]], [float(dd) for dd in fields[-3:-1]])
                    ]
                except ValueError:
                    print "[WARNING] :: Couldn't parse %s" % (fields[1])
    print "Intersection of LSVs between MAJIQ and RT-PCR %d" % len(set([ll[0][0].split("#")[0] for ll in lsvs]).intersection(rtpcr))

    for lsv in lsvs:
        for i, ran in enumerate(ranges):
            lsv_id = lsv[0][0].split("#")[0]
            if lsv_id in common_cov.keys():
                if common_cov[lsv_id] < ran:
                    read_bins[i].append(lsv[1]) # event[0][0][1]
                    if args.rtpcr_f and lsv_id in rtpcr:
                        rtpcr_bins[i].append(abs(rtpcr[lsv_id][0])) #+rtpcr[lsv_id][2]
                        print "%s\t%d\t%s" % (lsv_id, common_cov[lsv_id], rtpcr[lsv_id][1])
                    break
    for j, r in enumerate(read_bins):
        print "#Reads %s; from %d events, percentage reproduced %.2f%%" % (BINNAMES[j], len(r), (sum(r)*1.0/max(1.0, len(r))*100))

    BINNAMES = ["%s \n(n=%d)" % (a, b) for a,b in zip(BINNAMES, [len(aa) for aa in read_bins])]
    plot_bins(read_bins, BINNAMES, rtpcr_bins, args.output, args.V, extension=args.extension)


if __name__ == '__main__':
    main()
