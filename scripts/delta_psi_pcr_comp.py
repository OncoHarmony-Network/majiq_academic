from collections import defaultdict
import pickle
import ast
import os
import argparse
import numpy as np
import matplotlib.pyplot as ppl
import numpy.ma as ma
import scripts.utils
from scripts.utils import coverage_from_file
from voila.vlsv import collapse_matrix

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


def plot_rtpcr_majiq_var(rt_pcr, majiq, coverage, delta_delta_psi, plotpath, names=None):
    """Analyze:
        -   Delta(PSI),
        -   Prob(Delta(PSI))>V) and
        -   Num. of reads.
    Create a subplot comparing 2 of these magnitudes, leaving 1 out (codified by color).
    Notes:
        COV = Coverage
        DE = Delta of Expected psi
        PIED = Probability of Inaccurate Expected Delta psi
    """
    LOW_COV = 15
    MED_COV = 40
    low_cov_mask = (coverage <= LOW_COV)
    med_cov_mask = (coverage > LOW_COV) & (coverage <= MED_COV)
    high_cov_mask = (coverage > MED_COV)

    LOW_DE = .10
    MED_DE = .20

    delta_psi = majiq
    print "Delta Expected PSI: ", delta_psi
    low_de_mask = ma.mask_or((delta_psi >= -LOW_DE), (delta_psi <= LOW_DE))
    med_de_mask = ma.mask_or((delta_psi < -LOW_DE) & (delta_psi >= -MED_DE),  (delta_psi > LOW_DE) & (delta_psi <= MED_DE))
    high_de_mask = (delta_psi > MED_DE)

    print "Delta of Delta Expected PSI: ", delta_delta_psi
    LOW_PIED = .10
    MED_PIED = .20
    low_pied_mask = (delta_delta_psi <= LOW_PIED)
    med_pied_mask = (delta_delta_psi > LOW_PIED) & (delta_delta_psi <= MED_PIED)
    high_pied_mask = (delta_delta_psi > MED_PIED)


    fig, axx = ppl.subplots(3, 4, figsize=[22, 10], dpi=300)
    ppl.suptitle("Delta PSI comparison: RT-PCR Vs MAJIQ (N=%d) - Delta(PSI); Prob(Delta(PSI))>0.2); #Reads " % len(rt_pcr))
    diagonal = np.linspace(-1, 1, num=len(rt_pcr))

    masks = [
        [low_de_mask, med_de_mask, high_de_mask],
        [low_pied_mask, med_pied_mask, high_pied_mask],
        [low_cov_mask, med_cov_mask, high_cov_mask]
    ]
    titles = [
        "#Reads Vs Delta(Delta(PSI))",
        "#Reads Vs Delta(PSI)",
        "Delta(Delta(PSI)) Vs Delta(PSI)"
    ]
    pairs = [
        [coverage, delta_delta_psi],
        [coverage, delta_psi],
        [delta_delta_psi, delta_psi]
    ]

    pairs_labels = [
        ["#Reads", "Prob(E(Delta(PSI)>0.2)) [DD]"],
        ["#Reads", "E(Delta(PSI)) [E]"],
        ["Prob(E(Delta(PSI)>0.2)) [DD]", "E(Delta(PSI)) [E]"]
    ]

    color_labels = [
        "|E|",
        "DD",
        "#R"
    ]

    thres = [
        [LOW_DE, MED_DE],
        [LOW_PIED, MED_PIED],
        [LOW_COV, MED_COV]
    ]

    pairs_limits = [
        [[0, 200], [0, 1]],
        [[0, 200], [-1, 1]],
        [[0, 1], [-1, 1]]
    ]

    for n in xrange(3):
        if n == 1:
            print repr(np.array(names)[masks[n][2]])
            print repr(coverage[masks[n][2]])

        axx[n][0].scatter(pairs[n][0][masks[n][0]], pairs[n][1][masks[n][0]], marker='o', c=get_brewer_color(0), alpha=.8, s=50, label="%s < %.1f" % (color_labels[n], thres[n][0]))
        axx[n][0].scatter(pairs[n][0][masks[n][1]], pairs[n][1][masks[n][1]], marker='o', c=get_brewer_color(1), alpha=.8, s=50, label=" %.1f < %s < %.1f" % (thres[n][0], color_labels[n], thres[n][1]))
        axx[n][0].scatter(pairs[n][0][masks[n][2]], pairs[n][1][masks[n][2]], marker='o', c=get_brewer_color(2), alpha=.8, s=50, label="%s > %.1f" % (color_labels[n], thres[n][1]))
        axx[n][0].set_xlim(pairs_limits[n][0])
        axx[n][0].set_ylim(pairs_limits[n][1])
        axx[n][0].set_xlabel(pairs_labels[n][0])
        axx[n][0].set_ylabel(pairs_labels[n][1])
        # axx[n].set_title(titles[n])
        axx[n][0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., prop={'size':8})
        fit = np.polyfit(pairs[n][0], pairs[n][1], 1)
        fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
        axx[n][0].plot(pairs[n][0], fit_fn(pairs[n][0]), '--k')
        # axx[n].plot(diagonal, diagonal, '--', color="#cccccc")
        for j in xrange(1,4):
            axx[n][j].scatter(pairs[n][0][masks[n][j-1]], pairs[n][1][masks[n][j-1]], marker='o', c=get_brewer_color(j-1), alpha=1, s=50)
            axx[n][j].set_xlim(pairs_limits[n][0])
            axx[n][j].set_ylim(pairs_limits[n][1])
            axx[n][j].set_xlabel(pairs_labels[n][0])
            axx[n][j].set_ylabel(pairs_labels[n][1])
            # axx[n][j].set_title(titles[n])
            try:
                fit = np.polyfit(pairs[n][0][masks[n][j-1]], pairs[n][1][masks[n][j-1]], 1)
                fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
                axx[n][j].plot(pairs[n][0][masks[n][j-1]], fit_fn(pairs[n][0][masks[n][j-1]]), '--k')
            except TypeError, e:
                print e.message


    print majiq, rt_pcr
    scripts.utils.save_or_show(plotpath, "delta_psi_comparison_with_FN") #


def plot_rtpcr_majiq_reads_deltas(rt_pcr, majiq, coverage, delta_delta_psi, plotpath):

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
    axx[0][0].scatter(np.array(majiq)[low_cov_mask], np.array(rt_pcr)[low_cov_mask], c=get_brewer_color(0), alpha=.5 )
    axx[0][0].scatter(np.array(majiq)[med_cov_mask], np.array(rt_pcr)[med_cov_mask], c=get_brewer_color(1), alpha=.5)
    axx[0][0].scatter(np.array(majiq)[high_cov_mask], np.array(rt_pcr)[high_cov_mask], c=get_brewer_color(2), alpha=.5)
    axx[0][0].set_ylim([-1, 1])
    axx[0][0].set_xlim([-1, 1])
    axx[0][0].set_xlabel('MAJIQ')
    axx[0][0].set_ylabel('RT-PCR')
    axx[0][0].set_title('All')
    fit = np.polyfit(majiq, rt_pcr, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    axx[0][0].plot(majiq, fit_fn(majiq), '--k')
    axx[0][0].plot(diagonal, diagonal, '--', color="#cccccc")

    masks = [low_cov_mask, med_cov_mask, high_cov_mask]
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
        axx[x][y].scatter(np.array(majiq)[masks[n]], np.array(rt_pcr)[masks[n]], c=get_brewer_color(n), alpha=.5)
        axx[x][y].set_ylim([-1, 1])
        axx[x][y].set_xlim([-1, 1])
        axx[x][y].set_xlabel('MAJIQ')
        axx[x][y].set_ylabel('RT-PCR')
        axx[x][y].set_title(titles[n])
        fit = np.polyfit(np.array(majiq)[masks[n]], np.array(rt_pcr)[masks[n]], 1)
        fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
        axx[x][y].plot(np.array(majiq)[masks[n]], fit_fn(np.array(majiq)[masks[n]]), '--k')
        axx[x][y].plot(diagonal, diagonal, '--', color="#cccccc")

    print majiq, rt_pcr
    scripts.utils.save_or_show(plotpath, "delta_psi_comparison_rtpcr_majiq")


def plot_rtpcr_majiq_miso(rt_pcr, majiq, miso, plotpath):

    fig, (ax1, ax2) = ppl.subplots(1, 2, sharey=True, figsize=[12, 6], dpi=300)
    ppl.suptitle("Delta PSI comparison: RT-PCR Vs (MAJIQ N=%d; MISO N=%d)" % (len(majiq), len(miso)))
    diagonal = np.linspace(-1, 1, num=len(rt_pcr))

    # All
    ax1.scatter(majiq, rt_pcr, c=get_brewer_color(0), alpha=.5, s=50, label='MAJIQ')
    ax2.scatter(miso, rt_pcr, c=get_brewer_color(1), alpha=.5, s=50, label='MISO')
    ax1.set_ylim([-1, 1])
    ax2.set_ylim([-1, 1])
    ax1.set_xlim([-1, 1])
    ax2.set_xlim([-1, 1])
    ax1.set_ylabel('RT-PCR')
    ax1.set_xlabel('MAJIQ')
    ax2.set_xlabel('MISO')
    fit = np.polyfit(majiq, rt_pcr, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    ax1.plot(majiq, fit_fn(majiq), '--k')
    fit = np.polyfit(miso, rt_pcr, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    ax2.plot(miso, fit_fn(miso), '--k')

    ax2.plot(majiq, fit_fn(majiq), '--k')
    ax1.plot(diagonal, diagonal, '--', color="#cccccc")
    ax2.plot(diagonal, diagonal, '--', color="#cccccc")

    print majiq, miso, rt_pcr
    scripts.utils.save_or_show(plotpath, "delta_psi_comparison_rtpcr_majiq_miso")


def plot_rtpcr_majiq(rt_pcr, majiq, coverage, plotpath):
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
    scripts.utils.save_or_show(plotpath, "delta_psi_comparison_rtpcr_majiq")


def barchart_expected(expected_psis, plotpath, mfile):
    ppl.figure()
    ppl.hist(expected_psis, bins=40)
    name = os.path.basename(mfile)
    scripts.utils.save_or_show(plotpath, name +'_expected_dist')


def delta_expected_psi(bins, file_name):
    step = 2.0 / bins.size
    projection_prod = bins * np.arange(-1+step/2, 1, step)
    print "%s: %f" % (file_name, np.sum(projection_prod))

    return np.sum(projection_prod)


def avg_expected_delta_psi(bins_list):
    return np.mean([delta_expected_psi(collapse_matrix(bins[0]), bins[1]) for bins in bins_list])


def load_coverage(cov_file_list, cov_suffix, are_clean_reads):
    coverage_list = []
    try:
        for file in cov_file_list:
            if are_clean_reads:
                coverage_list.append(dict(pickle.load(open(file))))
            else:
                coverage_list.append(pickle.load(open(file + cov_suffix)))
    except IOError:
        print "[INFO]::Coverage not found for %s..." % file
        return []

    return coverage_list

def get_coverage(majiq_files_or_dir, are_clean_reads=False):

    cond1_flist = []
    cond2_flist = []
    cond_keyws = ['clean_reads.Cer_CT28_40_520', 'clean_reads.Liv_CT28_40_520']
    COV_SUFFIX = ".pkl"
    coverage_lists = []

    if are_clean_reads:
        if os.path.isdir(majiq_files_or_dir):
            for file in os.listdir(majiq_files_or_dir):
                if file.endswith(COV_SUFFIX) and cond_keyws[0] in file:
                    cond1_flist.append(majiq_files_or_dir+'/'+file)
                if file.endswith(COV_SUFFIX) and cond_keyws[1] in file:
                    cond2_flist.append(majiq_files_or_dir+'/'+file)
    else:
        if os.path.isdir(majiq_files_or_dir):
            for file in os.listdir(majiq_files_or_dir):
                if file.endswith(".majiq") and "toJuan.resting" in file:
                    cond1_flist.append(majiq_files_or_dir+'/'+file)
                if file.endswith(".majiq") and "toJuan.stim" in file:
                    cond2_flist.append(majiq_files_or_dir+'/'+file)
        else:
            for file in majiq_files_or_dir:
                if "toJuan.resting" in file:
                    cond1_flist.append(file)
                if "toJuan.stim" in file:
                    cond2_flist.append(file)

    coverage_lists.append(load_coverage(cond1_flist, COV_SUFFIX, are_clean_reads))
    coverage_lists.append(load_coverage(cond2_flist, COV_SUFFIX, are_clean_reads))

    if len(coverage_lists[0]) and len(coverage_lists[1]):
        return coverage_lists

    from multiprocessing import Pool
    pool = Pool()
    jobs_list = []

    for rest_file in cond1_flist:
        jobs_list.append(pool.apply_async(coverage_from_file, args=[rest_file, COV_SUFFIX]))

    for stim_file in cond2_flist:
        jobs_list.append(pool.apply_async(coverage_from_file, args=[stim_file, COV_SUFFIX]))

    pool.close()
    pool.join()

    coverage_lists.append(load_coverage(cond1_flist, COV_SUFFIX))
    coverage_lists.append(load_coverage(cond2_flist, COV_SUFFIX))
    return coverage_lists


def get_min_coverage(coverage_list, name, are_clean_reads=False):
    if are_clean_reads:
        name = str(name).split("#")[0]
        return min(np.mean([c[name] for c in coverage_list[0]]), np.mean([c[name] for c in coverage_list[1]]))
    return min(np.mean([c[name][0] for c in coverage_list[0]]), np.mean([c[name][0] for c in coverage_list[1]]))


def mean_prob_out_exp_psi(bins_list, rtpcr_delta):
    return np.mean([prob_out_expected_psi(collapse_matrix(bins[0]), rtpcr_delta) for bins in bins_list])


def prob_out_expected_psi(bins, rtpcr_delta, V=.2):
    """Calculate probability of delta psi outside the acceptable area"""
    step = 2.0 / bins.size
    left = 0
    right = bins.size*2 - 1
    for i, w in enumerate(np.arange(-1+step/2, 1, step)):
        if not left and w > (rtpcr_delta - abs(rtpcr_delta*V)):
            left = i-1
        if right == bins.size*2 and w > (rtpcr_delta + abs(rtpcr_delta*V)):
            right = i
    # print np.sum(bins[:left]), np.sum(bins[right:])
    return (np.sum(bins[:left]) + np.sum(bins[right:]))


def load_rtpcr_results(pcr_file):

    # Parse RT-PCR results
    pcr_rest_stim_delta = defaultdict(list)
    for i, pcr_line in enumerate(pcr_file):
        if i<1: continue  # headers
        pcr_fields = pcr_line.rstrip().split()
        if "target" in pcr_fields[0]:
            lsv_name = "%s" % (pcr_fields[0])
        else:
            lsv_name = "%s" % (pcr_fields[0])
        pcr_rest_stim_delta[lsv_name].append((float(pcr_fields[1]) - float(pcr_fields[2]))/100)

    return pcr_rest_stim_delta


def main():
    parser = argparse.ArgumentParser(description="Compare Delta PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--majiq-deltas", required=True, dest='majiq_deltas', nargs='+', help='MAJIQ Delta PSI predictions for resting and stimuli RNA-Seq data.')
    parser.add_argument("--miso-deltas", dest='miso_deltas', nargs='+', help='MISO Delta PSI predictions for resting and stimuli RNA-Seq data.')
    parser.add_argument("--names-map-file", dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')
    parser.add_argument("--builder-files", dest='majiq_builder_files', help='File containing MAJIQ builder output file.')
    parser.add_argument("--clean-reads", dest='clean_reads', action='store_true', default=False, help='Are you using clean reads?.')
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()


    # Process MISO files
    if args.miso_deltas:
        miso_delta_files = scripts.utils.list_files_or_dir(args.miso_deltas, ".miso_bf")
        miso_delta_dict = defaultdict(list)
        for mfile in miso_delta_files:
            scripts.utils.miso_delta_reader(mfile, complex_lsvs=False, result_dict=miso_delta_dict)

    ## Analyze coverage, saving it for later
    if args.majiq_builder_files:
        coverage_list = get_coverage(args.majiq_builder_files, args.clean_reads)  # List of lists of dicts
        coverage = []

    # Read name mapping file
    names_pcr2majiq_dict = {}  # key RT-PCR, value MAJIQ
    names_majiq2pcr_dict = {}
    names_junc_majiq = defaultdict(list)
    gene_names_counts = defaultdict(lambda: 0)

    if args.clean_reads:
        with open(args.pcr, 'r') as pcr_file:
            pcr_rest_stim_delta = load_rtpcr_results(pcr_file)

        # Process MAJIQ files for resting RNA-Seq data
        majiq_delta_dict = defaultdict(list)
        majiq_delta_files = scripts.utils.list_files_or_dir(args.majiq_deltas, "model2.pickle")

        for mfile in majiq_delta_files:
            # expected_delta_psi = []
            with open(mfile) as mfile_open:
                mpickle = pickle.load(mfile_open)
                for i, lsv_info in enumerate(mpickle[1]):
                    if len(mpickle[0][i])>2:
                        continue
                    if len(mpickle[0][i]) == 1:
                        mpickle[0][i].append(mpickle[0][i][-1].T)

                    # Find all junctions from that LSV that are included
                    for j, lsv_way in enumerate(mpickle[0][i]):
                        lsv_name = lsv_info[1] + "#" + str(j)
                        if lsv_name not in pcr_rest_stim_delta.keys(): continue # If the event is not in the list of events selected, skip
                        majiq_delta_dict[lsv_name].append([mpickle[0][i][j], mfile])  # add bins from delta psis
            # barchart_expected(expected_psis, args.plotpath, mfile)
            print "Number of events found in %s: %d" % (mfile, len(majiq_delta_dict.keys()))
        common_names_set = set(pcr_rest_stim_delta.keys()).intersection(set(majiq_delta_dict.keys()))

        rt_pcr = []
        majiq = []
        miso = []
        delta_delta_psi = []
        final_names = []


        with open("toJordi.txt", "w") as toJordi:
            flipped_thres = 0.6
            flipped_lsv_dict = defaultdict(str)  # List of strings with flipped LSVs info
            for name in common_names_set:
                print "Found: %s" % name

                # For Majiq, compute mean over Expected PSIs
                majiq_delta = avg_expected_delta_psi(majiq_delta_dict[name])
                min_delta = 10
                rtpcr_delta = pcr_rest_stim_delta[name][0]
                for rtpcr_delta_entry in pcr_rest_stim_delta[name]:
                    if abs(rtpcr_delta_entry - majiq_delta) < min_delta:
                        rtpcr_delta = rtpcr_delta_entry
                # check if event has expected PSIs suspicious of being flipped
                if abs(majiq_delta - rtpcr_delta) > .5:
                    try:
                        print "%s has an abs. expected difference in delta psi of %.2f. Coverage: %d " % (name, abs(majiq_delta - rtpcr_delta), get_min_coverage(coverage_list, name))
                    except UnboundLocalError, e:
                        print "%s has an abs. expected difference in delta psi of %.2f" % (name, abs(majiq_delta - rtpcr_delta))
                if abs(majiq_delta - rtpcr_delta) > flipped_thres:
                    flipped_lsv_dict[name] = "%s\t%f\t%f\t%d\n" % (name, rtpcr_delta, majiq_delta, int(gene_names_counts[name]<2))
                    continue

                rt_pcr.append(rtpcr_delta)
                majiq.append(majiq_delta)
                delta_delta_psi.append(mean_prob_out_exp_psi(majiq_delta_dict[name], rtpcr_delta))
                final_names.append(name)

                if args.miso_deltas:
                    miso.append(np.mean(miso_delta_dict[name]))

                try:
                    coverage.append(get_min_coverage(coverage_list, name, args.clean_reads))
                    if coverage[-1] > 15 and delta_delta_psi[-1]>.2:
                        toJordi.write("Suspicious guy: %s. Coverage=%d; Prob(Delta(PSI))>0.2=%.2f; Expected(Delta(PSI))=%.2f: RT-PCR=%.2f\n" % (name, coverage[-1], delta_delta_psi[-1], majiq[-1], rt_pcr[-1]))
                except NameError, e:
                    pass

    else:
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
        majiq_delta_files = scripts.utils.list_files_or_dir(args.majiq_deltas, "model2.pickle")

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
        print "Number of common names after intersection with MAJIQ: %d" % len(common_names_set)

    if args.miso_deltas:

        common_names_set = common_names_set.intersection(set([k for k in miso_delta_dict.keys() if len(miso_delta_dict[k]) > 0]))
        print "Number of common names after intersection with MISO: %d" % len(common_names_set)


        rt_pcr = []
        majiq = []
        miso = []
        delta_delta_psi = []
        final_names = []

        with open("toJordi.txt", "w") as toJordi:
            flipped_thres = 1.0
            flipped_lsv_dict = defaultdict(str)  # List of strings with flipped LSVs info
            for common_name in common_names_set:
                # for name_majiq, name_pcr in names_majiq2pcr_dict.iteritems():
                #     if names_majiq2pcr_dict[common_name] == name_pcr:
                # name = name_majiq

                name = common_name
                if not len(majiq_delta_dict[name]):
                    print "%s - %s: Not found in majiq_dict" % (common_name, common_name)
                    continue
                # print "Found: %s - %s" % (name, names_majiq2pcr_dict[name])

                # For Majiq, compute mean over Expected PSIs
                majiq_delta = avg_expected_delta_psi(majiq_delta_dict[name])
                min_delta = 10
                # rtpcr_delta = pcr_rest_stim_delta[names_majiq2pcr_dict[name]][0]
                rtpcr_delta = pcr_rest_stim_delta[name][0]
                # for rtpcr_delta_entry in pcr_rest_stim_delta[names_majiq2pcr_dict[name]]:
                for rtpcr_delta_entry in pcr_rest_stim_delta[name]:
                    if abs(rtpcr_delta_entry - majiq_delta) < min_delta:
                        rtpcr_delta = rtpcr_delta_entry
                # check if event has expected PSIs suspicious of being flipped
                if abs(majiq_delta - rtpcr_delta) > .5:
                    try:
                        print "%s has an abs. expected difference in delta psi of %.2f. Coverage: %d " % (name, abs(majiq_delta - rtpcr_delta), get_min_coverage(coverage_list, name))
                    except UnboundLocalError, e:
                        print "%s has an abs. expected difference in delta psi of %.2f" % (name, abs(majiq_delta - rtpcr_delta))
                if abs(majiq_delta - rtpcr_delta) > flipped_thres:
                    flipped_lsv_dict[name] = "%s\t%s\t%f\t%f\t%d\t%d\n" % (name, name, rtpcr_delta, majiq_delta, int(gene_names_counts[names_majiq2pcr_dict[name]]<2), int(len(names_junc_majiq[str(name).split('#')[0]])<2) )
                    continue

                rt_pcr.append(rtpcr_delta)
                majiq.append(majiq_delta)
                delta_delta_psi.append(mean_prob_out_exp_psi(majiq_delta_dict[name], rtpcr_delta))
                final_names.append(name)

                if args.miso_deltas:
                    if (np.count_nonzero(np.isnan(np.mean(miso_delta_dict[name])))):
                        print miso_delta_dict[name]
                    miso.append(np.mean(miso_delta_dict[name]))

                try:
                    coverage.append(get_min_coverage(coverage_list, name))
                    if coverage[-1] > 15 and delta_delta_psi[-1]>.2:
                        toJordi.write("Suspicious guy: %s. Coverage=%d; Prob(Delta(PSI))>0.2=%.2f; Expected(Delta(PSI))=%.2f: RT-PCR=%.2f\n" % (name, coverage[-1], delta_delta_psi[-1], majiq[-1], rt_pcr[-1]))
                except (NameError, KeyError), e:
                    pass

    if args.miso_deltas:
        plot_rtpcr_majiq_miso(rt_pcr, majiq, miso, args.plotpath)
    if args.majiq_builder_files:
        # plot_rtpcr_majiq(rt_pcr, majiq, coverage, args.plotpath)
        # plot_rtpcr_majiq_var(np.array(rt_pcr), np.array(majiq), np.array(coverage), np.array(delta_delta_psi), args.plotpath, final_names)
        plot_rtpcr_majiq_reads_deltas(rt_pcr, majiq, coverage, delta_delta_psi, args.plotpath)

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