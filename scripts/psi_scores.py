# from matplotlib import use
# use('Agg')

from collections import defaultdict
import argparse
from math import log

from scipy.stats import pearsonr
from pylab import *


"""
Calculate PSI values
"""
DEBUG = False
TESTBREAK = 50
BINS = linspace(0, 1, num=40)


def plot_PSIs1VsPSIs2(score1, score2, replica1_name, replica2_name, method1, method2, plotpath=None):
    """Compute 2 plots: one for the variance, one for the mean"""
    plotname = "%sVs%s -- Delta PSIs for %s - %s" % (method1, method2, replica1_name, replica2_name)

    total_psis = float(len(score1))

    print len(score1), len(score2)

    better_in_method1 = np.sum(array(score1) < array(score2))
    better_in_method2 = np.sum(array(score1) > array(score2))

    print "Better in method %s: %.2f%%" % (method1, (better_in_method1 / total_psis) * 100)
    print "Better in method %s: %.2f%%" % (method2, (better_in_method2 / total_psis) * 100)

    f = figure()
    ylabel(method2)
    xlabel(method1)

    title(plotname, fontsize=10)
    max_value = max(max(score1), max(score2))

    max_value = 1
    xlim(0, max_value)
    ylim(0, max_value)
    plot(score1, score2, 'b.')
    plot([0, max_value], [0, max_value], '--', color="#cccccc")

    pear, pvalue = pearsonr(score1, score2)
    r_squared = pear ** 2

    # axarr[name_i].text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, r'$R^2$: %.2f (p-value: %.2E)'%(r_squared, pvalue), fontsize=12, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})
    text(abs(max_value) * 0.1, max_value - abs(max_value) * 0.2, '%.2f%%' % ((better_in_method1 / total_psis) * 100),
         fontsize=16)
    text(abs(max_value) * 0.8, max_value - abs(max_value) * 0.8, '%.2f%%' % ((better_in_method2 / total_psis) * 100),
         fontsize=16)

    _save_or_show(plotpath, plotname)


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
        savefig("%s/%s.png" % (plotpath, plotname), bbox_inches='tight')
        clf()
    else:
        show()


def calculate_dkl(p, q):
    """
    Calculates distance between two distributions:

    Dkl(P|Q) = sum_i ln(P(i)/Q(i))*P(i)
    """
    pseudo = 0.001
    p += pseudo
    q += pseudo
    left = log(array(p / q))
    return (left * p).sum(axis=1)


def calculate_l1_expected(p, q):
    return sum(abs(p - q) * scripts.src.psi.BINS_CENTER)


def calculate_l1(p, q):
    return (abs(p - q)).sum()


def calculate_cov(psi_list1, psi_list2):
    return abs(np.array(psi_list1) - np.array(psi_list2)) / ((np.array(psi_list1) + np.array(psi_list2)) / 2.0)


def calculate_ead_simple(psi_list1, psi_list2):
    return abs(np.array(psi_list1) - np.array(psi_list2))


def main():
    """
    Script for testing MAJIQ against other algorithms for PSIs.

    - Notice that the first file is for MAJIQ results, whereas the second file is reserved for others (initially, MISO).
    - All LSVs 'ways' are equally considered. The order of each way within a LSV should be the same in MAJIQ and MISO.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('psivalues_met1', help='Path for psi pickles to evaluate')
    parser.add_argument('psivalues_met2', nargs='+', help='Path for psi pickles to evaluate')
    parser.add_argument('--name1', default='replica1')
    parser.add_argument('--name2', default='replica2')
    parser.add_argument('--plotpsidist', default=False, help="Plot the PSI distributions for ALL junctions. Slow.")
    parser.add_argument('--plotpath', default=None,
                        help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--output')
    parser.add_argument('--type', default='old_majiq-miso', choices=['old_majiq-miso', 'old_majiq-old_majiq'])
    args = parser.parse_args()

    # For debugging, method1 is MAJIQ, method2 is MISO
    psivalues = pickle.load(open(args.psivalues_met1))
    psi_met1_rep1 = psivalues[0][0]
    psi_met1_rep2 = psivalues[0][1]

    psi_names_met1 = defaultdict()

    lsv_types_dict = {
        's|1e1.1|1e2.1': 'SE',
        't|1e1.1|1e2.1': 'SE',
        's|1e1.1|1e1.2': 'A3SS',
        't|1e1.1|2e1.1': 'A3SS',
        't|1e1.1|1e1.2': 'A5SS',
        's|1e1.1|2e1.1': 'A5SS'
    }

    # Discard LSVs with only one PSI
    for i, psis_lsv_met1 in enumerate(psi_met1_rep1):
        if len(psis_lsv_met1) < 2 or len(psi_met1_rep2[i]) < 2:
            print "1-way LSV skipped..."
            continue  # TODO: check that skipping is not necessary. LSVs with only 1 PSI are wrong..
        if psivalues[1][i][2] not in lsv_types_dict.keys():
            continue
        psi_names_met1[psivalues[1][i][1]] = i

    # Method1 (MAJIQ) psi scores
    psi_list1_met1 = []
    psi_list2_met1 = []

    psi_lists_met2 = []

    debug_dict1 = {}
    debug_dict2 = {}

    debug_names_miso_list = defaultdict(list)

    miso_all = []
    if args.type == 'old_majiq-miso':
        for miso_file in args.psivalues_met2:
            miso_psis_dict = defaultdict()
            with open(miso_file, 'r') as miso_res:
                for miso_line in miso_res:
                    if miso_line.startswith("event"): continue
                    miso_fields = miso_line.split('\t')
                    if miso_fields[0] not in psi_names_met1:
                        continue
                    miso_psis_dict[miso_fields[0]] = miso_fields[1]
            miso_all.append(miso_psis_dict)

        miso_common_names = set(miso_all[0].keys()).intersection(miso_all[1].keys())
        for miso_psis_dict in miso_all:
            miso_psis = []
            for psi_name in sorted(psi_names_met1.keys()):
                if psi_name not in miso_common_names:
                    print "%s is not in all MISO replicates" % psi_name
                    del psi_names_met1[psi_name]
                    continue
                try:
                    miso_psis_values = [float(miso_psi) for miso_psi in miso_psis_dict[psi_name].split(",")]
                except KeyError, e:
                    print "LSV %s is in MAJIQ but not in MISO!" % e
                    del psi_names_met1[psi_name]
                    continue
                    # if len(miso_psis_values) < 2:
                    # del majiq_psi_names[psi_name]
                    # continue

                if len(miso_psis_values) == 1:
                    # print miso_psis_values[0], 1 - miso_psis_values[0]
                    miso_psis_values.append(1.0 - miso_psis_values[0])
                debug_dict1[psi_name] = len(miso_psis_values)
                miso_psis.extend(miso_psis_values)
                debug_names_miso_list[psi_name].append(miso_psis_values)

            psi_lists_met2.append(miso_psis)

    elif args.type == 'old_majiq-old_majiq':
        psi_names_met2 = defaultdict()
        psi_list1_met2 = []
        psi_list2_met2 = []

        psivalues_met2 = pickle.load(open(args.psivalues_met2[0]))
        psi_met2_rep1 = psivalues_met2[0][0]
        psi_met2_rep2 = psivalues_met2[0][1]

        set_common_names = set([n[1] for n in psivalues_met2[1]]).intersection(set([n[1] for n in psivalues[1]]))

        for i, psis_lsv_met2 in enumerate(psi_met2_rep1):
            if psivalues_met2[1][i][1] not in set_common_names:
                print "[WARNING] :: %s in replica 1 but not in replica 2" % psivalues_met2[1][i][1]
                continue
            if len(psis_lsv_met2) < 2 or len(psi_met2_rep2[i]) < 2:
                continue  # TODO: check that skipping is not necessary. LSVs with only 1 PSI are wrong..
            if psivalues_met2[1][i][2] not in lsv_types_dict.keys():
                continue
            psi_names_met2[psivalues_met2[1][i][1]] = i

        for psi_name in sorted(psi_names_met1.keys()):
            if psi_name not in set_common_names:
                del psi_names_met1[psi_name]
                continue
            for j, psi_lsv in enumerate(psi_met2_rep1[psi_names_met2[psi_name]]):
                # Try L1 distance
                psi_list1_met2.append(sum(psi_lsv * scripts.src.psi.BINS_CENTER))
                psi_list2_met2.append(sum(psi_met2_rep2[psi_names_met2[psi_name]][j] * scripts.src.psi.BINS_CENTER))

        psi_lists_met2.append(psi_list1_met2)
        psi_lists_met2.append(psi_list2_met2)

    list_l1_expected = []
    for psi_name in sorted(psi_names_met1.keys()):
        debug_dict2[psi_name] = len(psi_met1_rep1[psi_names_met1[psi_name]])
        for j, psi_lsv in enumerate(psi_met1_rep1[psi_names_met1[psi_name]]):
            # Try L1 distance
            # psi_list1_met1.append(sum(psi_lsv*analysis.psi.BINS_CENTER))
            # psi_list2_met1.append(sum(psi_met1_rep2[psi_names_met1[psi_name]][j]*analysis.psi.BINS_CENTER))
            psi_list1_met1.append(psi_lsv)
            psi_list2_met1.append(psi_met1_rep2[psi_names_met1[psi_name]][j])
            # list_l1_expected.append(calculate_l1_expected(psi_lsv, psi_values_lsv2[majiq_psi_names[psi_name]][j]))
            print "%s - MAJIQ:\t%f - %f" % (psi_name, sum(psi_lsv * scripts.src.psi.BINS_CENTER),
                                            sum(psi_met1_rep2[psi_names_met1[psi_name]][j] * scripts.src.psi.BINS_CENTER))
            # print "%s - MISO:\t%s - %s" % (psi_name, str(debug_names_miso_list[psi_name][0][j]), str(debug_names_miso_list[psi_name][1][j]))
            # print "MAJIQ L1 distance:\t%f" % (calculate_l1_expected(psi_lsv, psi_values_lsv2[majiq_psi_names[psi_name]][j]))
        print "-----"

    # print len(debug_dict1), len(debug_dict2), len(psi_list1), len(psi_list2)
    for k in sorted(debug_dict1):
        if k in debug_dict1 and k not in debug_dict2:
            print "This is the guy!! %s" % k
        if debug_dict2[k] - debug_dict1[k]:  # LSVs where the number of junctions in one method differs with the other
            print "Num junctions for LSV %s: MAJIQ - %d MISO: %d" % (k, debug_dict2[k], debug_dict1[k])
            print "MAJIQ LSV and juncs:"
            print "\t", psivalues[1][psi_names_met1[k]]
            print "\t", psivalues[0][0][psi_names_met1[k]]
    # plot_delta_expected_method1Vsmethod2(low_med_high, "MAJIQ", "MISO")

    # plot_PSIs1VsPSIs2(np.array(list_l1_expected), abs(np.array(miso_psis_list[0]) - np.array(miso_psis_list[1])), args.name1, args.name2, "MAJIQ", "MISO", args.plotpath)

    # names_duplicated = [x for x in sorted(debug_dict1) for _ in (0, 1)]
    # names_lsv_where_majiq_lose = np.array(names_duplicated)[abs(np.array(psi_list1) - np.array(psi_list2)) > abs(np.array(miso_psis_list[0]) - np.array(miso_psis_list[1]))][::2]
    # print names_lsv_where_majiq_lose
    # pickle.dump(names_lsv_where_majiq_lose, open('names_lsv_where_majiq_lose.pickle', 'w'))

    # f = figure()
    # ylabel('Number of Junctions')
    # xlabel('Expected PSI')
    # title('MISO and MAJIQ expected PSIs\nHippo1 Vs Hippo2', fontsize=10)
    #
    # hist(np.array(miso_psis_list[1]), bins=40, label="MISO", histtype='step', cumulative=True)
    # hist(np.array(psi_list2), bins=40, label="MAJIQ", histtype='step', cumulative=Truez)
    # legend()
    #
    # _save_or_show(plotpath=None, plotname='MISO')

    print len(psi_list1_met1), len(psi_list2_met1), len(psi_lists_met2[0]), len(psi_lists_met2[1])

    # plot_PSIs1VsPSIs2(calculate_cov(psi_list1_met1, psi_list2_met1), calculate_cov(psi_lists_met2[0], psi_lists_met2[1]), args.name1, args.name2, "MAJIQ", "MISO", args.plotpath)
    plot_PSIs1VsPSIs2(calculate_ead_simple(psi_list1_met1, psi_list2_met1),
                      calculate_ead_simple(psi_lists_met2[0], psi_lists_met2[1]), args.name1, args.name2,
                      "MAJIQ Empirical", "MAJIQ Binomial", args.plotpath)
    # plot_PSIs1VsPSIs2(calculate_dkl(np.array(psi_list1_met1), np.array(psi_list2_met1)), calculate_ead_simple(psi_lists_met2[0], psi_lists_met2[1]), args.name1, args.name2, "MAJIQ", "MISO", args.plotpath)


if __name__ == '__main__':
    main()

