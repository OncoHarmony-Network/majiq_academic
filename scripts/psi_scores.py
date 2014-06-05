from matplotlib import use
use('Agg')

from collections import defaultdict
import sys
import os
import pickle
from scipy.stats import pearsonr
import analysis.psi
import argparse
from math import log

from pylab import *
from matplotlib import rcParams
from scipy.spatial.distance import cityblock

"""
Calculate PSI values
"""
DEBUG = False
TESTBREAK = 50
BINS = linspace(0, 1, num=40)



def plot_PSIs1VsPSIs2(score1, score2, replica1_name, replica2_name, method1, method2, plotpath=None):
    """Compute 2 plots: one for the variance, one for the mean"""
    plotname="%sVs%s\nDelta PSIs for %s - %s" % (method1, method2, replica1_name, replica2_name)

    total_psis = float(len(score1))

    f = figure()
    ylabel(method2)
    xlabel(method1)

    title(plotname, fontsize=10)

    better_in_method1 = np.sum(array(score1) < array(score2))
    better_in_method2 = np.sum(array(score1) > array(score2))

    print "Better in method %s: %.2f%%" % (method1, (better_in_method1/total_psis)*100)
    print "Better in method %s: %.2f%%" % (method2, (better_in_method2/total_psis)*100)
    max_value = max(max(score1), max(score2))

    xlim(0, max_value)
    ylim(0, max_value)
    plot([0, max_value], [0, max_value])
    plot(score1, score2, '.')

    pear, pvalue = pearsonr(score1, score2)
    r_squared = pear**2

    # axarr[name_i].text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, r'$R^2$: %.2f (p-value: %.2E)'%(r_squared, pvalue), fontsize=12, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})
    text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, '%.2f%%' % ((better_in_method1/total_psis)*100), fontsize=16)
    text(abs(max_value)*0.8, max_value-abs(max_value)*0.8, '%.2f%%' % ((better_in_method2/total_psis)*100), fontsize=16)

    _save_or_show(plotpath, plotname)


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
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
    left = log(array(p/q))
    return (left*p).sum(axis=1)

def calculate_l1_expected(p, q):
    return sum(abs(p - q)*analysis.psi.BINS_CENTER)

def calculate_l1(p, q):
    return (abs(p - q)).sum()

def calculate_ead(psi_samples):
    """
    P(S = |PSI_1 - PSI_2|) = sum_psi_1 sum_psi_2 p(psi_1)*p(psi_2)*|psi_1 - psi_2| 

    Expected Absolute Difference = EAD 
    """
    sample1 = psi_samples[:, 0]
    sample2 = psi_samples[:, 1]
    score = 0
    for event_num in xrange(sample1.shape[0]):
        for i in xrange(sample1.shape[1]):
            for j in xrange(sample2.shape[1]):
                psi_val1 = BINS[i]
                psi_val2 = BINS[j]
                score += sample1[event_num][i]*sample2[event_num][j]*abs(psi_val1 - psi_val2)
            
    return score
    #cellcalc = sample1*sample2*abs(sample1 - sample2)
    #return (cellcalc.sum(axis=1)).sum(axis=0)


def main():
    """
    Script for testing MAJIQ against other algorithms for PSIs.

    - Notice that the first file is for MAJIQ results, whereas the second file is reserved for others (initially, MISO).
    - All LSVs 'ways' are equally considered. The order of each way within a LSV should be the same in MAJIQ and MISO.
    """
    parser = argparse.ArgumentParser() 
    parser.add_argument('psivalues_met1', help='Path for psi pickles to evaluate')
    parser.add_argument('psivalues_met2', nargs='+',  help='Path for psi pickles to evaluate')
    parser.add_argument('--name1', default='replica1')
    parser.add_argument('--name2', default='replica2')
    parser.add_argument('--plotpsidist', default=False, help="Plot the PSI distributions for ALL junctions. Slow.")
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--output')
    args = parser.parse_args()

    # For debugging, method1 is MAJIQ, method2 is MISO
    psivalues = pickle.load(open(args.psivalues_met1))
    psi_values_lsv1 = psivalues[0][0]
    psi_values_lsv2 = psivalues[0][1]

    # MAJIQ psi scores
    psi_list1 = []
    psi_list2 = []

    majiq_psi_names = defaultdict()

    lsv_types_dict = {
        's|1e1.1|1e2.1':'SE',
        't|1e1.1|1e2.1':'SE',
        's|1e1.1|1e1.2':'A3SS',
        't|1e1.1|2e1.1':'A3SS',
        't|1e1.1|1e1.2':'A5SS',
        's|1e1.1|2e1.1':'A5SS'
    }


    # Discard LSVs with only one PSI
    for i, psis_lsv in enumerate(psi_values_lsv1):
        if len(psis_lsv) < 2 or len(psi_values_lsv2[i]) < 2:
            continue  # TODO: check that skipping is not necessary. LSVs with only 1 PSI are wrong..
        if psivalues[1][i][2] not in lsv_types_dict.keys():
            continue
        majiq_psi_names[psivalues[1][i][1]] = i
        # print psivalues[1][i]

    miso_psis_list = []

    debug_dict1 = {}
    debug_dict2 = {}

    debug_names_miso_list = defaultdict(list)


    miso_all = []
    for miso_file in args.psivalues_met2:
        miso_psis_dict = defaultdict()
        with open(miso_file, 'r') as miso_res:
            for miso_line in miso_res:
                if miso_line.startswith("event"): continue
                miso_fields = miso_line.split('\t')
                if miso_fields[0] not in majiq_psi_names:
                    continue
                miso_psis_dict[miso_fields[0]] = miso_fields[1]
        miso_all.append(miso_psis_dict)

    miso_common_names = set(miso_all[0].keys()).intersection(miso_all[1].keys())

    for miso_psis_dict in miso_all:
        miso_psis = []
        for psi_name in sorted(majiq_psi_names.keys()):
            if psi_name not in miso_common_names:
                print "%s is not in all MISO replicates" % psi_name
                del majiq_psi_names[psi_name]
                continue
            try:
                miso_psis_values = [float(miso_psi) for miso_psi in miso_psis_dict[psi_name].split(",")]
            except KeyError, e:
                print "LSV %s is in MAJIQ but not in MISO!" % e
                del majiq_psi_names[psi_name]
                continue
                # if len(miso_psis_values) < 2:
                #     del majiq_psi_names[psi_name]
                #     continue

            if len(miso_psis_values) == 1:
                # print miso_psis_values[0], 1 - miso_psis_values[0]
                miso_psis_values.append(1.0 - miso_psis_values[0])
            debug_dict1[psi_name] = len(miso_psis_values)
            miso_psis.extend(miso_psis_values)
            debug_names_miso_list[psi_name].append(miso_psis_values)

        miso_psis_list.append(miso_psis)

    # set_final = set(miso_psis_names[0]).intersection(set(miso_psis_names[1]))
    # print set_final
    # for miso_list in debug_names_miso_list.keys():
    #     for miso_elem in miso_list:
    #         if miso_elem not in set_final:
    #             print "Remove %s" % miso_elem
    #             miso_list.remove(miso_elem)


    list_l1_expected = []
    for psi_name in sorted(majiq_psi_names.keys()):
        debug_dict2[psi_name] = len(psi_values_lsv1[majiq_psi_names[psi_name]])
        # print "%s: " % psi_name
        sys.stdout.flush()
        for j, psi_lsv in enumerate(psi_values_lsv1[majiq_psi_names[psi_name]]):

            # Try L1 distance
            psi_list1.append(sum(psi_lsv*analysis.psi.BINS_CENTER))
            psi_list2.append(sum(psi_values_lsv2[majiq_psi_names[psi_name]][j]*analysis.psi.BINS_CENTER))
            list_l1_expected.append(calculate_l1_expected(psi_lsv, psi_values_lsv2[majiq_psi_names[psi_name]][j]))
            # print "MAJIQ:\t%f - %f" % (sum(psi_lsv*analysis.psi.BINS_CENTER), sum(psi_values_lsv2[majiq_psi_names[psi_name]][j]*analysis.psi.BINS_CENTER))
            # print "MAJIQ L1 distance:\t%f" % (calculate_l1_expected(psi_lsv, psi_values_lsv2[majiq_psi_names[psi_name]][j]))
            # print "MISO:\t%s - %s" % (str(debug_names_miso_list[psi_name][0][j]), str(debug_names_miso_list[psi_name][1][j]))


    print len(debug_dict1), len(debug_dict2)
    for k in sorted(debug_dict1):
        if k in debug_dict1 and k not in debug_dict2:
            print "This is the guy!! %s" % k
        # if debug_dict2[k] - debug_dict1[k]:  # LSVs where the number of junctions in one method differs with the other
        #     print "Num junctions for LSV %s: MAJIQ - %d MISO: %d" % (k, debug_dict2[k], debug_dict1[k])
        #     print "MAJIQ LSV and juncs:"
        #     print "\t", psivalues[1][majiq_psi_names[k]]
        #     print "\t", psivalues[0][0][majiq_psi_names[k]]
    # plot_PSIs1VsPSIs2(np.array(list_l1_expected), abs(np.array(miso_psis_list[0]) - np.array(miso_psis_list[1])), args.name1, args.name2, "MAJIQ", "MISO", args.plotpath)
    plot_PSIs1VsPSIs2(abs(np.array(psi_list1) - np.array(psi_list2)), abs(np.array(miso_psis_list[0]) - np.array(miso_psis_list[1])), args.name1, args.name2, "MAJIQ", "MISO", args.plotpath)


if __name__ == '__main__':
    main()

