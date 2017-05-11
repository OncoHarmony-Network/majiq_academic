from matplotlib import use
from voila.vlsv import collapse_matrix

use('Agg')
from matplotlib import rcParams
rcParams.update({'font.size': 10})
from collections import defaultdict
import cPickle as pickle
import numpy as np
import argparse
import os
import matplotlib.pyplot as pyplot
import scripts.utils
import colorbrewer as cb
from scipy.stats import pearsonr


__author__ = 'abarrera'


def get_color(tissue):
    colors_dict = {'rest_stim': '#%02x%02x%02x' % cb.Paired[10][-1],
               'cer_liv': '#%02x%02x%02x' % cb.Paired[10][-3]
               }
    return colors_dict[tissue]


def scatterplot_rtpcr_majiq_miso(rt_pcr, majiq, miso, cov, plotpath, pcr_majiq_extra=None, pcr_miso_extra=None, majiq_extra=None, miso_extra=None):
    #figure out how many groups of events exist

    # majiq_rest_yerr = [var_expected_psi(dis) for dis in rt_pcr_majiq[0]]

    fig, axx = pyplot.subplots(2, 3, sharex=True, sharey=True, figsize=[9, 6], dpi=300)
    fig.suptitle("DeltaPSI comparison: RT-PCR Vs MAJIQ; RT-PCR Vs MISO")


    # Hogenesch data
    # pcr_majiq_extra_all = [a for b in pcr_majiq_extra for a in b]
    # pcr_miso_extra_all = [a for b in pcr_miso_extra for a in b]
    # majiq_extra_all = [a for b in majiq_extra for a in b]
    # miso_extra_all = [a for b in miso_extra for a in b]


    diagonal = np.linspace(-1, 1, num=len(rt_pcr))
    axx[0][0].plot(diagonal, diagonal, '--', color="#cccccc")

    fit = np.polyfit(np.append(majiq, majiq_extra), np.append(rt_pcr, np.array(pcr_majiq_extra)), 1)
    fit_fn = np.poly1d(fit)  # fit_fn is now a function which takes in x and returns an estimate for y

    axx[0][0].plot(np.append(majiq, majiq_extra), fit_fn(np.append(majiq, majiq_extra)), '--k')
    axx[0][0].plot(majiq, rt_pcr, '.', color=get_color('rest_stim'), label='Unstim Vs Stim')
    axx[0][0].plot(majiq_extra, pcr_majiq_extra, 'd', color=get_color('cer_liv'), label='Hogenesch Cer Vs Liv')
    pearson_r = pearsonr(np.append(majiq, majiq_extra), np.append(rt_pcr, pcr_majiq_extra))[0]
    axx[0][0].text(-.9, .7, 'R=%.2f' % pearson_r, fontsize=14)


    axx[0][0].set_xlabel('MAJIQ')
    axx[0][0].set_ylabel('RT-PCR')
    axx[0][0].set_title('All (N=%d)' % (len(majiq) + len(majiq_extra)))
    axx[0][0].set_xlim([-1,1])
    axx[0][0].set_ylim([-1,1])
    # axx[0][0].legend(loc=2, fontsize=8)

    diagonal = np.linspace(-1, 1, num=len(rt_pcr))
    fit = np.polyfit(np.append(miso, miso_extra), np.append(rt_pcr, np.array(pcr_miso_extra)), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[1][0].plot(np.append(miso, miso_extra), fit_fn(np.append(miso, miso_extra)), '--k')
    axx[1][0].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][0].plot(miso, rt_pcr, '.', color=get_color('rest_stim'), label='Unstim Vs Stim')
    axx[1][0].plot(miso_extra, pcr_miso_extra, 'd', color=get_color('cer_liv'), label='Hogenesch Cer Vs Liv')
    pearson_r = pearsonr(np.append(miso, miso_extra), np.append(rt_pcr, pcr_miso_extra))[0]
    axx[1][0].text(-.9, .7, 'R=%.2f' % pearson_r, fontsize=14)


    axx[1][0].set_xlabel('MISO')
    axx[1][0].set_ylabel('RT-PCR')
    axx[1][0].set_title('All (N=%d)' % (len(miso) + len(miso_extra)))
    axx[1][0].set_xlim([-1,1])
    axx[1][0].set_ylim([-1,1])
    # axx[1][0].legend(loc=2, fontsize=8)


    diagonal = np.linspace(-1, 1, num=len(rt_pcr))
    fit = np.polyfit(majiq, rt_pcr, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[0][1].plot(majiq, fit_fn(majiq), '--k')
    axx[0][1].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][1].plot(majiq, rt_pcr, '.', color=get_color('rest_stim'), label='Unstim Vs Stim')
    pearson_r = pearsonr(majiq, rt_pcr)[0]
    axx[0][1].text(-.9, .7, 'R=%.2f' % pearson_r, fontsize=14)

    axx[0][1].set_xlabel('MAJIQ')
    axx[0][1].set_ylabel('RT-PCR')
    axx[0][1].set_title('Unstim Vs Stim (N=%d)' % len(majiq))
    axx[0][1].set_xlim([-1, 1])
    axx[0][1].set_ylim([-1, 1])
    # axx[0][0].legend(loc=4)

    fit = np.polyfit(majiq_extra, pcr_majiq_extra, 1)
    fit_fn = np.poly1d(fit)
    axx[0][2].plot(majiq_extra, fit_fn(majiq_extra), '--k')
    axx[0][2].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][2].plot(majiq_extra, pcr_majiq_extra, 'd', color=get_color('cer_liv'), label='Hogenesch Cer Vs Liv')
    pearson_r = pearsonr(majiq_extra, pcr_majiq_extra)[0]
    axx[0][2].text(-.9, .7, 'R=%.2f' % pearson_r, fontsize=14)


    axx[0][2].set_xlabel('MAJIQ')
    axx[0][2].set_ylabel('RT-PCR')
    axx[0][2].set_title('Cerebellum Vs Liver (N=%d)' % len(majiq_extra))
    axx[0][2].set_ylim([-1,1])
    # axx[0][1].legend(loc=4)


    fit = np.polyfit(miso, rt_pcr, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[1][1].plot(miso, fit_fn(miso), '--k')
    axx[1][1].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][1].plot(miso, rt_pcr, '.', color=get_color('rest_stim'), label='Unstim Vs Stim')
    pearson_r = pearsonr(miso, rt_pcr)[0]
    axx[1][1].text(-.9, .7, 'R=%.2f' % pearson_r, fontsize=14)

    axx[1][1].set_xlabel('MISO')
    axx[1][1].set_ylabel('RT-PCR')
    axx[1][1].set_title('Unstim Vs Stim (N=%d)' % len(miso))
    axx[1][1].set_xlim([-1, 1])
    axx[1][1].set_ylim([-1, 1])
    # axx[1][0].legend(loc=4)


    fit = np.polyfit(miso_extra, pcr_miso_extra, 1)
    fit_fn = np.poly1d(fit)
    axx[1][2].plot(miso_extra, fit_fn(miso_extra), '--k')
    axx[1][2].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][2].plot(miso_extra, pcr_miso_extra, 'd', color=get_color('cer_liv'), label='Hogenesch Cer Vs Liv')
    pearson_r = pearsonr(miso_extra, pcr_miso_extra)[0]
    axx[1][2].text(-.9, .7, 'R=%.2f' % pearson_r, fontsize=14)



    axx[1][2].set_xlabel('MISO')
    axx[1][2].set_ylabel('RT-PCR')
    axx[1][2].set_title('Cerebellum Vs Liver (N=%d)' % len(miso_extra))
    axx[1][2].set_ylim([-1,1])
    # axx[0][1].legend(loc=4)

    scripts.utils.save_or_show(plotpath, "dpsi_rtpcr_majiq_miso", exten='pdf')


def scatterplot_rtpcr_simple(rt_pcr, method_epsis, cov,  plotpath, rt_pcr_extra=None, method_extra=None, cov_extra=None, plotname='psi_rtpcr', met_name='MAJIQ'):
    #figure out how many groups of events exist

    # majiq_rest_yerr = [var_expected_psi(dis) for dis in rt_pcr_majiq[0]]

    # Changing events
    CHANGE = .2
    rt_pcr_chg_mask = abs(rt_pcr) > CHANGE
    rt_pcr_extra_chg_mask = abs(rt_pcr_extra) > CHANGE
    method_chg_mask = abs(method_epsis) > CHANGE
    method_extra_chg_mask = abs(method_extra) > CHANGE

    kristen_false_pos = np.count_nonzero(~rt_pcr_chg_mask & method_chg_mask)
    hogenesch_false_pos = np.count_nonzero(~rt_pcr_extra_chg_mask & method_extra_chg_mask)
    print "%s [Kristen] false positive: %.2f%%" % (met_name, 100.*kristen_false_pos/rt_pcr_chg_mask.size)
    print "%s [Hogenesch] false positive: %.2f%%" % (met_name, 100.*hogenesch_false_pos/rt_pcr_extra_chg_mask.size)
    print "%s [all] false positive: %.2f%%" % (met_name, 100.*(kristen_false_pos + hogenesch_false_pos)/(rt_pcr_chg_mask.size + rt_pcr_extra_chg_mask.size))

    kristen_false_neg = np.count_nonzero(rt_pcr_chg_mask & ~method_chg_mask)
    hogenesch_false_neg = np.count_nonzero(rt_pcr_extra_chg_mask & ~method_extra_chg_mask)
    print "%s [Kristen] false negative: %.2f%%" % (met_name, 100.*kristen_false_neg / rt_pcr_chg_mask.size)
    print "%s [Hogenesch] false negative: %.2f%%" % (met_name, 100.*hogenesch_false_neg / rt_pcr_extra_chg_mask.size)
    print "%s [all] false negative: %.2f%%" % (met_name, 100.*(kristen_false_neg + hogenesch_false_neg) / (rt_pcr_chg_mask.size + rt_pcr_extra_chg_mask.size))

    # rt_pcr = rt_pcr[rt_pcr_chg_mask]
    # method_epsis = method_epsis[rt_pcr_chg_mask]
    # rt_pcr_extra = rt_pcr_extra[rt_pcr_extra_chg_mask]
    # method_extra = method_extra[rt_pcr_extra_chg_mask]

    # Coverage breakdown (used for coloring)
    cov_thres = [0, 15, 40, 100000]
    cov_colors = cb.Accent[3]

    fig = pyplot.figure(figsize=[6, 6], dpi=300)
    fig.suptitle("RT-PCR PSI comparison")


    diagonal = np.linspace(-1, 1, num=len(rt_pcr))
    fit = np.polyfit(np.concatenate((method_epsis[rt_pcr_chg_mask], method_extra[rt_pcr_extra_chg_mask])), np.concatenate((rt_pcr[rt_pcr_chg_mask], rt_pcr_extra[rt_pcr_extra_chg_mask])), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    pyplot.plot(np.concatenate((method_epsis, method_extra)), fit_fn(np.concatenate((method_epsis, method_extra))), '--k')
    pyplot.plot(diagonal, diagonal, '--', color="#cccccc")
    pearson_r = pearsonr(np.append(method_epsis[rt_pcr_chg_mask], method_extra[rt_pcr_extra_chg_mask]), np.append(rt_pcr[rt_pcr_chg_mask], rt_pcr_extra[rt_pcr_extra_chg_mask]))[0]
    pyplot.text(-.9, .8, 'R=%.2f%%' % pearson_r, fontsize=14)

    for cov_idx, cov_th in enumerate(cov_thres[1:]):
        cov_mask = (cov > cov_thres[cov_idx]) & (cov <= cov_th)
        cov_extra_mask = (cov_extra > cov_thres[cov_idx]) & (cov_extra <= cov_th)
        pyplot.plot(method_epsis[cov_mask & rt_pcr_chg_mask], rt_pcr[cov_mask & rt_pcr_chg_mask], '.',
                    color='#%02x%02x%02x' % cb.Paired[10][-1],
                    label='Unst Vs Stim - Cov (%d, %d],N=%d' % (cov_thres[cov_idx], cov_th, np.count_nonzero(cov_mask & rt_pcr_chg_mask)))
        pyplot.plot(method_extra[cov_extra_mask & rt_pcr_extra_chg_mask], rt_pcr_extra[cov_extra_mask & rt_pcr_extra_chg_mask], 'd',
                    color='#%02x%02x%02x' % cb.Paired[10][-3],
                    label='Cer Vs Liv - Cov (%d, %d],N=%d' % (cov_thres[cov_idx], cov_th, np.count_nonzero(cov_extra_mask & rt_pcr_extra_chg_mask)))

    pyplot.xlabel(met_name)
    pyplot.ylabel('RT-PCR')
    pyplot.title('All (N=%d)' % (np.count_nonzero(rt_pcr_chg_mask) + np.count_nonzero(rt_pcr_extra_chg_mask)))
    pyplot.xlim([-1, 1])
    pyplot.ylim([-1, 1])
    pyplot.legend(loc='lower right', fontsize=7)
    scripts.utils.save_or_show(plotpath, plotname, exten='pdf')


def barchart_expected(expected_psis, plotpath, mfile):
    pyplot.figure()
    pyplot.hist(expected_psis, bins=40)
    name = os.path.basename(mfile)
    scripts.utils.save_or_show(plotpath, name +'_expected_dist')




def get_expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1+1./len(bins), 1., 2./len(bins)))


def expected_psi(bins):
    bins = np.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * np.arange((1.*step) / 2, 1, step)
    return np.sum(projection_prod)


def avg_expected_psi(bins_list):
    return np.mean([expected_psi(bins) for bins in bins_list])


def get_variance(bins, mean):
    """Compute the variance = E[X^2] - (E[X])^2"""
    bins = np.array(bins)
    step_bins = 1 / bins.size
    projection_prod = bins * np.arange(step_bins / 2, 1, step_bins)**2
    return np.sum(projection_prod) - mean**2

def var_expected_psi(bins_list):
    return np.mean([expected_psi(bins) for bins in bins_list])


def expected_avg_psi(bins_list):
    return expected_psi(np.mean(bins_list, axis=0))


def expected_delta_psi(bins_list, rtpcr_psi):
    step = 1.0 / 40
    bins_index = np.arange(step / 2, 1, step)
    return np.mean([np.sum([bins[i]*abs(psi - rtpcr_psi) for i, psi in enumerate(bins_index)]) for bins in bins_list])


def parse_rtpcr_results(pcr_file, names_pcr2majiq_dict):
    pcr_rest_stim = defaultdict(list)
    with open(pcr_file, 'r') as pcr_file:
        for i, pcr_line in enumerate(pcr_file):
            if i<1 or pcr_line.startswith('#'): continue  # headers
            pcr_fields = pcr_line.rstrip().split('\t')
            # print pcr_fields
            if pcr_fields[0] not in names_pcr2majiq_dict.keys():
                continue  # If the event is not in the list of events selected, skip it
            try:
                pcr_rest_stim[pcr_fields[0]].append((float(pcr_fields[1]) - float(pcr_fields[2]))/100)
            except IndexError:
                print "%s value not found, assigned 0..." % pcr_fields[0]
                pcr_rest_stim[pcr_fields[0]].append(0)
            except ValueError, e:
                print e.message
                print "Event PSI could not be converted into float, wrong parsing??"
                pcr_rest_stim[pcr_fields[0]].append(0)
    print "Number of events found in RT-PCR file: %d" % len(pcr_rest_stim)
    return pcr_rest_stim


def parse_majiq_results(files_majiq, names_junc_majiq):
    majiq_dict = defaultdict(list)
    for mfile in files_majiq:
        with open(mfile) as mfile_open:
            mpickle = pickle.load(mfile_open)
            for i, lsv in enumerate(mpickle.lsvs):
                if lsv.get_id() not in names_junc_majiq.keys():
                    continue # If the event is not in the list of events selected, skip it

                jcoords_lsv = [jvisual.coords for jvisual in lsv.lsv_graphic.get_junctions()]
                # Find all junctions from that LSV that are included
                for jj, jcoords in enumerate(names_junc_majiq[lsv.get_id()]):
                    j = jcoords_lsv.index(jcoords[:2])
                    try:
                        majiq_dict[lsv.get_id()+"#"+str(jcoords[2])].append(collapse_matrix(np.array(lsv.get_bins()[j])))
                    except IndexError:
                        majiq_dict[lsv.get_id()+"#"+str(jcoords[2])].append(collapse_matrix(np.array(lsv.get_bins()[0]))[::-1])

    print "Number of events found in MAJIQ %s: %d" % ("; ".join([mf for mf in files_majiq]), len(majiq_dict.keys()))
    return majiq_dict


def parse_miso_results(files_miso, names_junc_majiq):

    # event_name	miso_posterior_mean	ci_low	ci_high	isoforms	counts	assigned_counts	chrom	strand	mRNA_starts	mRNA_ends
    # ENST00000301332:145697325-145697429:target	0.98	0.91	1.00	'ENST00000301332:145697325-145697429:target.0.ex_ENST00000301332:145697325-145697429:target.0.lsv','ENST00000301332:145697325-145697429:target.1.ex_ENST00000301332:145697325-145697429:target.1.lsv'	(0,0):36,(1,0):37	0:37	chr8	+	145694873,145695965	145697429,145697429

    miso_dict = defaultdict(list)
    files_miso = scripts.utils.list_files_or_dir(files_miso)

    for mfile in files_miso:
        with open(mfile) as mfile_open:
            for line in mfile_open:
                miso_fields = line.rstrip().split()
                if miso_fields[0] not in names_junc_majiq.keys(): continue # If the event is not in the list of events selected, skip it

                lsv_name = miso_fields[0]
                miso_psis = miso_fields[1]
                # Find all junctions from that LSV that are included

                for i, miso_psi in enumerate(miso_psis.split(',')):
                    miso_dict[lsv_name+"#"+str(i)].append(float(miso_psi))

    print "Number of events found in MISO %s: %d" % ("; ".join([mf for mf in files_miso]), len(miso_dict.keys()))
    return miso_dict

def plot_dpsi_rtpcr_majiq_miso(rt_pcr_majiq, rt_pcr_miso, majiq, miso, plotpath):

    fig, (ax1, ax2) = pyplot.subplots(1, 2, sharey=True, figsize=[12, 6], dpi=300)
    pyplot.suptitle("Delta PSI comparison: RT-PCR Vs (MAJIQ N=%d; MISO N=%d)" % (len(majiq), len(miso)))
    diagonal = np.linspace(-1, 1, num=len(rt_pcr_majiq))

    # All
    ax1.scatter(majiq, rt_pcr_majiq, c=cb.Blues[3][-1], alpha=.5, s=50, label='MAJIQ')
    ax2.scatter(miso, rt_pcr_miso, c=cb.Reds[3][-1], alpha=.5, s=50, label='MISO')
    ax1.set_ylim([-1, 1])
    ax2.set_ylim([-1, 1])
    ax1.set_xlim([-1, 1])
    ax2.set_xlim([-1, 1])
    ax1.set_ylabel('RT-PCR')
    ax1.set_xlabel('MAJIQ')
    ax2.set_xlabel('MISO')
    fit = np.polyfit(majiq, rt_pcr_majiq, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    ax1.plot(majiq, fit_fn(majiq), '--k')
    fit = np.polyfit(miso, rt_pcr_miso, 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    ax2.plot(miso, fit_fn(miso), '--k')

    ax2.plot(majiq, fit_fn(majiq), '--k')
    ax1.plot(diagonal, diagonal, '--', color="#cccccc")
    ax2.plot(diagonal, diagonal, '--', color="#cccccc")

    print majiq, rt_pcr_majiq, miso, rt_pcr_miso
    scripts.utils.save_or_show(plotpath, "dpsi_rtpcr_majiq_miso", exten='pdf')


def load_coverage(lcov_met1, lcov_met2, clean_reads_files, dataset='ho'):
    misannotated = {}
    misannotated['ENSMUSG00000031996:31153414-31153499:source'] = 'ENSMUSG00000031996:31153414-31153501:source'
    misannotated['ENSMUSG00000024830:4161068-4161215:target'] = 'ENSMUSG00000024830:4161064-4161215:target'
    cov1_ll = [[0, 0] for ii in xrange(len(lcov_met1))]
    cov2_ll = [[0, 0] for ii in xrange(len(lcov_met2))]
    for crf_ix, clean_reads_file in enumerate(clean_reads_files):
        with open(clean_reads_file) as clean_reads_f:
            for clean_reads in pickle.load(clean_reads_f):
                if dataset == 'ho' and clean_reads[0] in misannotated.keys():
                    clean_reads[0] = misannotated[clean_reads[0]]
                try:
                    # lcov_met1[lcov_met1.index(clean_reads[0])] = clean_reads[1].sum()
                    cov1_ll[lcov_met1.index(clean_reads[0])][crf_ix] = clean_reads[1].sum()
                except ValueError:
                    pass
                try:
                    # lcov_met2[lcov_met2.index(clean_reads[0])] = clean_reads[1].sum()
                    cov2_ll[lcov_met2.index(clean_reads[0])][crf_ix] = clean_reads[1].sum()
                except ValueError:
                    pass
    lcov_met1 = np.average(cov1_ll, axis=1)
    try:
        lcov_met2 = np.average(cov2_ll, axis=1)
    except IndexError:
        pass
    return lcov_met1, lcov_met2


def main():
    parser = argparse.ArgumentParser(description="Compare PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--old_majiq-rest-stim", required=True, dest='majiq_rest_stim', nargs=1, help='MAJIQ delta PSI predictions for Unstim Vs Stim RNA-Seq data.')
    parser.add_argument("--miso-rest", dest='miso_rest', nargs='*', help='MISO PSI predictions for resting RNA-Seq data.')
    parser.add_argument("--miso-stim", dest='miso_stim', nargs='*', help='MISO PSI predictions for stimulated RNA-Seq data.')
    parser.add_argument("--names-map-file", required=True, dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')
    parser.add_argument("--clean-reads", dest='clean_reads', nargs=2, help='Clean reads from MAJIQ quantifications (one file per condition).')

    parser.add_argument("--old_majiq-extra", dest='majiq_extra', nargs=1, help='MAJIQ PSI predictions for extra conditions.')
    parser.add_argument("--miso-extra", dest='miso_extra', nargs='+', help='MISO PSI predictions for extra conditions.')
    parser.add_argument("--pcr-extra", dest='pcr_extra', help='RT-PCR validations for extra conditions.')
    parser.add_argument("--clean-reads-extra", dest='clean_reads_extra', nargs=2, help='Clean reads from MAJIQ quantifications (one file per condition).')

    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    ## ------------------------ ##
    ## Hogenesch data
    ## ------------------------ ##
    lrtpcr_majiq_extra = []
    lrtpcr_miso_extra = []
    lmajiq_extra = []
    lmiso_extra = []
    lcov_majiq_extra = []
    lcov_miso_extra = []

    if args.majiq_extra:
        # 1. Read PCR results
        rtpcr_extra = defaultdict(list)
        with open(args.pcr_extra) as pcr_extra:
            for pcr_elem in pcr_extra:
                if pcr_elem.startswith("#"): continue
                pcr_elems=pcr_elem.rstrip().split()
                rtpcr_extra[pcr_elems[0].split("#")[0]] = [float(pcr_elems[2])/100, np.nan, [], [], float(pcr_elems[6])/100, float(pcr_elems[5])/100, pcr_elems[-4]]


        # Read Majiq results for the elements in PCR
        djunc_selected = {}

        rtpcr_extra['ENSMUSG00000031996:31153414-31153501:source'] = rtpcr_extra['ENSMUSG00000031996:31153414-31153499:source']
        del rtpcr_extra['ENSMUSG00000031996:31153414-31153499:source']
        rtpcr_extra['ENSMUSG00000024830:4161064-4161215:target'] = rtpcr_extra['ENSMUSG00000024830:4161068-4161215:target']
        del rtpcr_extra['ENSMUSG00000024830:4161068-4161215:target']

        for mfile in args.majiq_extra:
            majiq_found = []
            with open(mfile) as mfile_open:
                mpickle = pickle.load(mfile_open)
                for j, lsv in enumerate(mpickle.get_lsvs()):
                    if lsv.get_id() in rtpcr_extra:
                        majiq_found.append(lsv.get_id())
                        for lway_aux, lway_junc in enumerate(lsv.lsv_graphic.get_junctions()):
                            if set(lway_junc.get_coords()) == set([int(aa) for aa in rtpcr_extra[lsv.get_id()][-1].split('-')]):
                                print "Junction selected: %d" % lway_aux
                                lsv_aux = lsv.get_id() + "#%d"%lway_aux
                                djunc_selected[lsv.get_id()] = lway_aux
                                try:
                                    if len(lsv.get_bins())>1:
                                        print "Complex LSV; %s" % lsv_aux
                                    else:
                                        print "Binary LSV; %s" % lsv_aux
                                    rtpcr_extra[lsv.get_id()][1] = -get_expected_dpsi(collapse_matrix(lsv.get_bins()[lway_aux]))
                                except IndexError:
                                    rtpcr_extra[lsv.get_id()][1] = get_expected_dpsi(collapse_matrix(lsv.get_bins()[0]))
            print "Missing LSVs in %s:" % mfile
            print '\n'.join([kk for kk in rtpcr_extra if kk not in majiq_found])

        # Read MISO results for the elements in PCR
        for cn, miso_dir in enumerate(args.miso_extra):
            files_miso = scripts.utils.list_files_or_dir([miso_dir], containing='miso')
            for mfile in files_miso:
                with open(mfile) as mfile_open:
                    for line in mfile_open:
                        miso_fields = line.rstrip().split()
                        lsv_name = miso_fields[0]
                        miso_psis = miso_fields[1].split(',')
                        miso_starts = miso_fields[-2].split(',')
                        miso_ends = miso_fields[-1].split(',')
                        if lsv_name in rtpcr_extra:
                            # Find all junctions from that LSV that are included
                            try:
                                rtpcr_extra[lsv_name][2+cn].append(float(miso_psis[djunc_selected[lsv_name]]))
                            except IndexError:
                                rtpcr_extra[lsv_name][2+cn].append(1-float(miso_psis[0]))
                            except KeyError:
                                print "[WARNING] :: %s in MISO, but not in MAJIQ, skipped given the impossibility of determining which junction it is." % lsv_name

        with open('dpsi_hogenesch.txt', 'w') as psi_txt:
            headers=['LSV ID', 'RT-PCR', 'Majiq', 'Miso Cerebellum', 'Miso Liver', 'RT-PCR Cerebellum Avg. STD', 'RT-PCR Liver Avg. STD']
            psi_txt.write('\t'.join(headers))
            psi_txt.write('\n')
            for ilsv, vals in rtpcr_extra.iteritems():
                line = [ilsv]
                line.extend([repr(vv) for vv in vals[:-1]])
                psi_txt.write("\t".join(line))
                psi_txt.write('\n')
                if not np.isnan(vals[1]):
                    lrtpcr_majiq_extra.append(vals[0])
                    lmajiq_extra.append(vals[1])
                    lcov_majiq_extra.append(ilsv)
                if not np.isnan(np.mean(vals[2])) and not np.isnan(np.mean(vals[3])):
                    lrtpcr_miso_extra.append(vals[0])
                    lmiso_extra.append(np.mean(vals[2]) - np.mean(vals[3]))
                    lcov_miso_extra.append(ilsv)
        # Find out coverage (clean reads) for each LSV junction)
        lcov_majiq_extra, lcov_miso_extra = load_coverage(lcov_majiq_extra, lcov_miso_extra, args.clean_reads_extra)

    ## ------------------------ ##
    ## Kristen's data
    ## ------------------------ ##
    names_pcr2majiq_dict = {}  # key RT-PCR, value MAJIQ
    names_majiq2pcr_dict = {}
    names_junc_majiq = defaultdict(list)
    gene_names_counts = defaultdict(lambda: 0)
    with open(args.names_map_file) as names_map:
        for name_map in names_map:
            # MAJIQ
            mapped_name = name_map.rstrip().split()
            majiq_name, junc_idx_old = mapped_name[1].split('#')
            names_junc_majiq[majiq_name].append([int(jcoord) for jcoord in mapped_name[-1].split('-')] + [int(junc_idx_old)])

            # RT-PCR
            names_pcr2majiq_dict[mapped_name[0]] = mapped_name[1]
            names_majiq2pcr_dict[mapped_name[1]] = mapped_name[0]
            gene_names_counts[mapped_name[0]] += 1

    print "Number of events in names-map-file RT-PCR: %d" % len(names_pcr2majiq_dict.keys())
    print "Number of events in names-map-file MAJIQ: %d" % len(names_junc_majiq.keys())

    # Parse RT-PCR results
    pcr_rest_stim = parse_rtpcr_results(args.pcr, names_pcr2majiq_dict)

    # Process MAJIQ files for resting RNA-Seq data
    majiq_rest_stim_dict = parse_majiq_results(args.majiq_rest_stim, names_junc_majiq)

    # Process MAJIQ files for resting RNA-Seq data
    miso_rest_dict = parse_miso_results(args.miso_rest, names_junc_majiq)

    # Process MAJIQ files for stimuli RNA-Seq data
    miso_stim_dict = parse_miso_results(args.miso_stim, names_junc_majiq)

    ## Intersect names from RT-PCR and MAJIQ
    common_names_set = set([names_pcr2majiq_dict[k] for k in pcr_rest_stim.keys()]).intersection(set(majiq_rest_stim_dict.keys()))
    print "Common names afer intersection with MAJIQ: %d" % len(common_names_set)
    common_names_set = common_names_set.intersection(set(miso_rest_dict.keys())).intersection(set(miso_stim_dict.keys()))
    print "Common names afer intersection with MISO: %d" % len(common_names_set)

    kr_rt_pcr_dpsi= []
    kr_majiq_dpsi = []
    kr_miso_dpsi = []
    kr_cov = []

    with open('dpsi_kristens.txt', 'w') as psi_txt:
        headers=['LSV ID', 'RT-PCR', 'Majiq', 'Miso Stim.', 'Miso Unstim.', 'RT-PCR Stim. Avg. STD', 'RT-PCR Liver Avg. STD']
        psi_txt.write('\t'.join(headers))
        psi_txt.write('\n')

        for common_name in common_names_set:
            for name_majiq, name_pcr in names_majiq2pcr_dict.iteritems():
                if names_majiq2pcr_dict[common_name] == name_pcr:
                    name = name_majiq
                    psi_txt_line = [name]

                    # For Majiq, compute mean over Expected PSIs
                    majiq_rest_stim_stat = get_expected_dpsi(majiq_rest_stim_dict[name][0])

                    # For MISO, compute mean
                    miso_rest_stat = np.mean(miso_rest_dict[name])
                    miso_stim_stat = np.mean(miso_stim_dict[name])


                    print "%s - %s" % (name, names_majiq2pcr_dict[name])
                    print "---- RT-PCR ----"

                    kr_rt_pcr_dpsi.append(pcr_rest_stim[names_majiq2pcr_dict[name]][0])
                    kr_majiq_dpsi.append(majiq_rest_stim_stat)
                    kr_miso_dpsi.append(miso_rest_stat-miso_stim_stat)
                    kr_cov.append(name.split("#")[0])

                    psi_txt_line.append(str(np.mean(kr_rt_pcr_dpsi[-1])))
                    psi_txt_line.append(str(kr_majiq_dpsi[-1]))
                    psi_txt_line.append(str(miso_rest_dict[name]))
                    psi_txt_line.append(str(miso_stim_dict[name]))
                    psi_txt_line.append(str(np.std(kr_rt_pcr_dpsi[-1])))
                    psi_txt.write("\t".join(psi_txt_line))
                    psi_txt.write('\n')

                    print "---- RT-PCR ----"
                    print "[Resting Vs Stimulated - RTPCR]:\t%f" % kr_rt_pcr_dpsi[-1]

                    print "---- MAJIQ ----"
                    print "[Resting Vs Stimulated - MAJIQ]:\t%f" % (kr_majiq_dpsi[-1])

                    print "---- MISO -----"
                    print "[Resting Vs Stimulated - MISO]:\t%f" % (kr_miso_dpsi[-1])

    kr_cov, foo = load_coverage(kr_cov, [], args.clean_reads, dataset='kr')

    scatterplot_rtpcr_majiq_miso(np.array(kr_rt_pcr_dpsi), np.array(kr_majiq_dpsi), np.array(kr_miso_dpsi), np.array(kr_cov), args.plotpath,
                                 pcr_majiq_extra=np.array(lrtpcr_majiq_extra), pcr_miso_extra=np.array(lrtpcr_miso_extra), majiq_extra=np.array(lmajiq_extra),
                                 miso_extra=np.array(lmiso_extra))
    scatterplot_rtpcr_simple(np.array(kr_rt_pcr_dpsi), np.array(kr_majiq_dpsi), np.array(kr_cov), args.plotpath, rt_pcr_extra=np.array(lrtpcr_majiq_extra),
                             method_extra=np.array(lmajiq_extra), cov_extra=np.array(lcov_majiq_extra), plotname='dpsi_majiq_only', met_name='MAJIQ')
    scatterplot_rtpcr_simple(np.array(kr_rt_pcr_dpsi), np.array(kr_miso_dpsi), np.array(kr_cov), args.plotpath, rt_pcr_extra=np.array(lrtpcr_miso_extra),
                             method_extra=np.array(lmiso_extra), cov_extra=np.array(lcov_miso_extra), plotname='dpsi_miso_only', met_name='MISO')

if __name__ == '__main__':
    main()