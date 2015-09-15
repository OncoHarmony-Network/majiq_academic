from matplotlib import use
use('Agg')
from collections import defaultdict
import cPickle as pickle
import numpy as np
import argparse
import os
import matplotlib.pyplot as pyplot
import scripts.utils
import colorbrewer as cb
from scipy.stats.stats import pearsonr
from matplotlib import rcParams
rcParams.update({'font.size': 10})


__author__ = 'abarrera'


def get_color(tissue):
    colors_dict = {'rest': '#%02x%02x%02x' % cb.Paired[10][-1],
               'stim': '#%02x%02x%02x' % cb.Paired[10][-2],
               'cer': '#%02x%02x%02x' % cb.Paired[10][-3],
               'liv': '#%02x%02x%02x' % cb.Paired[10][-4]
               }
    return colors_dict[tissue]


def scatterplot_rtpcr_majiq_miso(rt_pcr_majiq, rt_pcr_miso, majiq, miso, plotpath, pcr_majiq_extra=None, pcr_miso_extra=None, majiq_extra=None, miso_extra=None):
    #figure out how many groups of events exist

    # majiq_rest_yerr = [var_expected_psi(dis) for dis in rt_pcr_majiq[0]]

    fig, axx = pyplot.subplots(2, 5, sharex=True, sharey=True, figsize=[15, 6], dpi=300)
    fig.suptitle("PSI comparison: RT-PCR Vs MAJIQ; RT-PCR Vs MISO")

    rt_pcr_majiq_all = [a for b in rt_pcr_majiq for a in b]
    rt_pcr_miso_all = [a for b in rt_pcr_miso for a in b]
    majiq_all = [a for b in majiq for a in b]
    miso_all = [a for b in miso for a in b]

    pcr_majiq_extra_all = [a for b in pcr_majiq_extra for a in b]
    pcr_miso_extra_all = [a for b in pcr_miso_extra for a in b]
    majiq_extra_all = [a for b in majiq_extra for a in b]
    miso_extra_all = [a for b in miso_extra for a in b]


    diagonal = np.linspace(0, 1, num=len(rt_pcr_majiq_all))
    fit = np.polyfit(np.append(majiq_all, majiq_extra_all), np.append(rt_pcr_majiq_all, np.array(pcr_majiq_extra_all)), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y


    pearson_majiq = pearsonr(np.append(majiq_all, majiq_extra_all), np.append(rt_pcr_majiq_all, pcr_majiq_extra_all))[0]
    print 'MAJIQ R=%.2f' % pearson_majiq

    axx[0][0].text(.1, .9, 'R=%.2f' % (pearson_majiq), fontsize=14)
    axx[0][0].plot(np.append(majiq_all, majiq_extra_all), fit_fn(np.append(majiq_all, majiq_extra_all)), '--k')
    axx[0][0].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][0].plot(majiq[0], rt_pcr_majiq[0], '.', color=get_color('rest'), label='Resting')
    axx[0][0].plot(majiq[1], rt_pcr_majiq[1], '.', color=get_color('stim'), label='Stimulating')
    axx[0][0].plot(majiq_extra[0], pcr_majiq_extra[0], 'd', color=get_color('cer'), label='Cerebellum')
    axx[0][0].plot(majiq_extra[1], pcr_majiq_extra[1], 'd', color=get_color('liv'), label='Liver')


    axx[0][0].set_xlabel('MAJIQ')
    axx[0][0].set_ylabel('RT-PCR')
    axx[0][0].set_title('All (N=%d)' % (len(majiq_all) + len(majiq_extra_all)))
    axx[0][0].set_ylim([0,1])
    # axx[0][0].legend(loc=2, fontsize=8)

    diagonal = np.linspace(0, 1, num=len(rt_pcr_miso_all))
    fit = np.polyfit(np.append(miso_all, miso_extra_all), np.append(rt_pcr_miso_all, np.array(pcr_miso_extra_all)), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    pearson_miso = pearsonr(np.append(miso_all, miso_extra_all), np.append(rt_pcr_miso_all, pcr_miso_extra_all))[0]
    print 'MISO R=%.2f' % pearson_miso

    axx[1][0].text(.1, .9, 'R=%.2f' % (pearson_miso), fontsize=14)
    axx[1][0].plot(np.append(miso_all, miso_extra_all), fit_fn(np.append(miso_all, miso_extra_all)), '--k')
    axx[1][0].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][0].plot(miso[0], rt_pcr_miso[0], '.', color=get_color('rest'), label='Resting')
    axx[1][0].plot(miso[1], rt_pcr_miso[1], '.', color=get_color('stim'), label='Stimulating')
    axx[1][0].plot(miso_extra[0], pcr_miso_extra[0], 'd', color=get_color('cer'), label='Cerebellum')
    axx[1][0].plot(miso_extra[1], pcr_miso_extra[1], 'd', color=get_color('liv'), label='Liver')


    axx[1][0].set_xlabel('MISO')
    axx[1][0].set_ylabel('RT-PCR')
    axx[1][0].set_title('All (N=%d)' % (len(miso_all) + len(miso_extra_all)))
    axx[1][0].set_ylim([0,1])
    # axx[1][0].legend(loc=2, fontsize=8)


    diagonal = np.linspace(0, 1, num=len(rt_pcr_majiq[0]))
    fit = np.polyfit(majiq[0], rt_pcr_majiq[0], 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    pearson_r = pearsonr(majiq[0], rt_pcr_majiq[0])[0]
    axx[0][1].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[0][1].plot(majiq[0], fit_fn(majiq[0]), '--k')
    axx[0][1].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][1].plot(majiq[0], rt_pcr_majiq[0], '.', color=get_color('rest'), label='Resting')

    axx[0][1].set_xlabel('MAJIQ')
    axx[0][1].set_ylabel('RT-PCR')
    axx[0][1].set_title('Unstim (N=%d)' % len(majiq[0]))
    axx[0][1].set_ylim([0,1])
    # axx[0][0].legend(loc=4)

    fit = np.polyfit(majiq[1], rt_pcr_majiq[1], 1)
    fit_fn = np.poly1d(fit)
    pearson_r = pearsonr(majiq[1], rt_pcr_majiq[1])[0]
    axx[0][2].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[0][2].plot(majiq[1], fit_fn(majiq[1]), '--k')
    axx[0][2].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][2].plot(majiq[1], rt_pcr_majiq[1], '.', color=get_color('stim'), label='Stimulating')


    axx[0][2].set_xlabel('MAJIQ')
    axx[0][2].set_ylabel('RT-PCR')
    axx[0][2].set_title('Stim (N=%d)' % len(majiq[1]))
    axx[0][2].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(majiq_extra[0], pcr_majiq_extra[0], 1)
    fit_fn = np.poly1d(fit)
    pearson_r = pearsonr(majiq_extra[0], pcr_majiq_extra[0])[0]
    axx[0][3].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[0][3].plot(majiq_extra[0], fit_fn(majiq_extra[0]), '--k')
    axx[0][3].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][3].plot(majiq_extra[0], pcr_majiq_extra[0], 'd', color=get_color('cer'), label='Cerebellum')


    axx[0][3].set_xlabel('MAJIQ')
    axx[0][3].set_ylabel('RT-PCR')
    axx[0][3].set_title('Cerebellum (N=%d)' % len(majiq_extra[0]))
    axx[0][3].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(majiq_extra[1], pcr_majiq_extra[1], 1)
    fit_fn = np.poly1d(fit)
    pearson_r = pearsonr(majiq_extra[1], pcr_majiq_extra[1])[0]
    axx[0][4].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[0][4].plot(majiq_extra[1], fit_fn(majiq_extra[1]), '--k')
    axx[0][4].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][4].plot(majiq_extra[1], pcr_majiq_extra[1], 'd', color=get_color('liv'), label='Liver')


    axx[0][4].set_xlabel('MAJIQ')
    axx[0][4].set_ylabel('RT-PCR')
    axx[0][4].set_title('Liver (N=%d)' % len(majiq_extra[1]))
    axx[0][4].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(miso[0], rt_pcr_miso[0], 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    pearson_r = pearsonr(miso[0], rt_pcr_miso[0])[0]
    axx[1][1].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[1][1].plot(miso[0], fit_fn(miso[0]), '--k')
    axx[1][1].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][1].plot(miso[0], rt_pcr_miso[0], '.', color=get_color('rest'), label='Resting')

    axx[1][1].set_xlabel('MISO')
    axx[1][1].set_ylabel('RT-PCR')
    axx[1][1].set_title('Unstim (N=%d)' % len(miso[0]))
    axx[1][1].set_ylim([0,1])
    # axx[1][0].legend(loc=4)

    fit = np.polyfit(miso[1], rt_pcr_miso[1], 1)
    fit_fn = np.poly1d(fit)
    pearson_r = pearsonr(miso[1], rt_pcr_miso[1])[0]
    axx[1][2].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[1][2].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][2].plot(miso[1], rt_pcr_miso[1], '.', color=get_color('stim'), label='Stimulating')
    axx[1][2].plot(miso[1], fit_fn(miso[1]), '--k')

    axx[1][2].set_xlabel('MISO')
    axx[1][2].set_ylabel('RT-PCR')
    axx[1][2].set_title('Stim (N=%d)' % len(miso[1]))
    axx[1][2].set_ylim([0,1])
    # axx[1][1].legend(loc=4)

    fit = np.polyfit(miso_extra[0], pcr_miso_extra[0], 1)
    fit_fn = np.poly1d(fit)
    pearson_r = pearsonr(miso_extra[0], pcr_miso_extra[0])[0]
    axx[1][3].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[1][3].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][3].plot(miso_extra[0], pcr_miso_extra[0], 'd', color=get_color('cer'), label='Cerebellum')
    axx[1][3].plot(miso_extra[0], fit_fn(miso_extra[0]), '--k')

    axx[1][3].set_xlabel('MISO')
    axx[1][3].set_ylabel('RT-PCR')
    axx[1][3].set_title('Cerebellum (N=%d)' % len(miso_extra[0]))
    axx[1][3].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(miso_extra[1], pcr_miso_extra[1], 1)
    fit_fn = np.poly1d(fit)
    pearson_r = pearsonr(miso_extra[1], pcr_miso_extra[1])[0]
    axx[1][4].text(.1, .9, 'R=%.2f' % (pearson_r), fontsize=14)
    axx[1][4].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][4].plot(miso_extra[1], pcr_miso_extra[1], 'd', color=get_color('liv'), label='Liver')
    axx[1][4].plot(miso_extra[1], fit_fn(miso_extra[1]), '--k')

    axx[1][4].set_xlabel('MISO')
    axx[1][4].set_ylabel('RT-PCR')
    axx[1][4].set_title('Liver (N=%d)' % len(miso_extra[1]))
    axx[1][4].set_ylim([0,1])

    scripts.utils.save_or_show(plotpath, "psi_comp_rtpcr_majiq_miso", exten='pdf')


def scatterplot_rtpcr_simple(rt_pcr, method_epsis, plotpath, pcr_majiq_extra=None, majiq_extra=None, plotname='psi_rtpcr', met_name='MAJIQ'):
    #figure out how many groups of events exist

    # majiq_rest_yerr = [var_expected_psi(dis) for dis in rt_pcr_majiq[0]]

    fig = pyplot.figure(figsize=[6, 6], dpi=300)
    fig.suptitle("RT-PCR PSI comparison")

    rt_pcr_set1_all = [a for b in rt_pcr for a in b]
    method_set1_all = [a for b in method_epsis for a in b]

    rt_pcr_set2_all = [a for b in pcr_majiq_extra for a in b]
    method_set2_all = [a for b in majiq_extra for a in b]


    diagonal = np.linspace(0, 1, num=len(rt_pcr_set1_all))
    fit = np.polyfit(np.append(method_set1_all, method_set2_all), np.append(rt_pcr_set1_all, np.array(rt_pcr_set2_all)), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    pyplot.plot(np.append(method_set1_all, method_set2_all), fit_fn(np.append(method_set1_all, method_set2_all)), '--k')
    pyplot.plot(diagonal, diagonal, '--', color="#cccccc")
    pyplot.plot(method_epsis[0], rt_pcr[0], '.', color=get_color('rest'), label='Resting')
    pyplot.plot(method_epsis[1], rt_pcr[1], '.', color=get_color('stim'), label='Stimulating')
    pyplot.plot(majiq_extra[0], pcr_majiq_extra[0], 'd', color=get_color('cer'), label='Cerebellum')
    pyplot.plot(majiq_extra[1], pcr_majiq_extra[1], 'd', color=get_color('liv'), label='Liver')

    pyplot.xlabel(met_name)
    pyplot.ylabel('RT-PCR')
    pyplot.title('All (N=%d)' % (len(method_set1_all) + len(method_set2_all)))
    pyplot.ylim([0, 1])
    scripts.utils.save_or_show(plotpath, plotname, exten='pdf')


def barchart_expected(expected_psis, plotpath, mfile):
    pyplot.figure()
    pyplot.hist(expected_psis, bins=40)
    name = os.path.basename(mfile)
    scripts.utils.save_or_show(plotpath, name +'_expected_dist')


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
                pcr_rest_stim[pcr_fields[0]].append([float(pcr_fields[1])/100, float(pcr_fields[2])/100])
            except IndexError:
                print "%s value not found, assigned 0..." % pcr_fields[0]
                pcr_rest_stim[pcr_fields[0]].append([0, 0])
            except ValueError, e:
                print e.message
                print "Event PSI could not be converted into float, wrong parsing??"
                pcr_rest_stim[pcr_fields[0]].append([0, 0])
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
                        majiq_dict[lsv.get_id()+"#"+str(jcoords[2])].append(lsv.get_bins()[j])
                    except IndexError:
                        majiq_dict[lsv.get_id()+"#"+str(jcoords[2])].append(lsv.get_bins()[0][::-1])

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


def main():
    parser = argparse.ArgumentParser(description="Compare PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--majiq-rest", required=True, dest='majiq_rest', nargs='+', help='MAJIQ PSI predictions for resting RNA-Seq data.')
    parser.add_argument("--majiq-stim", required=True, dest='majiq_stim', nargs='+', help='MAJIQ PSI predictions for stimulated RNA-Seq data.')
    parser.add_argument("--miso-rest", dest='miso_rest', nargs='*', help='MISO PSI predictions for resting RNA-Seq data.')
    parser.add_argument("--miso-stim", dest='miso_stim', nargs='*', help='MISO PSI predictions for stimulated RNA-Seq data.')
    parser.add_argument("--names-map-file", required=True, dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')

    parser.add_argument("--majiq-extra", dest='majiq_extra', nargs='+', help='MAJIQ PSI predictions for extra conditions.')
    parser.add_argument("--miso-extra", dest='miso_extra', nargs='+', help='MISO PSI predictions for extra conditions.')
    parser.add_argument("--pcr-extra", dest='pcr_extra', help='RT-PCR validations for extra conditions.')


    # parser.add_argument("--builder-file", required=True, dest='majiq_builder_file', help='File containing MAJIQ builder output file.')
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    # 6. Generate vectors/lists, saving our predictions in a txt file
    lrtpcr_majiq_extra=[[], []]
    lrtpcr_miso_extra=[[], []]
    lmajiq_extra = [[], []]
    lmiso_extra=[[], []]

    if args.majiq_extra:
        # 1. Read PCR results
        rtpcr_extra = defaultdict(list)
        with open(args.pcr_extra) as pcr_extra:
            for pcr_elem in pcr_extra:
                if pcr_elem.startswith("#"): continue
                pcr_elems=pcr_elem.rstrip().split()
                rtpcr_extra[pcr_elems[0].split("#")[0]]=[float(pcr_elems[4])/100, float(pcr_elems[3])/100, np.nan, np.nan, [], [], float(pcr_elems[6])/100, float(pcr_elems[5])/100, pcr_elems[-4]]


        # Read Majiq results for the elements in PCR
        djunc_selected = {}
        for cn, mfile in enumerate(args.majiq_extra):
            majiq_found = []
            with open(mfile) as mfile_open:
                mpickle = pickle.load(mfile_open)
                for j, lsv in enumerate(mpickle.get_lsvs()):
                    # if lsv.get_id().split(":")[0] in ensembl_tlb.keys():
                    #     for trans in set(ensembl_tlb[lsv.get_id().split(":")[0]]):
                    if lsv.get_id() in rtpcr_extra.keys():
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
                                    rtpcr_extra[lsv.get_id()][2+cn] = expected_psi(lsv.get_bins()[lway_aux])
                                except IndexError:
                                    rtpcr_extra[lsv.get_id()][2+cn] = 1-expected_psi(lsv.get_bins()[0])
            print "Missing LSVs in %s:" % mfile
            print '\n'.join([kk for kk in rtpcr_extra.keys() if kk not in majiq_found])

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
                        if lsv_name in rtpcr_extra.keys():
                            # Find all junctions from that LSV that are included
                            try:
                                rtpcr_extra[lsv_name][4+cn].append(float(miso_psis[djunc_selected[lsv_name]]))
                            except IndexError:
                                rtpcr_extra[lsv_name][4+cn].append(1-float(miso_psis[0]))
                            except KeyError:
                                print "[WARNING] :: %s in MISO, but not in MAJIQ, skipped given the impossibility of determining which junction it is." % lsv_name

        with open('psi_hogenesch.txt', 'w') as psi_txt:
            headers=['#', 'ID (transcript)', 'ID (gene)', 'RT-PCR Cerebellum', 'RT-PCR Liver', 'Majiq Cerebellum', 'Majiq Liver', 'Miso Cerebellum', 'Miso Liver', 'RT-PCR Cerebellum Avg. STD', 'RT-PCR Liver Avg. STD']
            psi_txt.write('\t'.join(headers))
            psi_txt.write('\n')
            for ilsv, vals in rtpcr_extra.iteritems():
                line = [ilsv]
                line.append(ilsv)
                line.extend([repr(vv) for vv in vals])
                psi_txt.write("\t".join(line))
                psi_txt.write('\n')
                if not np.isnan(vals[2]):
                    lrtpcr_majiq_extra[0].append(vals[0])
                    # Majiq
                    lmajiq_extra[0].append(vals[2])
                if not np.isnan(vals[3]):
                    lrtpcr_majiq_extra[1].append(vals[1])
                    # Majiq
                    lmajiq_extra[1].append(vals[3])
                if not np.isnan(np.mean(vals[4])):
                    lrtpcr_miso_extra[0].append(vals[0])
                    # Majiq
                    lmiso_extra[0].append(np.mean(vals[4]))
                if not np.isnan(np.mean(vals[5])):
                    lrtpcr_miso_extra[1].append(vals[1])
                    # Majiq
                    lmiso_extra[1].append(np.mean(vals[5]))

    # Read name mapping file
    names_pcr2majiq_dict = {}  # key RT-PCR, value MAJIQ
    names_majiq2pcr_dict = {}
    names_junc_majiq = defaultdict(list)
    gene_names_counts = defaultdict(lambda: 0)
    with open(args.names_map_file) as names_map:
        for name_map in names_map:
            # MAJIQ
            mapped_name = name_map.rstrip().split()
            majiq_name, junc_idx_old = mapped_name[1].split('#')
            names_junc_majiq[majiq_name].append([int(jcoord) for jcoord in  mapped_name[-1].split('-')] + [int(junc_idx_old)])

            # RT-PCR
            names_pcr2majiq_dict[mapped_name[0]] = mapped_name[1]
            names_majiq2pcr_dict[mapped_name[1]] = mapped_name[0]
            gene_names_counts[mapped_name[0]] += 1

    print "Number of events in names-map-file RT-PCR: %d" % len(names_pcr2majiq_dict.keys())
    print "Number of events in names-map-file MAJIQ: %d" % len(names_junc_majiq.keys())

    # Parse RT-PCR results
    pcr_rest_stim = parse_rtpcr_results(args.pcr, names_pcr2majiq_dict)

    # Process MAJIQ files for resting RNA-Seq data
    majiq_rest_dict = parse_majiq_results(args.majiq_rest, names_junc_majiq)

    # Process MAJIQ files for stimuli RNA-Seq data
    majiq_stim_dict = parse_majiq_results(args.majiq_stim, names_junc_majiq)

    # Process MAJIQ files for resting RNA-Seq data
    miso_rest_dict = parse_miso_results(args.miso_rest, names_junc_majiq)

    # Process MAJIQ files for stimuli RNA-Seq data
    miso_stim_dict = parse_miso_results(args.miso_stim, names_junc_majiq)

    ## Intersect names from RT-PCR and MAJIQ
    common_names_set = set([names_pcr2majiq_dict[k] for k in pcr_rest_stim.keys()]).intersection(set(majiq_rest_dict.keys())).intersection(set(majiq_stim_dict.keys()))
    print "Common names afer intersection with MAJIQ: %d" % len(common_names_set)
    common_names_set = common_names_set.intersection(set(miso_rest_dict.keys())).intersection(set(miso_stim_dict.keys()))
    print "Common names afer intersection with MISO: %d" % len(common_names_set)

    rt_pcr_majiq = [[], []]
    rt_pcr_miso = [[], []]
    majiq = [[], []]
    miso = [[], []]

    flipped_thres = 1
    flipped_majiq_dict = defaultdict(str)  # List of strings with flipped LSVs info
    flipped_miso_dict = defaultdict(str)  # List of strings with flipped LSVs info
    for common_name in common_names_set:
        for name_majiq, name_pcr in names_majiq2pcr_dict.iteritems():
            if names_majiq2pcr_dict[common_name] == name_pcr:
                name = name_majiq

                # For Majiq, compute mean over Expected PSIs
                majiq_rest_stat = avg_expected_psi(majiq_rest_dict[name])
                majiq_stim_stat = avg_expected_psi(majiq_stim_dict[name])

                # For MISO, compute mean
                miso_rest_stat = np.mean(miso_rest_dict[name])
                miso_stim_stat = np.mean(miso_stim_dict[name])

                
                print "%s - %s" % (name, names_majiq2pcr_dict[name])
                print "---- RT-PCR ----"

                rtpcr_rest = pcr_rest_stim[names_majiq2pcr_dict[name]][0][0]
                rtpcr_stim = pcr_rest_stim[names_majiq2pcr_dict[name]][0][1]
                min_rest = abs(rtpcr_rest - majiq_rest_stat)
                min_stim = abs(rtpcr_stim - majiq_stim_stat)
                for rtpcr_psi_value in pcr_rest_stim[names_majiq2pcr_dict[name]]:
                    if abs(rtpcr_psi_value[0] - majiq_rest_stat) < min_rest:
                        rtpcr_rest = rtpcr_psi_value[0]
                        min_rest = abs(rtpcr_rest - majiq_rest_stat)

                    if abs(rtpcr_psi_value[1] - majiq_rest_stat) < min_stim:
                        rtpcr_stim = rtpcr_psi_value[1]
                        min_stim = abs(rtpcr_stim - majiq_stim_stat)

                rt_pcr_majiq[0].append(rtpcr_rest)
                rt_pcr_majiq[1].append(rtpcr_stim)

                rtpcr_rest = pcr_rest_stim[names_majiq2pcr_dict[name]][0][0]
                rtpcr_stim = pcr_rest_stim[names_majiq2pcr_dict[name]][0][1]
                min_rest = abs(rtpcr_rest - miso_rest_stat)
                min_stim = abs(rtpcr_stim - miso_stim_stat)
                for rtpcr_psi_value in pcr_rest_stim[names_majiq2pcr_dict[name]]:
                    if abs(rtpcr_psi_value[0] - miso_rest_stat) < min_rest:
                        rtpcr_rest = rtpcr_psi_value[0]
                        min_rest = abs(rtpcr_rest - miso_rest_stat)

                    if abs(rtpcr_psi_value[1] - miso_rest_stat) < min_stim:
                        rtpcr_stim = rtpcr_psi_value[1]
                        min_stim = abs(rtpcr_stim - miso_stim_stat)

                rt_pcr_miso[0].append(rtpcr_rest)
                rt_pcr_miso[1].append(rtpcr_stim)

                # check if event has expected PSIs suspicious of being flipped
                if abs(majiq_rest_stat - rt_pcr_majiq[0][-1]) > flipped_thres or abs(majiq_stim_stat - rt_pcr_majiq[1][-1]) > flipped_thres:
                    flipped_majiq_dict[name] = "%s\t%s\t%f\t%f\t%f\t%f\t%d\t%d\n" % (names_majiq2pcr_dict[name], name, rt_pcr_majiq[0][-1], majiq_rest_stat, rt_pcr_majiq[1][-1], majiq_stim_stat, int(gene_names_counts[names_majiq2pcr_dict[name]]<2), int(len(names_junc_majiq[str(name).split('#')[0]])<2) )
                    del rt_pcr_majiq[0][-1]
                    del rt_pcr_majiq[1][-1]
                else:
                    majiq[0].append(majiq_rest_stat)
                    majiq[1].append(majiq_stim_stat)

                if abs(miso_rest_stat - rt_pcr_miso[0][-1]) > flipped_thres or abs(miso_stim_stat - rt_pcr_miso[1][-1]) > flipped_thres:
                    flipped_miso_dict[name] = "%s\t%s\t%f\t%f\t%f\t%f\t%d\t%d\n" % (names_majiq2pcr_dict[name], name, rt_pcr_miso[0][-1], miso_rest_stat, rt_pcr_miso[1][-1], miso_stim_stat, int(gene_names_counts[names_majiq2pcr_dict[name]]<2), int(len(names_junc_majiq[str(name).split('#')[0]])<2) )
                    del rt_pcr_miso[0][-1]
                    del rt_pcr_miso[1][-1]
                elif np.isnan(miso_rest_stat) or np.isnan(miso_stim_stat):
                    del rt_pcr_miso[0][-1]
                    del rt_pcr_miso[1][-1]
                else:
                    miso[0].append(miso_rest_stat)
                    miso[1].append(miso_stim_stat)


                print "[Resting - MAJIQ]:\t%f" % rtpcr_rest
                print "[Stimuli - MAJIQ]:\t%f" % rtpcr_stim

                print "---- MAJIQ ----"
                print "[Resting]: Mean of expected:\t%f" % (float(majiq_rest_stat ))
                print "[Stimuli]: Mean of expected:\t%f" % (float(majiq_stim_stat))

                print "---- MISO -----"
                print "[Resting]: Mean of psi:\t%f" % (float(miso_rest_stat))
                print "[Stimuli]: Mean of psi:\t%f" % (float(miso_stim_stat))

    scatterplot_rtpcr_majiq_miso(rt_pcr_majiq, rt_pcr_miso, majiq, miso, args.plotpath, pcr_majiq_extra=lrtpcr_majiq_extra, pcr_miso_extra=lrtpcr_miso_extra, majiq_extra=lmajiq_extra, miso_extra=lmiso_extra)
    scatterplot_rtpcr_simple(rt_pcr_majiq, majiq, args.plotpath, pcr_majiq_extra=lrtpcr_majiq_extra,majiq_extra=lmajiq_extra, plotname='psi_majiq_only', met_name='MAJIQ')
    scatterplot_rtpcr_simple(rt_pcr_miso, miso, args.plotpath, pcr_majiq_extra=lrtpcr_miso_extra,majiq_extra=lmiso_extra, plotname='psi_miso_only', met_name='MISO')

if __name__ == '__main__':
    main()