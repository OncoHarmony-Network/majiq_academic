from collections import defaultdict
import cPickle as pickle
import numpy as np
import argparse
import os
import matplotlib.pyplot as pyplot
import scripts.utils
import colorbrewer as cb

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


def scatterplot_rtpcr_majiq(rt_pcr_majiq, rt_pcr_miso, majiq, miso, plotpath, pcr_majiq_extra=None, pcr_miso_extra=None, majiq_extra=None, miso_extra=None):
    #figure out how many groups of events exist

    # majiq_rest_yerr = [var_expected_psi(dis) for dis in rt_pcr_majiq[0]]
    import matplotlib
    matplotlib.rcParams.update({'font.size': 10})

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
    # print majiq[0], rt_pcr_majiq[0], miso[0]
    # fit = np.polyfit(majiq_all, rt_pcr_majiq_all, 1)
    fit = np.polyfit(np.append(majiq_all, majiq_extra_all), np.append(rt_pcr_majiq_all, np.array(pcr_majiq_extra_all)), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[0][0].plot(np.append(majiq_all, majiq_extra_all), fit_fn(np.append(majiq_all, majiq_extra_all)), '--k')
    axx[0][0].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][0].plot(majiq[0], rt_pcr_majiq[0], '.', color='r', label='Resting')
    axx[0][0].plot(majiq[1], rt_pcr_majiq[1], '.', color='b', label='Stimulating')
    axx[0][0].plot(majiq_extra[0], pcr_majiq_extra[0], '*', color='#%02x%02x%02x' % cb.Dark2[4][2], label='Cerebellum')
    axx[0][0].plot(majiq_extra[1], pcr_majiq_extra[1], '*', color='#%02x%02x%02x' % cb.Dark2[4][3], label='Liver')


    axx[0][0].set_xlabel('MAJIQ')
    axx[0][0].set_ylabel('RT-PCR')
    axx[0][0].set_title('All (N=%d)' % (len(majiq_all) + len(majiq_extra_all)))
    axx[0][0].set_ylim([0,1])
    # axx[0][0].legend(loc=2, fontsize=8)

    diagonal = np.linspace(0, 1, num=len(rt_pcr_miso_all))
    # fit = np.polyfit(miso_all, rt_pcr_miso_all, 1)
    fit = np.polyfit(np.append(miso_all, miso_extra_all), np.append(rt_pcr_miso_all, np.array(pcr_miso_extra_all)), 1)
    # fit = np.polyfit(np.append(miso_all, np.array(miso_extra)[-np.isnan(miso_extra)]), np.append(rt_pcr_miso_all, np.array(pcr_extra)[-np.isnan(miso_extra)]), 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[1][0].plot(np.append(miso_all, miso_extra_all), fit_fn(np.append(miso_all, miso_extra_all)), '--k')
    axx[1][0].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][0].plot(miso[0], rt_pcr_miso[0], '.', color='r', label='Resting')
    axx[1][0].plot(miso[1], rt_pcr_miso[1], '.', color='b', label='Stimulating')
    axx[1][0].plot(miso_extra[0], pcr_miso_extra[0], '*', color='#%02x%02x%02x' % cb.Dark2[4][2], label='Cerebellum')
    axx[1][0].plot(miso_extra[1], pcr_miso_extra[1], '*', color='#%02x%02x%02x' % cb.Dark2[4][3], label='Liver')


    axx[1][0].set_xlabel('MISO')
    axx[1][0].set_ylabel('RT-PCR')
    axx[1][0].set_title('All (N=%d)' % (len(miso_all) + len(miso_extra_all)))
    axx[1][0].set_ylim([0,1])
    # axx[1][0].legend(loc=2, fontsize=8)


    diagonal = np.linspace(0, 1, num=len(rt_pcr_majiq[0]))
    fit = np.polyfit(majiq[0], rt_pcr_majiq[0], 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[0][1].plot(majiq[0], fit_fn(majiq[0]), '--k')
    axx[0][1].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][1].plot(majiq[0], rt_pcr_majiq[0], '.', color='r', label='Resting')

    axx[0][1].set_xlabel('MAJIQ')
    axx[0][1].set_ylabel('RT-PCR')
    axx[0][1].set_title('Resting (N=%d)' % len(majiq[0]))
    axx[0][1].set_ylim([0,1])
    # axx[0][0].legend(loc=4)

    fit = np.polyfit(majiq[1], rt_pcr_majiq[1], 1)
    fit_fn = np.poly1d(fit)
    axx[0][2].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][2].plot(majiq[1], rt_pcr_majiq[1], '.', color='b', label='Stimulating')
    axx[0][2].plot(majiq[1], fit_fn(majiq[1]), '--k')

    axx[0][2].set_xlabel('MAJIQ')
    axx[0][2].set_ylabel('RT-PCR')
    axx[0][2].set_title('Stimulating (N=%d)' % len(majiq[1]))
    axx[0][2].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(majiq_extra[0], pcr_majiq_extra[0], 1)
    fit_fn = np.poly1d(fit)
    axx[0][3].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][3].plot(majiq_extra[0], pcr_majiq_extra[0], '*', color='#%02x%02x%02x' % cb.Dark2[4][2], label='Cerebellum')
    axx[0][3].plot(majiq_extra[0], fit_fn(majiq_extra[0]), '--k')

    axx[0][3].set_xlabel('MAJIQ')
    axx[0][3].set_ylabel('RT-PCR')
    axx[0][3].set_title('Cerebellum (N=%d)' % len(majiq_extra[0]))
    axx[0][3].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(majiq_extra[1], pcr_majiq_extra[1], 1)
    fit_fn = np.poly1d(fit)
    axx[0][4].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[0][4].plot(majiq_extra[1], pcr_majiq_extra[1], '*', color='#%02x%02x%02x' % cb.Dark2[4][3], label='Liver')
    axx[0][4].plot(majiq_extra[1], fit_fn(majiq_extra[1]), '--k')

    axx[0][4].set_xlabel('MAJIQ')
    axx[0][4].set_ylabel('RT-PCR')
    axx[0][4].set_title('Liver (N=%d)' % len(majiq_extra[1]))
    axx[0][4].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(miso[0], rt_pcr_miso[0], 1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y

    axx[1][1].plot(miso[0], fit_fn(miso[0]), '--k')
    axx[1][1].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][1].plot(miso[0], rt_pcr_miso[0], '.', color='r', label='Resting')

    axx[1][1].set_xlabel('MISO')
    axx[1][1].set_ylabel('RT-PCR')
    axx[1][1].set_title('Resting (N=%d)' % len(miso[0]))
    axx[1][1].set_ylim([0,1])
    # axx[1][0].legend(loc=4)

    fit = np.polyfit(miso[1], rt_pcr_miso[1], 1)
    fit_fn = np.poly1d(fit)
    axx[1][2].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][2].plot(miso[1], rt_pcr_miso[1], '.', color='b', label='Stimulating')
    axx[1][2].plot(miso[1], fit_fn(miso[1]), '--k')

    axx[1][2].set_xlabel('MISO')
    axx[1][2].set_ylabel('RT-PCR')
    axx[1][2].set_title('Stimulating (N=%d)' % len(miso[1]))
    axx[1][2].set_ylim([0,1])
    # axx[1][1].legend(loc=4)

    fit = np.polyfit(miso_extra[0], pcr_miso_extra[0], 1)
    fit_fn = np.poly1d(fit)
    axx[1][3].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][3].plot(miso_extra[0], pcr_miso_extra[0], '*', color='#%02x%02x%02x' % cb.Dark2[4][2], label='Cerebellum')
    axx[1][3].plot(miso_extra[0], fit_fn(miso_extra[0]), '--k')

    axx[1][3].set_xlabel('MISO')
    axx[1][3].set_ylabel('RT-PCR')
    axx[1][3].set_title('Cerebellum (N=%d)' % len(miso_extra[0]))
    axx[1][3].set_ylim([0,1])
    # axx[0][1].legend(loc=4)

    fit = np.polyfit(miso_extra[1], pcr_miso_extra[1], 1)
    fit_fn = np.poly1d(fit)
    axx[1][4].plot(diagonal, diagonal, '--', color="#cccccc")
    axx[1][4].plot(miso_extra[1], pcr_miso_extra[1], '*', color='#%02x%02x%02x' % cb.Dark2[4][3], label='Liver')
    axx[1][4].plot(miso_extra[1], fit_fn(miso_extra[1]), '--k')

    axx[1][4].set_xlabel('MISO')
    axx[1][4].set_ylabel('RT-PCR')
    axx[1][4].set_title('Liver (N=%d)' % len(miso_extra[1]))
    axx[1][4].set_ylim([0,1])

    _save_or_show(plotpath, "psi_comp_rtpcr_majiq_miso")


def barchart_expected(expected_psis, plotpath, mfile):
    pyplot.figure()
    pyplot.hist(expected_psis, bins=40)
    name = os.path.basename(mfile)
    _save_or_show(plotpath, name +'_expected_dist')


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
            if i<1: continue  # headers
            pcr_fields = pcr_line.rstrip().split()

            if pcr_fields[0] not in names_pcr2majiq_dict.keys(): continue  # If the event is not in the list of events selected, skip it
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
    files_majiq = scripts.utils.list_files_or_dir(files_majiq, containing='majiq')
    for mfile in files_majiq:
        expected_psis = []
        with open(mfile) as mfile_open:
            mpickle = pickle.load(mfile_open)
            for i, lsv_info in enumerate(mpickle[1]):
                if lsv_info[1] not in names_junc_majiq.keys():
                    continue # If the event is not in the list of events selected, skip it
                if len(mpickle[0][i])<3:
                    expected_psis.append(expected_psi(mpickle[0][i][0]))

                # Find all junctions from that LSV that are included
                for j, lsv_way in enumerate(mpickle[0][i]):
                    if j in names_junc_majiq[lsv_info[1]]:
                        majiq_dict[lsv_info[1]+"#"+str(j)].append(mpickle[0][i][j])
        # barchart_expected(expected_psis, args.plotpath, mfile)
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


def main():
    parser = argparse.ArgumentParser(description="Compare PSIs computed with MAJIQ against the RT-PCR results.")
    parser.add_argument("pcr", help="Tab-delimted file with the RT-PCR scores")
    parser.add_argument("--majiq-rest", required=True, dest='majiq_rest', nargs='+', help='MAJIQ PSI predictions for resting RNA-Seq data.')
    parser.add_argument("--majiq-stim", required=True, dest='majiq_stim', nargs='+', help='MAJIQ PSI predictions for stimulated RNA-Seq data.')
    parser.add_argument("--miso-rest", required=True, dest='miso_rest', nargs='*', help='MISO PSI predictions for resting RNA-Seq data.')
    parser.add_argument("--miso-stim", required=True, dest='miso_stim', nargs='*', help='MISO PSI predictions for stimulated RNA-Seq data.')
    parser.add_argument("--names-map-file", required=True, dest='names_map_file', help='File containing the mapping for events names used in MAJIQ and RT-PCR files.')

    parser.add_argument("--majiq-extra", dest='majiq_extra', nargs='+', help='MAJIQ PSI predictions for extra conditions.')
    parser.add_argument("--miso-extra", dest='miso_extra', nargs='+', help='MISO PSI predictions for extra conditions.')
    parser.add_argument("--pcr-extra", dest='pcr_extra', help='RT-PCR validations for extra conditions.')
    parser.add_argument("--tlb-extra", dest='tlb_extra', help='TLB for ENSEMBL gene and transcript names.')
    parser.add_argument("--tlb-extra-pkl", dest='tlb_extra_pkl', help='TLB for ENSEMBL gene and transcript names.')
    parser.add_argument("--juncs-extra", dest='junc_extra', help='Junction coordinates txt file.')
    parser.add_argument("--voila-extra", dest='voila_extra', nargs=2, help='Junction coordinates txt file.')


    # parser.add_argument("--builder-file", required=True, dest='majiq_builder_file', help='File containing MAJIQ builder output file.')
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    if args.majiq_extra:
        # 1. Read PCR results
        rtpcr_extra = defaultdict(list)
        with open(args.pcr_extra) as pcr_extra:
            for pcr_elem in pcr_extra:
                if pcr_elem.startswith("#"): continue
                pcr_elems=pcr_elem.rstrip().split()
                rtpcr_extra[pcr_elems[0]]=[float(pcr_elems[2])/100, float(pcr_elems[1])/100, np.nan, np.nan, [], []]

        ensembl_tlb = defaultdict(list)
        if args.tlb_extra_pkl:
            ensembl_tlb = pickle.load(open(args.tlb_extra_pkl))
        else:
            with open(args.tlb_extra) as tlb_f:
                rtpcr_keys = [rk.split(":")[0] for rk in rtpcr_extra.keys()]
                for tlb_line in tlb_f:
                    fields=tlb_line.split()
                    if len(fields)<2: continue
                    transcript=fields[1].strip(";").strip("\"")
                    gene=fields[0].strip(";").strip("\"")
                    if transcript in rtpcr_keys:
                        ensembl_tlb[gene].append(transcript)
            pickle.dump(ensembl_tlb, open("ensembl.tlb.pickle", "w"))
        # for k, v in ensembl_tlb.iteritems():
        #     print "%s: %s" % (k,repr(v))

        # 2. Read junctions coordinates file
        djunc_extra = defaultdict(list)
        with open(args.junc_extra) as junc_extra:
            for junc in junc_extra:
                if junc.startswith("#"): continue
                juncs=junc.rstrip().split()
                djunc_extra[juncs[0]]=juncs[1:]

        # 3. Read voila Tab-delimited output file
        djuncs = {}
        for fvoila in args.voila_extra:
            with open(fvoila) as voilaf:
                for voilal in voilaf:
                    if voilal.startswith("#"): continue
                    vfields=voilal.rstrip().split()
                    djuncs[vfields[0]] = vfields[1].split(";")  # Store LSV ID; junction coordinates list

        # 4. Read Majiq results for the elements in PCR
        for cn, mfile in enumerate(args.majiq_extra):
            with open(mfile) as mfile_open:
                mpickle = pickle.load(mfile_open)
                for i, lsv_info in enumerate(mpickle[1]):
                    for j, lsv_way in enumerate(mpickle[0][i]):
                        if lsv_info[1].split(":")[0] in ensembl_tlb.keys():
                            for trans in set(ensembl_tlb[lsv_info[1].split(":")[0]]):
                                for lway_aux in range(2):
                                    lsv_aux = trans+":"+":".join(lsv_info[1].split(":")[1:]) + "#%d"%lway_aux
                                    if lsv_aux in rtpcr_extra.keys():
                                        for ji, jcoords in enumerate(djuncs[lsv_info[1]]):
                                            if djunc_extra[lsv_aux.split("#")[0]][lway_aux] == jcoords:
                                                try:
                                                    rtpcr_extra[lsv_aux][2+cn] = expected_psi(mpickle[0][i][ji])
                                                    print "Complex LSV %s" % lsv_aux
                                                except IndexError:
                                                    rtpcr_extra[lsv_aux][2+cn] = 1-expected_psi(mpickle[0][i][0])
                                                    print "Binary LSV %s" % lsv_aux
                                                continue

        # 5. Read MISO results for the elements in PCR
        for cn, miso_dir in enumerate(args.miso_extra):
            files_miso = scripts.utils.list_files_or_dir([miso_dir], containing='miso')
            for mfile in files_miso:
                with open(mfile) as mfile_open:
                    for line in mfile_open:
                        miso_fields = line.rstrip().split()
                        lsv_name = miso_fields[0]
                        miso_psis = miso_fields[1]
                        # Find all junctions from that LSV that are included
                        for i, miso_psi in enumerate(miso_psis.split(',')):
                            if lsv_name+"#"+str(i) in rtpcr_extra.keys():
                                rtpcr_extra[lsv_name+"#"+str(i)][4+cn].append(float(miso_psi))
                            if lsv_name+"#"+str(1) in rtpcr_extra.keys():
                                if len(miso_psis.split(',')) > 1:
                                    rtpcr_extra[lsv_name+"#"+str(1)][4+cn].append(float(miso_psi))
                                else:
                                    rtpcr_extra[lsv_name+"#"+str(1)][4+cn].append(1-float(miso_psi))



        # 6. Generate vectors/lists, saving our predictions in a txt file
        lrtpcr_majiq_extra=[[], []]
        lrtpcr_miso_extra=[[], []]
        lmajiq_extra = [[], []]
        lmiso_extra=[[], []]

        with open('psi_hogenesch.txt', 'w') as psi_txt:
            headers=['#', 'ID (transcript)', 'ID (gene)', 'RT-PCR Cerebellum', 'RT-PCR Liver', 'Majiq Cerebellum', 'Majiq Liver', 'Miso Cerebellum', 'Miso Liver']
            psi_txt.write('\t'.join(headers))
            psi_txt.write('\n')
            for ilsv, vals in rtpcr_extra.iteritems():
                line = [ilsv]
                line.append("%s:%s" % ([ensg for ensg, enst in ensembl_tlb.items() if ilsv.split(":")[0] in enst][0], ":".join(ilsv.split(":")[1:])))
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
    gene_names_counts = defaultdict(lambda : 0)
    import ast
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
                majiq_rest_dist = majiq_rest_dict[name]
                majiq_stim_dist = majiq_stim_dict[name]

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

                # rt_pcr[0].append(pcr_rest_stim[names_majiq2pcr_dict[name]][0])
                # rt_pcr[1].append(pcr_rest_stim[names_majiq2pcr_dict[name]][1])
                rt_pcr_majiq[0].append(rtpcr_rest)
                rt_pcr_majiq[1].append(rtpcr_stim)
                # print "[Resting]:\t%f" % (pcr_rest_stim[names_majiq2pcr_dict[name]][0])
                # print "[Stimuli]:\t%f" % (pcr_rest_stim[names_majiq2pcr_dict[name]][1])

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
                    # majiq[0].append(majiq_rest_dist)
                    # majiq[1].append(majiq_stim_dist)

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
                # print "[Resting]: Expected delta psi statistic: %f" % (float(expected_delta_psi(majiq_rest_dict[name], pcr_rest_stim[names_majiq2pcr_dict[name]][0])))
                print "[Stimuli]: Mean of expected:\t%f" % (float(majiq_stim_stat))
                # print "[Stimuli]: Expected delta psi statistic: %f" % (float(expected_delta_psi(majiq_stim_dict[name], pcr_rest_stim[names_majiq2pcr_dict[name]][1])))

                print "---- MISO -----"
                print "[Resting]: Mean of psi:\t%f" % (float(miso_rest_stat))
                print "[Stimuli]: Mean of psi:\t%f" % (float(miso_stim_stat))


                # print "[Resting]: Mean of %s:\t%f" % (majiq_rest_dict[name], float(majiq_rest_stat))
                # print "[Stimuli]: Mean of %s:\t%f" % (majiq_stim_dict[name], float(majiq_stim_stat))


    # print repr(rt_pcr), repr(majiq)



    scatterplot_rtpcr_majiq(rt_pcr_majiq, rt_pcr_miso, majiq, miso, args.plotpath, pcr_majiq_extra=lrtpcr_majiq_extra, pcr_miso_extra=lrtpcr_miso_extra, majiq_extra=lmajiq_extra, miso_extra=lmiso_extra)

    # Save presumably flipped events in majiq
    # flipped_lsv_names = [str(name_str).split("#")[0] for name_str in flipped_majiq_dict.keys()]
    with open('flipped_majiq_thres_%.2f.txt' % flipped_thres, 'w') as flipped_file:
        flipped_file.write("name_rtpcr\tname_majiq\trtpcr_rest%f\tmajiq_rest\trtpcr_stim\tmajiq_stim\tunique_rtpcr_name?\tlsv_1_way?\ttarget_junction_coord\n")
        # with open(args.majiq_builder_file) as majiq_builder_file:
        #     majiq_builder = pickle.load(majiq_builder_file)
        #     for lsv in majiq_builder[1]:
        #         if lsv.id in flipped_lsv_names:
        #             for flip_key in flipped_majiq_dict:
        #                 if lsv.id in flip_key:
        #                     flipped_majiq_dict[flip_key] = flipped_majiq_dict[flip_key] + "\t%s\n" % str(lsv.junction_id[int(str(flip_key).split('#')[1])])
        for flipped_lsv in flipped_majiq_dict:
            flipped_file.write(flipped_majiq_dict[flipped_lsv])

    # Save presumably flipped events in miso
    # flipped_lsv_names = [str(name_str).split("#")[0] for name_str in flipped_miso_dict.keys()]
    with open('flipped_miso_thres_%.2f.txt' % flipped_thres, 'w') as flipped_file:
        flipped_file.write("name_rtpcr\tname_majiq\trtpcr_rest%f\tmajiq_rest\trtpcr_stim\tmajiq_stim\tunique_rtpcr_name?\tlsv_1_way?\ttarget_junction_coord\n")
        # with open(args.majiq_builder_file) as majiq_builder_file:
        #     majiq_builder = pickle.load(majiq_builder_file)
        #     for lsv in majiq_builder[1]:
        #         if lsv.id in flipped_lsv_names:
        #             for flip_key in flipped_majiq_dict:
        #                 if lsv.id in flip_key:
        #                     flipped_majiq_dict[flip_key] = flipped_majiq_dict[flip_key] + "\t%s\n" % str(lsv.junction_id[int(str(flip_key).split('#')[1])])
        for flipped_lsv in flipped_miso_dict:
            flipped_file.write(flipped_miso_dict[flipped_lsv])


if __name__ == '__main__':
    main()