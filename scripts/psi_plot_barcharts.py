from __future__ import division
from matplotlib import use

use('Agg')
from scripts import utils as utils_scripts
from collections import defaultdict
import analysis
from scripts.psi_scores import calculate_ead_simple, _save_or_show
import argparse
from pylab import *
from itertools import izip
import colorbrewer as cb

LSV_TYPES_DICT = {
    's|1e1.1|1e2.1':'SE',
    't|1e1.1|1e2.1':'SE',
    's|1e1.1|1e1.2':'A3SS',
    't|1e1.1|2e1.1':'A3SS',
    't|1e1.1|1e1.2':'A5SS',
    's|1e1.1|2e1.1':'A5SS'
}

BUILDER_OUT_FILE_TEMPLATE = "data/builder_output/toJuan.%s.old_majiq"

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)


def get_low_med_high_psis(psi1_met1, psi2_met1, psi1_met2, psi2_met2, low_thres = .1):
    low_med_high_psis_set = [[0, 0, 0, 0] , [0, 0, 0, 0]]

    better_in_majiq_mask = array(calculate_ead_simple(psi1_met1, psi2_met1)) < array(calculate_ead_simple(psi1_met2, psi2_met2))
    better_in_miso_mask = array(calculate_ead_simple(psi1_met1, psi2_met1)) > array(calculate_ead_simple(psi1_met2, psi2_met2))

    low_med_high_psis_set[0][0] = np.count_nonzero((array(psi1_met1) <= low_thres) & better_in_majiq_mask)
    low_med_high_psis_set[1][0] = np.count_nonzero((array(psi1_met1) <= low_thres) & better_in_miso_mask)
    low_med_high_psis_set[0][1] = np.count_nonzero((array(psi1_met1) < 1 - low_thres) & (array(psi1_met1) > low_thres) & better_in_majiq_mask)
    low_med_high_psis_set[1][1] = np.count_nonzero((array(psi1_met1) < 1 - low_thres) & (array(psi1_met1) > low_thres) & better_in_miso_mask)
    low_med_high_psis_set[0][2] = np.count_nonzero((array(psi1_met1) >= 1 - low_thres) & better_in_majiq_mask)
    low_med_high_psis_set[1][2] = np.count_nonzero((array(psi1_met1) >= 1 - low_thres) & better_in_miso_mask)
    low_med_high_psis_set[0][3] = np.count_nonzero(better_in_majiq_mask)
    low_med_high_psis_set[1][3] = np.count_nonzero(better_in_miso_mask)

    low_med_high_psi_diff = []
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[(array(psi1_met1) <= low_thres) ]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[(array(psi1_met2) <= low_thres) ])))
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[(array(psi1_met1) > low_thres) & (array(psi1_met1) < 1-low_thres) ]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[(array(psi1_met2) > low_thres) & (array(psi1_met2) < 1-low_thres) ])))
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[(array(psi1_met1) >= 1- low_thres)]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[(array(psi1_met2) >= 1- low_thres)])))
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2)))))

    return low_med_high_psis_set, len(psi1_met1), low_med_high_psi_diff


def get_low_med_high_cov(psi1_met1, psi2_met1, psi1_met2, psi2_met2, coverage_mat):
    low_med_high_psis_set = [[0, 0, 0, 0], [0, 0, 0, 0]]

    LOW_BOUND = 15
    HIGH_BOUND = 40

    better_in_majiq_mask = array(calculate_ead_simple(psi1_met1, psi2_met1)) < array(calculate_ead_simple(psi1_met2, psi2_met2))
    better_in_miso_mask = array(calculate_ead_simple(psi1_met1, psi2_met1)) > array(calculate_ead_simple(psi1_met2, psi2_met2))

    low_med_high_psis_set[0][0] = np.count_nonzero([(coverage_mat <= LOW_BOUND) & (better_in_majiq_mask)])
    low_med_high_psis_set[1][0] = np.count_nonzero([(coverage_mat <= LOW_BOUND) & better_in_miso_mask])
    low_med_high_psis_set[0][1] = np.count_nonzero([(coverage_mat >= LOW_BOUND) & (coverage_mat < HIGH_BOUND) & better_in_majiq_mask])
    low_med_high_psis_set[1][1] = np.count_nonzero([(coverage_mat >= LOW_BOUND) & (coverage_mat < HIGH_BOUND) & better_in_miso_mask])
    low_med_high_psis_set[0][2] = np.count_nonzero([(coverage_mat >= HIGH_BOUND) & better_in_majiq_mask])
    low_med_high_psis_set[1][2] = np.count_nonzero([(coverage_mat >= HIGH_BOUND) & better_in_miso_mask])
    low_med_high_psis_set[0][3] = np.count_nonzero(better_in_majiq_mask)
    low_med_high_psis_set[1][3] = np.count_nonzero(better_in_miso_mask)

    low_med_high_psi_diff = []
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[coverage_mat <= LOW_BOUND]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[(coverage_mat <= LOW_BOUND)])))
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[(coverage_mat > LOW_BOUND) & (coverage_mat < HIGH_BOUND)]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[(coverage_mat >= LOW_BOUND) & (coverage_mat < HIGH_BOUND)])))
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[coverage_mat >= HIGH_BOUND]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[coverage_mat >= HIGH_BOUND])))
    low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2)))))

    return low_med_high_psis_set, len(psi1_met1), low_med_high_psi_diff

def get_low_med_high_cov_psi(psi1_met1, psi2_met1, psi1_met2, psi2_met2, coverage_mat):
    low_med_high_psis_set = [
        [0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0],
        [0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0]
    ]
    low_med_high_psi_diff = []

    LOW_THRES = .1
    tmp_sorted_coverage = np.array(coverage_mat)
    tmp_sorted_coverage.sort()
    LOW_BOUND = tmp_sorted_coverage[coverage_mat.size/3]
    HIGH_BOUND = tmp_sorted_coverage[2*coverage_mat.size/3]
    print "LOW BOUND: %d, HIGH_BOUND %d" % (LOW_BOUND, HIGH_BOUND)
    print coverage_mat.size, len(tmp_sorted_coverage[:coverage_mat.size/3]), len(tmp_sorted_coverage[coverage_mat.size/3: 2*coverage_mat.size/3]), len(tmp_sorted_coverage[2*coverage_mat.size/3:])
    COV_CONDS = [
        (coverage_mat < LOW_BOUND),
        (coverage_mat >= LOW_BOUND) & (coverage_mat < HIGH_BOUND),
        (coverage_mat >= HIGH_BOUND),
        True
    ]
    PSI_COND = [
        (array(psi1_met1) <= LOW_THRES),
        (array(psi1_met1) < 1 - LOW_THRES) & (array(psi1_met1) > LOW_THRES),
        (array(psi1_met1) >= (1- LOW_THRES)),
        True
    ]

    better_in_majiq_mask = array(calculate_ead_simple(psi1_met1, psi2_met1)) < array(calculate_ead_simple(psi1_met2, psi2_met2))
    better_in_miso_mask = array(calculate_ead_simple(psi1_met1, psi2_met1)) > array(calculate_ead_simple(psi1_met2, psi2_met2))

    for i, cov_cond in enumerate(COV_CONDS):
        for j, psi_cond in enumerate(PSI_COND):
            low_med_high_psis_set[0][i*len(COV_CONDS)+j] = np.count_nonzero([cov_cond & psi_cond & better_in_majiq_mask])
            low_med_high_psis_set[1][i*len(COV_CONDS)+j] = np.count_nonzero([cov_cond & psi_cond & better_in_miso_mask])
            low_med_high_psi_diff.append(-(np.mean(array(calculate_ead_simple(psi1_met1, psi2_met1))[cov_cond & psi_cond]) - np.mean(array(calculate_ead_simple(psi1_met2, psi2_met2))[cov_cond & psi_cond])))
            if np.count_nonzero(cov_cond & psi_cond) == 1:
                sys.stdout.write(str(coverage_mat.size)+"\t")
            else:
                sys.stdout.write(str(np.count_nonzero(cov_cond & psi_cond))+"\t")
        sys.stdout.write("\n")
    sys.stdout.flush()

    return low_med_high_psis_set, len(psi1_met1), low_med_high_psi_diff, (COV_CONDS[1] & better_in_miso_mask)


def plot_delta_bars_perc_all_in_grid(cov_psi_ndmatrix, majiq_num_events_list, majiq_diff_ndmatrix, replica_names, plotpath=None, alpha=0.5):

    replica_names_joint = '; '.join(["%s%s" % (name1, name2) for name1, name2 in grouped(replica_names, 2)])
    plotname="MAJIQ_Vs_MISO_delta_expected_psi_grid_alpha_%.1f. \nReplicates %s" % (float(alpha), replica_names_joint)

    print cov_psi_ndmatrix
    print np.reshape(np.repeat(np.array(majiq_num_events_list), cov_psi_ndmatrix.size/len(majiq_num_events_list)), cov_psi_ndmatrix.shape)
    psi_ndmatrix_perc = cov_psi_ndmatrix/np.reshape(np.repeat(np.array(majiq_num_events_list), cov_psi_ndmatrix.size/len(majiq_num_events_list)), cov_psi_ndmatrix.shape)

    means_psis_meth1 = np.apply_along_axis( np.mean, axis=2, arr=psi_ndmatrix_perc.T)[:, 0]
    means_psis_meth2 = np.apply_along_axis( np.mean, axis=2, arr=psi_ndmatrix_perc.T)[:, 1]

    stds_psis_meth1 = np.apply_along_axis( np.std, axis=2, arr=psi_ndmatrix_perc.T)[:, 0]
    stds_psis_meth2 = np.apply_along_axis( np.std, axis=2, arr=psi_ndmatrix_perc.T)[:, 1]

    n_groups = 4
    fig, axxs = plt.subplots(n_groups, n_groups, sharex=True, sharey=True, figsize=[12, 12], dpi=300)

    for ii in xrange(n_groups):
        for jj in xrange(n_groups):
            ax1 = axxs[ii][jj]

            index = np.arange(1)
            bar_width = 0.35
            # bar_width = 0.8

            opacity = 0.4
            error_config = {'ecolor': '0.3'}
            ax1.bar(index, np.array([means_psis_meth1[ii*n_groups+jj]-means_psis_meth2[ii*n_groups+jj]]), bar_width,
                     alpha=opacity,
                     color='b',
                     yerr=stds_psis_meth1[ii*n_groups+jj]-stds_psis_meth2[ii*n_groups+jj],
                     error_kw=error_config,
                     label='MAJIQ')
                     # label='Better in MAJIQ')
            # ax1.set_xlabel('PSIs values')
            # if is_coverage:
            #     ax1.set_xlabel('Coverage')
            if jj==0:
                ax1.set_ylabel('% of improved Events (LSVs)')
            ax1.set_ylim([-.1, .25])
            ax2 = ax1.twinx()
            ax2.bar(index + bar_width, np.mean(majiq_diff_ndmatrix, axis=0)[ii*n_groups+jj], bar_width,
                     alpha=opacity,
                     color='r',
                     yerr=np.std(majiq_diff_ndmatrix, axis=0)[ii*n_groups+jj],
                     error_kw=error_config,
                     label='Better in MISO')

            if jj == n_groups-1:
                ax2.set_ylabel('-E(Delta(Delta(PSI)))')
            else:
                ax2.axes.get_yaxis().set_ticks([])
            ax2.set_ylim([-.1, .25])
            plt.xlabel('PSIs values')


    # plt.ylabel('% of improved # Events (LSVs)')
    # plt.title(plotname)
    # plt.title('Delta in expected PSI between %s.' % (replica_names_joint))
    # plt.xticks(index + bar_width, psi_ranges)
    # plt.legend(loc=2)

    plt.tight_layout()
    _save_or_show(plotpath, plotname.replace('\n', ' - ') + '_percentages')


def plot_delta_bars_percentages(psi_ndmatrix, majiq_num_events_list, majiq_diff_ndmatrix, replica_names, plotpath=None, is_coverage=False):

    replica_names_joint = '; '.join(["%s%s" % (name1, name2) for name1, name2 in grouped(replica_names, 2)])
    plotname="MAJIQ_Vs_MISO_delta_expected_psi. \nReplicates %s" % (replica_names_joint)

    psi_ranges = ['low', 'med', 'high', 'all']
    psi_ndmatrix_perc = psi_ndmatrix/np.reshape(np.repeat(np.array(majiq_num_events_list), psi_ndmatrix.size/len(majiq_num_events_list)), psi_ndmatrix.shape)

    means_psis_meth1 = np.apply_along_axis( np.mean, axis=2, arr=psi_ndmatrix_perc.T)[:, 0]
    means_psis_meth2 = np.apply_along_axis( np.mean, axis=2, arr=psi_ndmatrix_perc.T)[:, 1]

    stds_psis_meth1 = np.apply_along_axis( np.std, axis=2, arr=psi_ndmatrix_perc.T)[:, 0]
    stds_psis_meth2 = np.apply_along_axis( np.std, axis=2, arr=psi_ndmatrix_perc.T)[:, 1]


    n_groups = 4


    fig, ax1 = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.35
    # bar_width = 0.8

    opacity = 0.4
    error_config = {'ecolor': '0.3'}

    ax1.bar(index, means_psis_meth1-means_psis_meth2, bar_width,
             alpha=opacity,
             color='b',
             yerr=stds_psis_meth1-stds_psis_meth2,
             error_kw=error_config,
             label='MAJIQ')
             # label='Better in MAJIQ')
    ax1.set_xlabel('PSIs values')
    if is_coverage:
        ax1.set_xlabel('Coverage')
    ax1.set_ylabel('% of improved Events (LSVs)')
    ax1.set_ylim([-.1, .35])
    ax2 = ax1.twinx()
    ax2.bar(index + bar_width, np.mean(majiq_diff_ndmatrix, axis=0), bar_width,
             alpha=opacity,
             color='r',
             yerr=np.std(majiq_diff_ndmatrix, axis=0),
             error_kw=error_config,
             label='Better in MISO')

    ax2.set_ylabel('E(Delta(Delta(PSI)))')
    ax2.set_ylim([-.1, .35])
    # plt.xlabel('PSIs values')
    # plt.ylabel('% of improved # Events (LSVs)')
    plt.title(plotname)
    # plt.title('Delta in expected PSI between %s.' % (replica_names_joint))
    plt.xticks(index + bar_width, psi_ranges)
    # plt.legend(loc=2)

    plt.tight_layout()
    _save_or_show(plotpath, plotname.replace('\n', ' - ') + '_percentages')



def plot_delta_expected_method1Vsmethod2(psi_list_data_set, replica_names, plotpath=None):

    replica_names_joint = '; '.join(["%s%s" % (name1, name2) for name1, name2 in grouped(replica_names, 2)])
    plotname="MAJIQ_Vs_MISO_delta_expected_psi. \nReplicates %s" % (replica_names_joint)

    psi_ranges = ['low', 'med', 'high', 'all']

    means_psis_meth1 = np.apply_along_axis( np.mean, axis=2, arr=array(psi_list_data_set).T)[:,0] #psi_list_data_set[0][0]
    means_psis_meth2 = np.apply_along_axis( np.mean, axis=2, arr=array(psi_list_data_set).T)[:,1] #psi_list_data_set[1][1]

    stds_psis_meth1 = np.apply_along_axis( np.std, axis=2, arr=array(psi_list_data_set).T)[:,0]
    stds_psis_meth2 = np.apply_along_axis( np.std, axis=2, arr=array(psi_list_data_set).T)[:,1]


    n_groups = 4


    fig, ax = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.35

    opacity = 0.4
    error_config = {'ecolor': '0.3'}

    plt.bar(index, means_psis_meth1, bar_width,
             alpha=opacity,
             color='b',
             yerr=stds_psis_meth1,
             error_kw=error_config,
             label='Better in MAJIQ')

    plt.bar(index + bar_width, means_psis_meth2, bar_width,
             alpha=opacity,
             color='r',
             yerr=stds_psis_meth2,
             error_kw=error_config,
             label='Better in MISO')

    plt.xlabel('PSIs values')
    plt.ylabel('# Events (LSVs)')
    plt.title(plotname)
    # plt.title('Delta in expected PSI between %s.' % (replica_names_joint))
    plt.xticks(index + bar_width, psi_ranges)
    plt.legend(loc=2)

    plt.tight_layout()
    _save_or_show(plotpath, plotname.replace('\n', ' - '))


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def plot_delta_expected_majiq_others(psi_dict_lists, replica_names, plotpath=None):

    colors_dict = {
        'miso': 'blue'
    }
    colors_list=['blue', 'red', 'green', 'purple', 'orange', 'yellow']
    replica_names_joint = '; '.join(["%s%s" % (name1, name2) for name1, name2 in grouped(replica_names, 2)])
    plotname="MAJIQ_Vs_Others_delta_expected_psi. \nReplicates %s" % (replica_names_joint)


    fig = plt.figure(figsize=[8, 6]) # In inches
    #figure out how many groups of events exist
    max_difference=0
    font = {'size': 10} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)

    lthres = [0.05, 0]
    for thres in lthres:
        win_cdf_all = []
        win_psi_all = []
        lose_cdf_all = []
        lose_psi_all = []
        nbins=5000

        for met_key, met_diff_list in psi_dict_lists.iteritems():
            for jj, methods_diff in enumerate(met_diff_list):
                palette=cb.Blues[7][::-1]

                # methods_diff=np.mean(met_diff_list)
                methods_diff.sort()
                win_elems=methods_diff[np.where(methods_diff>thres)]
                uu, ii = np.unique(win_elems, return_inverse=True)

                win_cdf=[0]
                for w_freq in np.bincount(ii):
                    win_cdf.append(win_cdf[-1]+((1.*w_freq)/len(win_elems)))
                win_cdf_final = np.dot(win_cdf[1:], len(win_elems)/len(methods_diff)) # np.dot(win_cdf[1:],len(win_elems))

                aux=[]
                aux_psi=[]
                for aa in xrange(0,nbins,1):
                    aux.append(win_cdf_final[int((aa/(nbins*1.0))*len(win_cdf_final))])
                    aux_psi.append(uu[int((aa/(nbins*1.0))*len(uu))])
                win_cdf_all.append(aux)
                win_psi_all.append(aux_psi)

                # lose_elems=methods_diff[np.where(methods_diff<0)][::-1]
                lose_elems=-methods_diff[np.where(methods_diff<-thres)]
                uul, iil = np.unique(lose_elems, return_inverse=True)

                _cdf=[0]
                for w_freq in np.bincount(iil):
                    _cdf.append(_cdf[-1]+((1.*w_freq)/len(lose_elems)))
                lose_cdf = np.dot(_cdf[1:],len(lose_elems)/len(methods_diff))  # np.dot(_cdf[1:][::-1],-len(lose_elems))

                aux=[]
                aux_psi=[]
                for aa in xrange(0,nbins,1):
                    aux.append(lose_cdf[int((aa/(nbins*1.0))*len(lose_cdf))])
                    aux_psi.append(uul[int((aa/(nbins*1.0))*len(uul))])
                lose_cdf_all.append(aux)
                lose_psi_all.append(aux_psi)


                max_difference = max(max_difference, max(len(win_elems), len(lose_elems)))

        lbl_majiq="MAJIQ winning |ddPSI|>%.0f%% - mean %.0f%%" % (100*thres, 100*np.mean(win_cdf_all, axis=0)[-1])
        lbl_miso="MISO winning |ddPSI|>%.0f%% - mean %.0f%%" % (100*thres, 100*np.mean(lose_cdf_all, axis=0)[-1])
        _ls = '-'
        if thres:
            _ls = 'dashed'
            plt.axvline(x=0.05, ymin=0, ymax=1, linewidth=2, color='grey', ls='dashed')

        # MAJIQ wins
        plt.plot(np.append(np.mean(win_psi_all, axis=0), 1),
                 np.append(np.mean(win_cdf_all, axis=0), np.mean(win_cdf_all, axis=0)[-1]),
                 # label="MAJIQ wins - Mean (%.2f) & STD (%.2f)" % (np.mean(win_cdf_all, axis=0)[-1], np.std(win_cdf_all, axis=0)[-1]),
                 label=lbl_majiq,
                 color=rgb_to_hex(cb.Blues[3][-1]), lw=2, ls=_ls)

        # MISO wins
        plt.plot(np.append(np.mean(lose_psi_all, axis=0), 1),
                 np.append(np.mean(lose_cdf_all, axis=0), np.mean(lose_cdf_all, axis=0)[-1]),
                 # label="MISO wins - Mean (%.2f) & STD (%.2f)" % (np.mean(lose_cdf_all, axis=0)[-1], np.std(lose_cdf_all, axis=0)[-1]),
                 label=lbl_miso,
                 color=rgb_to_hex(cb.Reds[3][-1]), lw=2, ls=_ls)

        if not thres:
            plt.fill_between(np.append(np.mean(win_psi_all, axis=0), 1),
                             np.append(np.mean(win_cdf_all, axis=0)+np.std(win_cdf_all, axis=0), np.mean(win_cdf_all, axis=0)[-1]+np.std(win_cdf_all, axis=0)[-1]),
                             np.append(np.mean(win_cdf_all, axis=0)-np.std(win_cdf_all, axis=0), np.mean(win_cdf_all, axis=0)[-1]-np.std(win_cdf_all, axis=0)[-1]),
                             facecolor=rgb_to_hex(cb.Blues[3][-1]),
                             alpha=0.4,
                             linewidth=0.0)

            plt.fill_between(np.append(np.mean(lose_psi_all, axis=0), 1),
                             np.append(np.mean(lose_cdf_all, axis=0)+np.std(lose_cdf_all, axis=0), np.mean(lose_cdf_all, axis=0)[-1]+np.std(lose_cdf_all, axis=0)[-1]),
                             np.append(np.mean(lose_cdf_all, axis=0)-np.std(lose_cdf_all, axis=0), np.mean(lose_cdf_all, axis=0)[-1]-np.std(lose_cdf_all, axis=0)[-1]),
                            facecolor=rgb_to_hex(cb.Reds[3][-1]),
                            alpha=0.4,
                            linewidth=0.0)

    plt.xlabel("Delta Delta PSI", fontsize=11)
    plt.ylabel("Number of LSVs", fontsize=11)
    # plt.xlim(-1, 1)
    # plt.ylim(-(max_difference+1000), max_difference+1000) # TODO: revise
    #
    # plt.plot([-1,1], [0,0], color='black', lw=1, alpha=.5)
    # plt.plot([0,0], [-(max_difference+1000), max_difference+1000], color='black', lw=1, alpha=.5)

    plt.xlim(0, .45)
    plt.ylim(0, 1)

    plt.title(plotname, fontsize=13)
    plt.legend(loc=2, fontsize=11)
    plt.tight_layout()
    _save_or_show(plotpath, plotname.replace('\n', ' - '))



def intersect_sets(majiq1, majiq2):
    """Find common names and return their psi values in 2 lists"""
    names1 = [m[1] for m in majiq1[1]]
    names2 = [m[1] for m in majiq2[1]]
    common_names = set(names1).intersection(set(names2))
    return np.array(majiq1[0])[np.array([name in common_names for name in names1])], np.array(majiq2[0])[np.array([name in common_names for name in names2])], np.array(majiq1[1])[np.array([name in common_names for name in names1])], np.array(majiq2[1])[np.array([name in common_names for name in names2])]


def main():
    """
    Script for testing MAJIQ against other algorithms for PSIs across several replicates

    - All LSVs 'ways' are equally considered. The order of each way within a LSV should be the same in MAJIQ and MISO.
    """

    #python ~/Projects/old_majiq/scripts/psi_plot_barcharts.py --old_majiq old_majiq/ --miso miso/ --names Hip1 Hip2 Hip1 Hip4 Hip5 Hip6 Liv1 Liv2 Liv1 Liv4 Liv4 Liv5 --plotpath ./output/
    parser = argparse.ArgumentParser()
    parser.add_argument('--old_majiq', dest='majiq_dir', type=str, help='Path for MAJIQ psi pickles to evaluate')
    parser.add_argument('--miso', dest='miso_dir', type=str,  help='Path for MISO psi pickles to evaluate')
    parser.add_argument('--nb', dest='nb_dir', required=False, type=str,  help='Path for Naive Bootstrapping psi pickles to evaluate')
    parser.add_argument('--names', dest='rep_names', nargs='+', required=True, help='Replicate names used to identify each pair [NOTE that the order in which the names are provided defines the pairs]')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--builder-files-dir', dest='builder_files_dir')
    parser.add_argument('--cached-cov', dest='cached_cov')
    parser.add_argument('--alpha', default=0.5, help='Alpha used in MAJIQ (for naming purposes only).')
    args = parser.parse_args()

    majiq_better_worse_list = []
    majiq_num_events_list = []
    majiq_low_med_high_diff_list = []
    majiq_coverage_list = []
    better_worse_dict = defaultdict(list)

    if args.cached_cov:
        majiq_coverage_list = pickle.load(open(args.cached_cov, 'r'))

    majiq_files = []
    miso_files = []
    foo_dict = {
        'Hip1': 'Hippocampus1',
        'Hip2': 'Hippocampus2',
        'Liv1': 'Liver1',
        'Liv3': 'Liver3',
        'Hip4': 'Hippocampus4',
        'Hip5': 'Hippocampus5',
        'Liv4': 'Liver4',
        'Liv5': 'Liver5',
        'Spl1': 'Spleen1',
        'Spl2': 'Spleen2',

    }
    for rep_name in args.rep_names:
        # majiq_files.extend(utils_scripts.list_files_or_dir([args.majiq_dir], suffix='majiq_psi.pickle', containing=foo_dict[rep_name])) #psigroup.pickle #'psigroup.pickle'
        majiq_files.extend(utils_scripts.list_files_or_dir([args.majiq_dir], suffix='psigroup.pickle', containing=rep_name)) #psigroup.pickle #
        miso_files.extend(utils_scripts.list_files_or_dir([args.miso_dir], suffix='miso_summary', containing=rep_name))

    for ii, majiq_file in enumerate(majiq_files):
        if ii % 2:
            majiq1 = pickle.load(open(majiq_files[ii-1]))
            majiq2 = pickle.load(open(majiq_file))
            print "Events in MAJIQ: rep1=%d; rep2=%d" % (len(majiq1[1]), len(majiq2[1]))
            psi_met1_rep1, psi_met1_rep2, majiq1_names, majiq2_names = intersect_sets(majiq1, majiq2)
            # psi_met1_rep2 = psivalues[0]
            psi_names_met1 = defaultdict()

            print "Events after intersection in MAJIQ: rep1=%d; rep2=%d" % (len(psi_met1_rep1), len(psi_met1_rep2))
            # Discard LSVs with only one PSI
            for i, psis_lsv_met1 in enumerate(psi_met1_rep1):
                if len(psis_lsv_met1) > 1 or len(psi_met1_rep2[i]) > 1:
                    continue
                # if majiq2[1][i][2] not in LSV_TYPES_DICT.keys():
                #     continue
                psi_names_met1[majiq2_names[i][1]] = i

            # Method1 (MAJIQ) psi scores
            psi_list1_met1 = []
            psi_list2_met1 = []

            psi_lists_met2 = []

            majiq_vs = 'miso'
            miso_all = []
            for miso_file_index in [ii-1, ii]:
                miso_file = miso_files[miso_file_index]
                miso_psis_dict = defaultdict()
                num_events = 0
                with open(miso_file, 'r') as miso_res:
                    for miso_line in miso_res:
                        if miso_line.startswith("event"): continue
                        num_events+=1
                        miso_fields = miso_line.split('\t')
                        if miso_fields[0] not in psi_names_met1:
                            continue
                        miso_psis_dict[miso_fields[0]] = miso_fields[1]

                print "Events in MISO: rep%d=%d" % ((miso_file_index % 2) + 1, num_events)
                miso_all.append(miso_psis_dict)

            print "Events in MISO after intersection with MAJIQ: rep1=%d; rep2=%d" % (len(miso_all[0]), len(miso_all[1]))
            miso_common_names = set(miso_all[0].keys()).intersection(miso_all[1].keys())
            print "Events common in MISO: %d" % len(miso_common_names)

            for miso_psis_dict in miso_all:
                miso_psis = []
                for psi_name in sorted(psi_names_met1.keys()):
                    if psi_name not in miso_common_names:
                        print "%s is not in all MISO replicates" % psi_name
                        del psi_names_met1[psi_name]
                        continue
                    try:
                        if len(miso_psis_dict[psi_name].split(",")) > 1:
                            print "[WARNING] %s LSV is multiway, shouldn't be here..." % psi_name
                        miso_psis_values = [float(miso_psi) for miso_psi in miso_psis_dict[psi_name].split(",")]
                    except KeyError, e:
                        print "LSV %s is in MAJIQ but not in MISO!" % e
                        del psi_names_met1[psi_name]
                        continue
                        # if len(miso_psis_values) < 2:
                        #     del majiq_psi_names[psi_name]
                        #     continue

                    # if len(miso_psis_values) == 1:
                    #     miso_psis_values.append(1.0 - miso_psis_values[0])
                    miso_psis.extend(miso_psis_values)

                psi_lists_met2.append(miso_psis)

            print "Events after intersection of MAJIQ and MISO: %d" % len(psi_names_met1.keys())
            for psi_name in sorted(psi_names_met1.keys()):
                # for j, psi_lsv in enumerate(psi_met1_rep1[psi_names_met1[psi_name]]):
                #     psi_list1_met1.append(sum(psi_lsv*analysis.psi.BINS_CENTER))
                #     psi_list2_met1.append(sum(psi_met1_rep2[psi_names_met1[psi_name]][j]*analysis.psi.BINS_CENTER))
                psi_list1_met1.append(sum(psi_met1_rep1[psi_names_met1[psi_name]][0]*analysis.psi.BINS_CENTER))
                psi_list2_met1.append(sum(psi_met1_rep2[psi_names_met1[psi_name]][0]*analysis.psi.BINS_CENTER))

            # lo_me_hi_mat, num_events, lo_med_high_diff, suspicous_guys = get_low_med_high_cov_psi(psi_list1_met1, psi_list2_met1, psi_lists_met2[0], psi_lists_met2[1], np.array(coverage_list))
            # lo_me_hi_mat, num_events, lo_med_high_diff = get_low_med_high_cov(psi_list1_met1, psi_list2_met1, psi_lists_met2[0], psi_lists_met2[1], np.array(coverage_list))
            # lo_me_hi_mat, num_events, lo_med_high_diff = get_low_med_high_psis(psi_list1_met1, psi_list2_met1, psi_lists_met2[0], psi_lists_met2[1])
            # pickle.dump(np.array(psi_names_met1.keys())[suspicous_guys], open('suspicious.pkl', 'w'))
            # for sus_guy in np.array(psi_names_met1.keys())[suspicous_guys]:
            #     print str(sus_guy).split(':')[0]

            better_worse_dict[majiq_vs].append(-calculate_ead_simple(psi_list1_met1, psi_list2_met1) +  calculate_ead_simple(psi_lists_met2[0], psi_lists_met2[1]))
            # majiq_better_worse_list.append(lo_me_hi_mat)
    #         majiq_num_events_list.append(num_events)
    #         majiq_low_med_high_diff_list.append(lo_med_high_diff)
    #
    # if not args.cached_cov:
    #     pickle.dump(majiq_coverage_list, open('coverage_tmp.pickle', 'w'))

    plot_delta_expected_majiq_others(better_worse_dict, args.rep_names, args.plotpath)
    # plot_delta_expected_method1Vsmethod2(majiq_better_worse_list, args.rep_names, args.plotpath)
    # plot_delta_bars_percentages(np.array(majiq_better_worse_list), majiq_num_events_list, np.array(majiq_low_med_high_diff_list), args.rep_names, args.plotpath, is_coverage=True)
    # plot_delta_bars_perc_all_in_grid(np.array(majiq_better_worse_list), majiq_num_events_list, np.array(majiq_low_med_high_diff_list), args.replica_names, plotpath=args.plotpath, alpha=args.alpha)

if __name__ == '__main__':
    main()