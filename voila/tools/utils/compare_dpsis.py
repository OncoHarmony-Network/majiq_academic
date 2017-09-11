import itertools
from matplotlib import pyplot as plt
import numpy as np
import pylab as P
import scipy.stats

# Matthew Gazzara
# Made py-3 compatible by Caleb Radens


def plot_scatter_2expt(txt_f1, txt_f2, thresh=0.2, prob_thresh=0, display_ev_type='intersection', set1_name='set1',
                       set2_name='set2', title='Title', which_max=1):
    # dic1 = get_dpsi_dic(txt_f1, thresh = thresh)
    # dic2 = get_dpsi_dic(txt_f2, thresh = thresh)

    dic1 = get_dpsi_dic(txt_f1, thresh=0)
    dic2 = get_dpsi_dic(txt_f2, thresh=0)

    set1_ids = set(get_ids(txt_f1, thresh=thresh, prob=prob_thresh))
    set2_ids = set(get_ids(txt_f2, thresh=thresh, prob=prob_thresh))

    if display_ev_type == 'union':
        ev_list = set1_ids | set2_ids
    elif display_ev_type == 'intersection':
        ev_list = set1_ids & set2_ids
    else:
        raise RuntimeError(
            'incorrect overlapping event type...please choose union or intersection of LSVs between the conditions')

    plot_corr(dic1, dic2, ev_list, title=title, xax='dPSI %s' % set1_name, yax='dPSI %s' % set2_name, f_idx=0,
              which_max=which_max)
    plt.show()


def get_dpsi_dic(txt_file_path, thresh=0.2, restrict_list=False, CE_only=False, No_IR=False, out_file='None'):
    if not restrict_list == False:
        # res_IDs = set(open(restrict_list_file,'r').read().splitlines())
        res_IDs = set(restrict_list)
    dic = {}
    fd = open(txt_file_path, 'r')
    for line in fd:
        if line.startswith('#'): continue
        l = line.strip().split('\t')
        ID = l[2]
        lsv_type = l[7]
        if No_IR == True:
            if lsv_type.endswith('i'): continue
        if CE_only == True:
            if not lsv_type in ['s|1e1.1o1|1e2.1o1', 't|1e1.1o1|1e2.1o1']: continue
        if not restrict_list == False:
            if not ID in res_IDs: continue
        deltas = [abs(float(x)) for x in l[3].split(';')]
        if not max(deltas) >= thresh:
            continue
        try:
            dic[ID]
        except:
            dic[ID] = {}
            dic[ID]['info'] = {}
            juncs = l[16].split(';')
            junc_lens = [abs(int(x.split('-')[0]) - int(x.split('-')[1])) for x in juncs]
            dic[ID]['info']['junc_lens'] = junc_lens
            dic[ID]['info']['junc_coords'] = juncs
        # deltas = [abs(float(x)) for x in l[3].split(';')]
        delta_order = np.argsort(deltas)
        deltas = [float(x) for x in l[3].split(';')]
        psis1 = [float(x) for x in l[5].split(';')]
        psis2 = [float(x) for x in l[6].split(';')]
        # dic[ID][labels[n]] = {}
        dic[ID]['delta_order'] = delta_order
        dic[ID]['deltas'] = deltas
        dic[ID]['psis1'] = psis1
        dic[ID]['psis2'] = psis2
    fd.close()
    return dic


def get_ids(voila_file, thresh=0.2, prob=0):
    ids = []
    lines_num = 0
    try:
        fd = open(voila_file, 'r')
    except:
        print('voila file %s not found....returning zero changing IDs' % voila_file)
        return []
    for line in fd:
        if line.startswith('#'): continue
        lines_num += 1
        l = line.strip().split('\t')
        ID = l[2]
        deltas = [abs(float(x)) for x in l[3].split(';')]
        # deltas = []
        # for x in l[3].split(';'):
        #    dpsi = float(x)
        #    abs_dpsi = abs(dpsi)
        #    deltas.append(abs_dpsi)
        max_delta = max(deltas)
        probs = [float(x) for x in l[4].split(';')]
        max_prob = max(probs)

        if max_delta >= thresh:
            if max_prob >= prob:
                ids.append(ID)
    fd.close()
    print('for %s:\n\tlines processed: %s\n\tIDs added: %s' % (voila_file, lines_num, len(ids)))
    return ids


def changing_heatmaps(dic_1, dic_2, ev_list, conds_list=['cont1', 'exp1', 'cond2', 'exp2'], title='insert title',
                      out_file='output.txt'):
    skipped_evs = []
    data = []
    fw = open(out_file + '_PSI.txt', 'w')
    fw2 = open(out_file + '_dPSI.txt', 'w')
    fw.write('#Event')
    fw2.write('#Event\tComp1\tComp2\n')
    for c in conds_list:
        fw.write('\t%s' % c)
    fw.write('\n')
    for ev in ev_list:
        try:
            psi1_cont = dic_1[ev]['psis1']
            psi1_exp = dic_1[ev]['psis2']
            psi2_cont = dic_2[ev]['psis1']
            psi2_exp = dic_2[ev]['psis2']
        except:
            skipped_evs.append(ev)
            continue
        deltas1 = dic_1[ev]['deltas']
        deltas2 = dic_2[ev]['deltas']
        order1 = dic_1[ev]['delta_order']
        order2 = dic_2[ev]['delta_order']
        max_junc = order1[-1]
        row = [psi1_cont[max_junc], psi1_exp[max_junc], psi2_cont[max_junc], psi2_exp[max_junc]]
        data.append(row)
        fw.write('%s' % ev)
        for val in row:
            fw.write('\t%s' % val)
        fw.write('\n')
        fw2.write('%s\t%s\t%s\n' % (ev, deltas1[max_junc], deltas2[max_junc]))
    fw.close()
    fw2.close()


def plot_corr(dic_1, dic_2, ev_list, title='Title', xax='dPSI set1', yax='dPSI set2', f_idx=0, which_max=1):
    skipped_evs = []
    data1 = []
    data2 = []
    for ev in ev_list:
        try:
            deltas1 = dic_1[ev]['deltas']
            deltas2 = dic_2[ev]['deltas']
        except:
            skipped_evs.append(ev)
            continue
        order1 = dic_1[ev]['delta_order']
        order2 = dic_2[ev]['delta_order']
        if which_max == 1:
            max_junc = order1[-1]
        elif which_max == 2:
            max_junc = order2[-1]
        else:  # If you do not want to focus on 1 comp vs. the other...just pick the largest change observed in either
            maxJ1 = order1[-1]
            maxJ2 = order2[-1]
            abs_d1 = abs(deltas1[maxJ1])
            abs_d2 = abs(deltas2[maxJ2])
            if abs_d1 > abs_d2:
                max_junc = order1[-1]
            else:
                max_junc = order2[-1]
        data1.append(deltas1[max_junc])
        data2.append(deltas2[max_junc])
        # for n in range(len(deltas1)):
        #    data1.append(deltas1[n])
        #    data2.append(deltas2[n])

    plt.figure(f_idx)
    fig, ax = plt.subplots(1)
    ax.set_aspect('equal')
    plt.scatter(data1, data2, facecolors='#262626', edgecolors='#262626', alpha=0.95)
    plt.axvline(0, c='#262626', linewidth=1.5, alpha=0.9)
    plt.axhline(0, c='#262626', linewidth=1.5, alpha=0.9)
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    fit = plot_fit(data1, data2)
    plt.plot(data1, fit(data1), '#262626', linewidth=1.5, alpha=0.9)
    plt.xlabel(xax)
    plt.ylabel(yax)
    plt.title(title)

    print('%s total events input' % len(ev_list))
    print('%s events were not found in both comparisons...' % len(skipped_evs))


def plot_fit(x_data, y_data, return_stats=False):
    r = scipy.stats.pearsonr(x_data, y_data)
    fit = P.polyfit(x_data, y_data, 1)
    fit_fn = P.poly1d(fit)
    # P.plot(x_data,fit_fn(x_data),'#262626',linewidth=1.0,alpha=0.7)
    # P.annotate('Pearson\'s r:\n%.3f\nR^2:%.3f'%(r[0],r[0]*r[0]),fontsize=10, xy=(0.15, 0.9), xycoords='axes fraction')
    if return_stats == False:
        print(r, r[0] * r[0])
        return fit_fn
    if return_stats == True:
        return r[0], r[0] * r[0], r[1]


def pairwise_sample_heatmap(ordered_sample_names, voila_txt_directory, dpsi_thresh=0.2, prob_thresh=0,
                            voila_suffix='.txt'):
    vals = np.zeros(shape=(len(ordered_sample_names), len(ordered_sample_names)), dtype=np.int)
    combs = list(itertools.combinations(ordered_sample_names, 2))
    dpsi_prefixes = ['%s_%s' % (x[0], x[1]) for x in combs]
    alt_dpsi_prefixes = ['%s_%s' % (x[1], x[0]) for x in combs]
    for n in range(len(dpsi_prefixes)):
        dpsi = dpsi_prefixes[n]
        alt_dpsi = alt_dpsi_prefixes[n]

        print(dpsi)
        x = ordered_sample_names.index(combs[n][0])
        y = ordered_sample_names.index(combs[n][1])
        voila_dpsi_file = '%s/%s%s' % (voila_txt_directory, dpsi, voila_suffix)
        try:
            fd = open(voila_dpsi_file, 'r')
            fd.close()
        except:
            voila_dpsi_file = '%s/%s%s' % (voila_txt_directory, alt_dpsi, voila_suffix)
        changingIDs = get_ids(voila_dpsi_file, thresh=dpsi_thresh, prob=prob_thresh)
        vals[x, y] += len(changingIDs)
        vals[y, x] += len(changingIDs)
    plot_heatmap(vals, ordered_sample_names)


def plot_heatmap(vals, grps):
    # vals is np array: NxN. vals is num changing, vals2 is #complex
    from scipy.cluster.hierarchy import linkage
    from scipy.cluster.hierarchy import dendrogram
    row_clusters = linkage(vals, method='complete')
    row_dendr = dendrogram(row_clusters, labels=grps, count_sort=True)

    fig = plt.figure(0)
    ax = fig.add_subplot(111)

    df_rowclust = np.zeros(shape=vals.shape, dtype=np.int)
    df_rowclust2 = np.zeros(shape=vals.shape, dtype=np.float)

    outfp1 = open('./dPSI_pairs.txt', 'w+')

    header = "Sample\t"
    group_order = row_dendr['leaves']
    group_order = range(len(grps))
    for hh in group_order:
        header += "\t%s" % grps[hh]
    outfp1.write("%s\n" % header)
    for xid, xx in enumerate(reversed(group_order)):  # make sep. group orders
        outfp1.write('\t%s' % grps[xx])
        for yid, yy in enumerate(group_order):
            if xx == yy:
                t = 0.0
            # else:
            #    t = float(vals2[xx, yy])/vals[xx, yy]
            df_rowclust[xid, yid] = vals[xx, yy]
            # df_rowclust2[xid, yid] = t
            outfp1.write('\t%s' % df_rowclust[xid, yid])
        outfp1.write('\n')
    outfp1.close()
    cax = ax.matshow(df_rowclust, interpolation='nearest', cmap='Purples')
    for x in range(df_rowclust.shape[0]):
        for y in range(df_rowclust.shape[1]):
            if x == vals.shape[0] - y - 1:
                continue
            plt.text(y, x, '%d' % df_rowclust[x, y],
                     horizontalalignment='center',
                     verticalalignment='center',
                     )

    rev_df_grps = [grps[xx] for xx in group_order]
    df_grps = [grps[xx] for xx in reversed(group_order)]
    ax.set_xlim(-0.5, len(grps) - 0.5)
    ax.set_ylim(-0.5, len(grps) - 0.5)
    tks = ax.get_xticks()

    ticks = np.arange(0, len(grps))
    ax.set_xticks(ticks)
    ax.set_xticklabels(rev_df_grps)
    tks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels(df_grps)
    ax.set_title('Changing events')
    plt.tick_params(axis='both', which='both', top='off', left='off', right='off', bottom='off')
    fig.colorbar(cax)
    plt.show()
