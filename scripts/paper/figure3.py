# -*- coding: utf-8 -*-
import sys
grimoire_path = '/Users/Jordi/software/majiq'
if grimoire_path not in sys.path:
    sys.path.append(grimoire_path)

from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.image as mpimg
import cPickle as pickle
from os import listdir
from os.path import isfile, join


fidx = 0
lims = [500, 1000]

BREWER_PALETTE = [(228, 26, 28),
                  (55, 126, 184),
                  (77, 175, 74),
                  (152, 78, 163),
                  (255, 127, 0)]


def autolabel(rects, ax, msg, size, b=False):
    # attach some text labels
    for ridx, rect in enumerate(rects):
        width = rect.get_width()
        x = rect.get_x()

        if b:
            x += width/2.
        else:
            x -= width/2.

        for midx, mm in enumerate(msg[ridx]):
            ax.text(x, 20000 + midx*200 - ridx * 1000, '%s' % mm,
                    ha='left', va='bottom', fontsize=size)


def autolabel_horz(rects, ax, offs=8, msg='', size=8):
    # attach some text labels
    for idx in np.arange(0, len(rects), offs+1):
        rect = rects[idx]
        width = rect.get_width()
        ypos = [+4, +2]
        for midx, mm in enumerate(msg[idx]):
            x = rect.get_x()+width + 150
            y = idx+ypos[midx % 2]
            if midx >= 2:
                x += 200
            ax.text(x, y, mm, ha='center', va='bottom', fontsize=size)


def print_message(bars, extrabars=None):
    msg = []
    nbars = bars / bars.sum()
    if not extrabars is None:
        nextrabars = extrabars/extrabars.sum()
    for bix, bb in enumerate(bars):
        total = "Abs: %s" % bb
        norm = "Nrm: %.2e" % nbars[bix]
        p1 = "P(X<=S): %.2e" % np.sum(nbars[:bix+1])
        p2 = "P(X>=S): %.2e" % np.sum(nbars[bix:])
        if not extrabars is None:
            total += "| %s" % extrabars[bix]
            norm += "| %.2e" % nextrabars[bix]
            p1 += "| %.2e" % np.sum(nextrabars[:bix+1])
            p2 += "| %.2e" % np.sum(nextrabars[bix:])
        msg.append((total, norm, p1, p2))
    return msg


def plot_fdrheatmap(vals, grps):

    global fidx

    fig, ax = pyplot.subplots(1, 1)
    pyplot.pcolor(vals, cmap='bwr')
    ax.set_xticklabels(grps)
    ax.set_yticklabels(grps)
    tks = ax.get_xticks()
    ax.set_xticks(tks+1)
    tks = ax.get_yticks()
    ax.set_yticks(tks+0.5)
    pyplot.show()
    fidx += 1


def plot_countings(vals, typs, counts=[0,p]):

    global fidx

    offs = [0, 0.8, 1.60, 2.40]

    pyplot.figure(fidx)
    junc_bars = [0] * 20
    ex_bars = [0] * 20
    bars_s = [0] * 20
    bars_t = [0] * 20

    for idx, tt in enumerate(typs):
        juncs = [xx for xx in tt.split('|')[1:] if not xx.endswith('0')]
        nums = len(juncs)
        junc_bars[nums] += vals[idx]

        exnum = []
        ssnum = []
        for jj in juncs:
            exnum.append(int(jj.split('e')[1].split('.')[0]))
            ssnum.append(int(jj.split('e')[0]))

        nums_ex = max(exnum)
        ex_bars[nums_ex] += vals[idx]

        if tt[0] == 's':
            bars = bars_s
        else:
            bars = bars_t

        nums_ss = max(ssnum)
        bars[nums_ss] += vals[idx]

    dummy_ex = [idx for idx, xx in enumerate(junc_bars) if xx > 0]
    dummy_junc = [idx for idx, xx in enumerate(ex_bars) if xx > 0]
    dummys = [idx for idx, xx in enumerate(bars_s) if xx > 0]
    dummyt = [idx for idx, xx in enumerate(bars_t) if xx > 0]

    lim_total = max([max(dummys), max(dummyt), max(dummy_junc), max(dummy_ex)])
    x = np.arange(1, 4*lim_total + 1, 4)

    lim = max(dummy_junc)
    junc_bars = np.array(junc_bars, dtype=float)
    nbars = junc_bars / junc_bars.sum()
    col = [xx/float(255) for xx in BREWER_PALETTE[0]]
    pyplot.bar(x[:lim] + offs[0], nbars[1:lim+1], edgecolor="none", label="#LSV per #juncs",
               width=0.8, color=col, alpha=0.5)

    lim = max(dummy_ex)
    ex_bars = np.array(ex_bars, dtype=float)
    nbars = ex_bars / ex_bars.sum()
    col = [xx/float(255) for xx in BREWER_PALETTE[1]]
    pyplot.bar(x[:lim] + offs[1], nbars[1:lim+1], edgecolor="none", label="#LSV per #exons", width=0.8,
               color=col, alpha=0.5)

    lim = max(max(dummys), max(dummyt))
    bars_s = np.array(bars_s, dtype=float)
    nbars_s = bars_s/bars_s.sum()
    bars_t = np.array(bars_t, dtype=float)
    nbars_t = bars_t/bars_t.sum()
    pyplot.bar(x[:lim] + offs[2], nbars_s[1:lim+1], label="5' splicesites in Source LSV", edgecolor="none",
               color=(1.0, 0.5, 0.0), width=0.8, alpha=0.5)
    pyplot.bar(x[:lim] + offs[3], nbars_t[1:lim+1], label="3' splicesites in Target LSV", edgecolor="none",
               color=(0.92, 0.0, 0.55), width=0.8, alpha=0.5)

    pyplot.ylabel(' # datapoints')
    pyplot.xlabel(' # elements')
    pyplot.title('Total Number of datapoints %s ( %d different LSVs)' % (counts[0], counts[1]))
    pyplot.legend(loc='best')
    pyplot.show()
    fidx += 1


def plot_lsv_types_hist(vals, typs, img_path=None, lim_val=None, extra_title=""):

    global fidx

    fig = pyplot.figure(fidx)

    if not img_path is None:
        offs = 6
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 10])
        gs.update(left=0.01, right=0.99, hspace=0.00, wspace=0.0001)
        gs0 = gridspec.GridSpecFromSubplotSpec(len(vals), 1, subplot_spec=gs[0], wspace=0.0, hspace=0.05)
        thumb = []
        for ii, yy in enumerate(gs0):
            thumb.append(fig.add_subplot(yy))
        pl = fig.add_subplot(gs[1])

        vals2 = []

        for ii, xx in enumerate(vals):
            vals2.append(xx)
            for ii2 in range(offs):
                vals2.append(0)

        lsv_typ = np.array(range(len(vals2)))
              
        bars_o = pl.barh(lsv_typ, vals2, height=4, edgecolor="none")
        
        pl.set_ylim((0, len(vals2)))
        pl.set_yticks([])
        msg = print_message(np.array(vals2, dtype=np.float))
        autolabel_horz(bars_o, pl, offs, msg)
        for xidx, xx in enumerate(reversed(typs)):
            bar_idx = (len(typs)-xidx-1) * (offs+1)

            if xx[0] == 's':
                bars_o[bar_idx].set_color((1.0, 0.5, 0.0))
            else:
                bars_o[bar_idx].set_color((0.92, 0.0, 0.55))

            xx_remp = xx.replace('|', '-')
            img = mpimg.imread('%s/%s.png' % (img_path, xx_remp))
            pp = thumb[xidx].imshow(img)
            pp.set_interpolation('bicubic')
            thumb[xidx].set_yticks([])
            thumb[xidx].set_xticks([])
            thumb[xidx].axis('off')
            #thumb[xidx].set_ylim((50,150))
    else:
        lsv_typ = np.array(range(len(vals)))
        pyplot.barh(lsv_typ, vals)
        pyplot.yticks(lsv_typ, typs)
        pl = pyplot
        offs = 1
        
    if not lim_val is None:
        for lidx, ll in enumerate(lim_val):
            print ll, lims[lidx]
            pl.hlines(ll * offs, 1, lims[lidx]+3000)
            pl.annotate('%s events' % (lims[lidx]), xy=(lims[lidx]+3200, ll * offs))
            
    pyplot.title(extra_title)
    pyplot.xlim((0, 8000))
    pyplot.show()
    fidx += 1


def plot_dominant_exons(dom_dict, name=''):
    global fidx

    fig, ax = pyplot.subplots(1)
    pos = [-1, -0.5, 0, 0.5]

    labels = ['2', '3', '4', '5']
    labels = [xx+' %s' % name for xx in labels]
    totalbins = 0

    for didx, vals in enumerate(dom_dict[1:]):
        nbins = len(set(vals))
        if nbins == 0:
            continue
        totalbins = max(totalbins, nbins)
        bins, attr = np.histogram(vals, nbins)
        x = np.arange(0, nbins)*3 + pos[didx]
        if didx == len(dom_dict[1:])-1:
            labels[didx] += '+'
        lab = labels[didx] + "(%d lsv)" % bins.sum()
        col = [xx/float(255) for xx in BREWER_PALETTE[didx]]
        pyplot.bar(x, bins/float(bins.sum()), width=0.5, label=lab, color=col, edgecolor="none")
#        pyplot.bar(x, bins, width=0.5, label=labels[didx], color=colors[didx])

    xminorlocations = np.arange(0, totalbins, 3) * 3
    # #
    ax.set_xticks(xminorlocations, minor=True)
    # pyplot.grid(True, which='minor', linestyle='-')

    pyplot.xlabel('%s position' % name)
    pyplot.ylabel('Normalized number of LSVs')
    pyplot.title('Dominant exon per LSV. Based on its relative position to the LSV')
    pyplot.legend(loc='best')
    # labs = ["exon %d" % xx for xx in np.arange(0, totalbins+1)]
    # ax.set_xticklabels(labs)
    pyplot.show()
    fidx += 1


def get_types(direc, list_exp, grps):
    d_types = dict()
    g_types = dict()
    for gp in grps:
        g_types[gp] = {}
    for ll in list_exp:
        grp_ll = ""
        for grp in grps:
            if ll.startswith(grp):
                grp_ll = grp
                break
        pp = pickle.load(open("%s/%s" % (direc, ll)))
        kk = []
        for lsv in pp[1]:
            kk.append(lsv.id)
            d_types[lsv.id] = lsv.type
            g_types[grp_ll][lsv.id] = lsv.type

    return d_types, g_types


def all_plots_wrapper(types):

    global fidx
    histo = sorted(types.iteritems(), key=lambda (k, v): (v, k))
    num_ev = np.sum([xx[1] for xx in histo])
    s_keys = [xx[0] for xx in histo if xx[1] > 100]
    s_vals = [types[xx] for xx in s_keys]

    # lim_val = [0]*len(lims)
    # for lidx, l in enumerate(lims):
    #     for vidx, v in enumerate(s_vals):
    #         #print l, vidx, v
    #         if v < l:
    #             lim_val[lidx] = vidx
    #             break
    
    out = open('lsv.short.types', 'w+')
    for tt in s_keys:
        out.write('%s\n' % tt)
    out.close()
    
    impath = './thumbs/'
    extra_title = " %s events" % num_ev
    plot_lsv_types_hist(s_vals, s_keys, img_path=impath, lim_val=None, extra_title=extra_title)

    s_keys = [xx[0] for xx in histo]
    s_vals = [types[xx] for xx in s_keys]

    plot_countings(s_vals, s_keys, [num_ev,0])


def psi_dominant(filename_list):
    res = [[], [], [], [], []]
    for filename in filename_list:
        fp = open(filename)
        lines = fp.readlines()
        fp.close()

        for ll in lines:
            if ll[0] == '#':
                continue
            tab = ll.strip().split('\t')
            lsv_id = tab[1]
            lsv_ex_id = tab[1].split(':')[1]
            strand = tab[12]
            typ = tab[4][0]
            psi_list = [float(xx) for xx in tab[2].split(';')]
            exon_list = [xx for xx in tab[14].split(';')]
            num_exons = len(exon_list) - 1
            if num_exons >= 5:
                num_exons = 5

            exon_psi = {}
            ex_idx = 0
            for ex in exon_list:
                if ex[0] == ' ':
                    ex = ex[1:]
                exshort = int(ex.split('-')[0])

                if ex == lsv_ex_id:
                    continue
                if not exshort in exon_psi:
                    exon_psi[exshort] = 0.0

                exon_psi[exshort] += psi_list[ex_idx]
                ex_idx += 1

            if (strand == '-' and typ == 's') or (strand == '+' and typ == 't'):
                rev = True
            else:
                rev = False

            ex_order = sorted(exon_psi.keys(), reverse=rev)

            max_idx = 0
            max_val = 0

            for ex_idx, ex in enumerate(ex_order):
                if exon_psi[ex] > max_val:
                    max_idx = ex_idx+1
                    max_val = exon_psi[ex]

            if max_val < 0.6:
                max_idx = 0

            res[num_exons-1].append(max_idx)

    return res


def ss_dominant(filename_list):
    res = [[[], [], [], []], [[], [], [], []]]

    for filename in filename_list:
        fp = open(filename)
        lines = fp.readlines()
        fp.close()

        for ll in lines:
            if ll[0] == '#':
                continue
            tab = ll.strip().split('\t')
            lsv_id = tab[1]
            lsv_ex_id = tab[1].split(':')[1]
            strand = tab[12]
            typ = tab[4][0]
            psi_list = [float(xx) for xx in tab[2].split(';')]
            junc_list = [xx for xx in tab[13].split(';')]

            if (strand == '-' and typ == 's') or (strand == '+' and typ == 't'):
                rev = False
                index = 1
            else:
                rev = True
                index = 0

            ss_psi = {}
            for psiidx, jj in enumerate(junc_list):
                if jj[0] == ' ':
                    jj = jj[1:]
                ss = int(jj.split('-')[index])

                if not ss in ss_psi:
                    ss_psi[ss] = 0.0

                ss_psi[ss] += psi_list[psiidx]

            num_ss = len(ss_psi)
            if num_ss >= 4:
                num_ss = 4

            ss_order = sorted(ss_psi.keys(), reverse=rev)
            max_idx = 0
            max_val = 0
            for ss_idx, ss in enumerate(ss_order):
                if ss_psi[ss] > max_val:
                    max_idx = ss_idx+1
                    max_val = ss_psi[ss]

            if max_val < 0.6:
                max_idx = 0

            res[index][num_ss-1].append(max_idx)

    return res


def fdr_parse(fname, group_list):
    fp = open(fname)
    lines = fp.readlines()
    fp.close()

    res = np.zeros(shape=(len(groups), len(groups)*2), dtype=int)

    for ll in lines:
        tab = ll.strip().split()
        fdr = tab[1]
        pair = tab[0].split('_')
        ii = group_list.index(pair[0])
        jj = group_list.index(pair[1])
        res[ii, jj*2] = fdr
        res[ii, jj*2 + 1] = -1 * int(fdr)

    return res


if __name__ == '__main__':

    dire = '/Users/Jordi/notebooks/figure3'
    onlyfiles = [f for f in listdir(dire) if isfile(join(dire, f)) and f.endswith('majiq')]
    groups = ['Heart', 'Hippocampus', 'Liver', 'Lung', 'Spleen', 'Thymus']
    onlyfiles = ['Heart1.mm10.sorted.majiq', 'Hippocampus1.mm10.sorted.majiq']
    groups = ['Heart', 'Hippocampus']
    list_types, group_types = get_types(dire, onlyfiles, groups)

    stypes = {}
    gtypes = dict()
    for grp in groups:
        gtypes[grp] = {}

    for grp in groups:
        for kk, tyt in group_types[grp].items():
            if not tyt in gtypes[grp]:
                gtypes[grp][tyt] = 0
            gtypes[grp][tyt] += 1

    for kk, tyt in list_types.items():
        if not tyt in stypes:
            stypes[tyt] = 0
        stypes[tyt] += 1

    all_plots_wrapper(stypes)

    #read psi values
    groups_vals = []
    filename_list = ['./psi/'+grp+'_psigroup_psi_gene.txt' for grp in groups]
    values = psi_dominant(filename_list)
    plot_dominant_exons(values, 'exons')
    values = ss_dominant(filename_list)
    plot_dominant_exons(values[0], ' 5\'splice sites')
    plot_dominant_exons(values[1], ' 3\'splice sites')

    #heatmap
    fdr_filename = 'allfdr.txt'
    values = fdr_parse(fdr_filename, groups)
    plot_fdrheatmap(values, groups)
