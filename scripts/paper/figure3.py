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
import colorbrewer as cb
import os

fidx = 0
lims = [500, 1000]

BREWER_PALETTE = [(228, 26, 28),
                  (55, 126, 184),
                  (77, 175, 74),
                  (152, 78, 163),
                  (255, 127, 0)]


def rgb_to_hex(rgb):
    print rgb
    color = [int(xx) for xx in rgb[4:-1].split(',')]
    return '#%02x%02x%02x' % (color[0], color[1], color[2])


def use_intron(lsvtype):
    res = 'i' in lsvtype
    if ir_plots:
        res = not res

    return res


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
                x += 600
            ax.text(x, y, mm, ha='center', va='bottom', fontsize=size)


def print_message(bars, extrabars=None, filename=None):
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

    if not filename is None:
        with open(filename, 'w+') as fp:
            for mm in msg:
                sms = " ".join(mm)
                fp.write("%s\n" % sms)
    return msg


def plot_fdrheatmap(vals, vals2, grps, output='.'):

    global fidx

    from scipy.cluster.hierarchy import linkage
    from scipy.cluster.hierarchy import dendrogram
    row_clusters = linkage(vals, method='complete')
    row_dendr = dendrogram(row_clusters)


    fig = pyplot.figure(fidx)
    ax = fig.add_subplot(111)

    df_rowclust = np.zeros(shape=vals.shape, dtype=np.int)
    df_rowclust2 = np.zeros(shape=vals.shape, dtype=np.float)


    outfp1 = open('%s/news/fig3f1.vals' % output, 'w+')
    outfp2 = open('%s/news/fig3f2.vals' % output, 'w+')

    header = "Tiss\t"

    group_order = row_dendr['leaves']
    group_order = [10, 5, 9, 11, 0, 1, 2, 8, 7, 3, 6, 4]




    for hh in group_order:
        header += "\t%s" % grps[hh]

    outfp1.write("%s\n" % header)
    outfp2.write("%s\n" % header)
    for xid, xx in enumerate(reversed(group_order)):
        outfp1.write('\t%s' % grps[xx])
        outfp2.write('\t%s' % grps[xx])
        for yid, yy in enumerate(group_order):
            if xx == yy:
                t = 0.0
            else:
                t = float(vals2[xx, yy])/vals[xx, yy]
            df_rowclust[xid, yid] = vals[xx, yy]
            df_rowclust2[xid, yid] = t
            outfp1.write('\t%s' % df_rowclust[xid, yid])
            outfp2.write('\t%.3f' % df_rowclust2[xid, yid])
        outfp1.write('\n')
        outfp2.write('\n')

    outfp1.close()
    outfp2.close()

    cax = ax.matshow(df_rowclust, interpolation='nearest', cmap='Purples')

    for x in range(df_rowclust.shape[0]):
        for y in range(df_rowclust.shape[1]):
            if x == vals.shape[0]-y - 1:
                continue
            pyplot.text(y, x, '%d' % df_rowclust[x, y],
                        horizontalalignment='center',
                        verticalalignment='center',
                        )

    rev_df_grps = [grps[xx] for xx in group_order]
    df_grps = [grps[xx] for xx in reversed(group_order)]
    ax.set_xlim(-0.5, 12)
    ax.set_ylim(-0.5, 12)
    tks = ax.get_xticks()


    ticks = np.arange(0, 12)
    ax.set_xticks(ticks)
    ax.set_xticklabels(rev_df_grps)
    tks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels(df_grps)
    ax.set_title('Changing events')
    #
    fig.colorbar(cax)
    pyplot.show()
    fidx += 1


    fig = pyplot.figure(fidx)
    ax = fig.add_subplot(111)

    for x in range(vals2.shape[0]):
        for y in range(vals2.shape[1]):
            if x == vals2.shape[0] - y - 1:
                continue
            pyplot.text(y, x, '%.2f' % df_rowclust2[x, y],
                        horizontalalignment='center',
                        verticalalignment='center',
                        )

    cax = ax.matshow(df_rowclust2, interpolation='nearest', cmap='Greens')
    fig.colorbar(cax)
    ax.set_xlim(-0.5, 12)
    ax.set_ylim(-0.5, 12)
    ax.set_xticks(ticks)
    ax.set_xticklabels(rev_df_grps)
    tks = ax.get_yticks()
    ax.set_yticks(ticks)
    ax.set_yticklabels(df_grps)
    ax.set_title('Complex changing events')

    # if output == '.':
    pyplot.show()
    # else:
    #     pyplot.savefig(output+'/hog.pdf')


    # plot

    # ax.set_xticklabels(grps)
    # ax.set_yticklabels(grps)
    # tks = ax.get_xticks()
    # ax.set_xticks(tks+1)
    # tks = ax.get_yticks()
    # ax.set_yticks(tks+0.5)
    # pyplot.show()
    fidx += 1


def plot_countings(vals, typs, counts, output='.'):

    global fidx

    offs = [0, 0.8, 1.60, 2.40]

    pyplot.figure(fidx)
    junc_bars = [0] * 20
    ex_bars = [0] * 20
    bars_s = [0] * 20
    bars_t = [0] * 20

    for idx, tt in enumerate(typs):
        if use_intron(tt):
            continue
        juncs = [xx for xx in tt.split('|')[1:] if not xx.endswith('0')]
        nums = len(juncs)
        if nums >= 20:
            nums = 19
        junc_bars[nums] += vals[idx]

        exnum = []
        ssnum = []
        for jj in juncs:

            exnum.append(int(jj.split('e')[1].split('.')[0]))
            ssnum.append(int(jj.split('e')[0]))

        nums_ex = max(exnum)
        if nums_ex >= 20:
            nums_ex = 19
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
    msg = print_message(junc_bars, filename='%s/news/fig3b.junctions.vals' % output)

    lim = max(dummy_ex)
    ex_bars = np.array(ex_bars, dtype=float)
    nbars = ex_bars / ex_bars.sum()
    col = [xx/float(255) for xx in BREWER_PALETTE[1]]
    pyplot.bar(x[:lim] + offs[1], nbars[1:lim+1], edgecolor="none", label="#LSV per #exons", width=0.8,
               color=col, alpha=0.5)

    msg = print_message(ex_bars, filename='%s/news/fig3b.exons.vals' % output)


    lim = max(max(dummys), max(dummyt))
    bars_s = np.array(bars_s, dtype=float)
    nbars_s = bars_s/bars_s.sum()
    bars_t = np.array(bars_t, dtype=float)
    nbars_t = bars_t/bars_t.sum()
    pyplot.bar(x[:lim] + offs[2], nbars_s[1:lim+1], label="5' splicesites in Source LSV", edgecolor="none",
               color=(1.0, 0.5, 0.0), width=0.8, alpha=0.5)
    pyplot.bar(x[:lim] + offs[3], nbars_t[1:lim+1], label="3' splicesites in Target LSV", edgecolor="none",
               color=(0.92, 0.0, 0.55), width=0.8, alpha=0.5)

    msg = print_message(bars_s, filename='%s/news/fig3b.ss5prime.vals' % output)
    msg = print_message(bars_t, filename='%s/news/fig3b.ss3prime.vals' % output)

    pyplot.ylabel(' # LSVs')
    pyplot.xlabel(' # elements')
    pyplot.title('Total Number of LSVs %s' % counts)
    pyplot.legend(loc='best')
    if output == '.':
        pyplot.show()
    else:
        pyplot.savefig(output+'.pdf')
    fidx += 1


def plot_lsv_types_hist(vals, typs, img_path=None, lim_val=None, extra_title="", output='.'):

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
              
        bars_o = pl.barh(lsv_typ, vals2, height=4, edgecolor="none", alpha=0.5)
        
        pl.set_ylim((0, len(vals2)))
        pl.set_yticks([])

        msg = print_message(np.array(vals2, dtype=np.float), filename='%s/news/fig3a.vals' % output)



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
    #pyplot.xlim((0, 8000))
    if output == '.':
        pyplot.show()
    else:
        pyplot.show()
        pyplot.savefig(output+'types.pdf')
    fidx += 1


def plot_dominant_exons(dom_dict, name='', color=cb.Blues[9], output='.'):
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
        pyplot.bar(x, bins/float(bins.sum()), width=0.5, label=lab, color=rgb_to_hex(color[-1*(didx*2+1)]), edgecolor="none")
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
    if output == '.':
        pyplot.savefig('%s.hist.pdf' % name)
        pyplot.show()

    else:
        pyplot.show()
        pyplot.savefig(output+'hist.pdf')
    fidx += 1


def collapse(type_str):

    tab = type_str.split('|')
    res = tab[0]

    dd = {}

    for jj in tab[1:]:
        if jj == 'i':
            #dd[dest] = 0
            continue
        dest = jj.split('e')
        tab2 = dest[1].split('.')
        if dest[1] == '0':
            continue
        if not tab2[0] in dd:
            dd[tab2[0]] = []

        dd[tab2[0]].append(tab2[1].split('o')[0])

    translate = dict()
    translate['i'] = 'i'
    for kk, vv in dd.items():
        unique_vv = list(set(vv))
        unique_vv.sort()

        for vidx, v in enumerate(unique_vv):
            translate["%s.%s" % (kk, v)] = "%s.%d" % (kk, vidx+1)

    for jj in tab[1:]:
        dest = jj.split('e')
        if dest[0] != 'i' and dest[1] == '0':
            res += "|%s" % jj
            continue
        if dest[0] == 'i':
            res += '|i'
        else:
            res += '|%se%s' % (dest[0], translate[dest[1].split('o')[0]])

    return res


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
            if use_intron(lsv.type):
                continue
            kk.append(lsv.id)
            d_types[lsv.id] = collapse(lsv.type)
            g_types[grp_ll][lsv.id] = collapse(lsv.type)

    return d_types, g_types


def all_plots_wrapper(types, nlsv=0, output='.'):

    global fidx
    histo = sorted(types.iteritems(), key=lambda (k, v): (v, k))
    num_ev = np.sum([xx[1] for xx in histo])
    #s_keys = [xx[0] for xx in histo if xx[1] > 100]
    s_keys = [xx[0] for xx in histo if xx[1] > 100]
    s_vals = [types[xx] for xx in s_keys]

    all_types = [xx[0] for xx in histo]

    total = 0
    complex = 0
    for xx in all_types:
        if len(xx.split('|')[1:]) > 2:
            complex += types[xx]
        total += types[xx]

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
    
    impath = '%s/thumbs/' % output
    extra_title = " %s events" % num_ev
    percent = float(complex) / total
    print "Complex %d/%d (%.3f)" % (complex, total, percent)

    #output = '.'
    plot_lsv_types_hist(s_vals, s_keys, img_path=impath, lim_val=None, extra_title=extra_title, output=output)
    #plot_lsv_types_hist(s_vals, s_keys, img_path=None, lim_val=None, extra_title=extra_title)


    s_keys = [xx[0] for xx in histo]
    s_vals = [types[xx] for xx in s_keys]

    plot_countings(s_vals, s_keys, num_ev, output=output)


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
            lsv_id = tab[2]
            lsv_ex_id = tab[2].split(':')[1]
            strand = tab[13]
            typ = tab[5][0]

            if use_intron(tab[5]):
                continue

            psi_list = [float(xx) for xx in tab[3].split(';')]
            exon_list = [xx for xx in tab[15].split(';')]
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


def ss_dominant(file_list):
    res = [[[], [], [], []], [[], [], [], []]]

    for filename in file_list:
        fp = open(filename)
        lines = fp.readlines()
        fp.close()

        for ll in lines:
            if ll[0] == '#':
                continue
            tab = ll.strip().split('\t')
            lsv_id = tab[2]
            lsv_ex_id = tab[2].split(':')[1]
            strand = tab[13]
            typ = tab[5][0]
            if use_intron(tab[5]):
                continue
            psi_list = [float(xx) for xx in tab[3].split(';')]
            junc_list = [xx for xx in tab[14].split(';')]

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


def fdr_parse(direc, file_list, group_list, intronic_junc=False):

    changing = np.zeros(shape=(len(group_list), len(group_list)), dtype=np.int)
    complx = np.zeros(shape=(len(group_list), len(group_list)), dtype=np.int)

    total_dict = set()
    total_dict_complx = set()

    per_tiss = {}
    per_tiss_complx = {}

    for xx in group_list:
        per_tiss[xx] = set()
        per_tiss_complx[xx] = set()

    for filename in file_list:
        fp = open("%s/%s" % (direc, filename))
        lines = fp.readlines()
        fp.close()

        pair = filename.split('.')[0].split('_')
        x = group_list.index(pair[0])
        y = group_list.index(pair[1])

        for ll in lines:
            if ll[0] == '#':
                continue
            tab = ll.strip().split('\t')
            if use_intron(tab[7]):
                continue
            typ = tab[7].split('|')[1:]

            ntyp = 0
            #print tab[7], typ
            for tt in typ:
                if tt.endswith('e0'):
                    continue
                ntyp += 1
            psi_list = [float(xx) for xx in tab[3].split(';')]

            if intronic_junc:
                psi_list_check = psi_list[-1:]
            else:
                # if abs(psi_list[-1]) < 0.2:
                #     psi_list_check = psi_list[:-1]
                # else:
                #     continue
                psi_list_check = psi_list

            for pp in psi_list_check:
                if abs(pp) > 0.2:
                    changing[x, y] += 1
                    changing[y, x] += 1
                    total_dict.add(tab[2])
                    per_tiss[pair[0]].add(tab[2])
                    per_tiss[pair[1]].add(tab[2])
                    if ntyp > 2:
                        total_dict_complx.add(tab[2])
                        per_tiss_complx[pair[0]].add(tab[2])
                        per_tiss_complx[pair[1]].add(tab[2])
                        complx[x, y] += 1
                        complx[y, x] += 1
                    break

    print "TOTAL %s" % len(total_dict)
    print "COMPLX %s" % len(total_dict_complx)

    fract = float(len(total_dict_complx))/len(total_dict)

    print "FRACTION", fract

    print "TOTAL PER TISSUE"
    for kk, vv in per_tiss.items():
        print kk, len(vv)

    for kk, vv in per_tiss_complx.items():
        print kk, len(vv)

    return changing, complx



if __name__ == '__main__':

    dire = sys.argv[1]
    onlyfiles = [f for f in listdir(dire) if isfile(join(dire, f)) and f.endswith('majiq')]
    groups = ['Adr', 'Aor', 'BFat', 'Bstm', 'Cer', 'Hrt', 'Hyp', 'Kid', 'Liv', 'Lun', 'Mus', 'WFat']
    # onlyfiles = ['Adr_CT22.mm10.sorted.majiq', 'Aor_CT22.mm10.sorted.majiq']

    global ir_plots
    ir_plots = False

    output = sys.argv[2]
    #groups = sys.argv[3:]

    print "Parse files"
    if not os.path.exists(output):
        os.makedirs(output)
        os.makedirs('%s/news' % output)

    # list_types, group_types = get_types(dire, onlyfiles, groups)
    # count_lsv = len(set(list_types.keys()))
    # stypes = {}
    # for kk, tyt in list_types.items():
    #     if not tyt in stypes:
    #         stypes[tyt] = 0
    #     stypes[tyt] += 1
    # print "Plot 3.a 3.b"
    # all_plots_wrapper(stypes, count_lsv, output=output)


    #
    #read psi values
    # groups_vals = []
    # filename_list = ['./psi/%s_psigroup_psi.txt' % grp for grp in groups]
    # print "Plot Dominant exons"
    # values = psi_dominant(filename_list)
    # plot_dominant_exons(values, 'exons', color=cb.Blues[9])
    # values = ss_dominant(filename_list)
    # print "Plot Dominant splicesites"
    # plot_dominant_exons(values[0], ' 5\'splice sites', color=cb.Oranges[9])
    # plot_dominant_exons(values[1], ' 3\'splice sites', color=cb.PuRd[9])

    #heatmap
    print "Plot dpsi changing events heatmap"
    dire = './dpsi_voila'
    filename_list = [f for f in listdir(dire) if isfile(join(dire, f)) and f.endswith('txt')]
    chg_lsv, complx_lsv = fdr_parse(dire, filename_list, groups, intronic_junc=ir_plots)
    plot_fdrheatmap(chg_lsv, complx_lsv, groups, output)


