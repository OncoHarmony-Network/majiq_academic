# This file contains old functions that are not useanymore, but we would prefer not remove for now.

def plot_exons(vals, typs, extra_title="", subplt=None):
    global fidx
    bars = [0] * 20
    if subplt is None:
        fig = pyplot.figure(fidx)
        plt = pyplot
    else:
        fig = subplt
        plt = subplt

    for idx, tt in enumerate(typs):
        juncs = [xx for xx in tt.split('|')[1:] if not xx.endswith('0')]
        exnum = []
        for jj in juncs:
            exnum.append(int(jj.split('e')[1].split('.')[0]))
        nums = max(exnum)
        bars[nums] += vals[idx]

    dummy = [idx for idx, xx in enumerate(bars) if xx > 0]
    lim = max(dummy)
    bars = np.array(bars, dtype=float)
    nbars = bars / bars.sum()
    bar_o = plt.bar(np.arange(2, lim+1), bars[2:lim+1])
    msg = print_message(bars[2:lim+1])
    autolabel(bar_o, plt, msg, size=8)

    if subplt is None:
        pyplot.ylabel(' #lsv')
        pyplot.xlabel(' #exons')
        pyplot.title('Number lsv per #exons')
        pyplot.show()
        fidx += 1
    else:
        plt.set_ylabel(' #lsv')
        plt.set_xlabel(' #exons')
        plt.set_title('Number lsv per #exons')


def plot_splicesites(vals, typs, extra_title="", subplt=None):
    global fidx
    if subplt is None:
        fig = pyplot.figure(fidx)
        plt = pyplot
    else:
        fig = subplt
        plt = subplt

    bars_S = [0] * 20
    bars_T = [0] * 20
    for idx, tt in enumerate(typs):
        if tt[0] == 's':
            bars = bars_S
        else:
            bars = bars_T
        juncs = [xx for xx in tt.split('|')[1:] if not xx.endswith('0')]
        exnum = []
        for jj in juncs:
            exnum.append(int(jj.split('e')[0]))
        nums = max(exnum)
        bars[nums] += vals[idx]

    dummys = [idx for idx, xx in enumerate(bars_S) if xx > 0]
    dummyt = [idx for idx, xx in enumerate(bars_T) if xx > 0]



    lim = max(max(dummys), max(dummyt))
    bars_S = np.array(bars_S, dtype=float)
    bar_oS = plt.bar(np.arange(1, 2*lim + 1, 2)-0.6, bars_S[1:lim+1], label="5' splicesites in Source LSV",
                     edgecolor="none", color=(1.0, 0.5, 0.0), width=0.8)

    bars_T = np.array(bars_T, dtype=float)
    bar_oT = plt.bar(np.arange(1, 2*lim + 1, 2)+0.2, bars_T[1:lim+1], label="3' splicesites in Target LSV",
                     edgecolor="none", color=(0.92, 0.0, 0.55), width=0.8)

    msg = print_message(bars_S[1:lim+1], extrabars=bars_T[1:lim+1])
    autolabel(bar_oS, plt, msg, size=8)
    #
    # msg = print_message(bars_T[1:lim+1])
    # autolabel(bar_oT, plt, msg, size=6, b=True)


    plt.set_xticklabels(np.arange(1, lim+1))
    ylab = '#lsv'
    xlab = '#splicesites'
    tit = 'Number lsv per #splicesites'

    if subplt is None:
        pyplot.ylabel(ylab)
        pyplot.xlabel(xlab)
        pyplot.title(tit)
        pyplot.show()
        fidx += 1
    else:
        plt.set_ylabel(ylab)
        plt.set_xlabel(xlab)
        plt.set_title(tit)


def plot_juncs(vals, typs, extra_title="", subplt=None):
    global fidx
    bars = [0] * 20
    if subplt is None:
        fig = pyplot.figure(fidx)
    else:
        fig = subplt
        plt = subplt
    for idx, tt in enumerate(typs):
        juncs = [xx for xx in tt.split('|')[1:] if not xx.endswith('0')]
        nums = len(juncs)
        bars[nums] += vals[idx]

    dummy = [idx for idx, xx in enumerate(bars) if xx > 0]
    lim = max(dummy)
    bars = np.array(bars, dtype=float)
    nbars = bars / bars.sum()
    bar_o = plt.bar(np.arange(2, lim+1), bars[2:lim+1])
    msg = print_message(bars[2:lim+1])
    autolabel(bar_o, plt, msg, size=8)

    if subplt is None:
        pyplot.ylabel(' #lsv')
        pyplot.xlabel(' #junctions')
        pyplot.title('Number lsv per #junctions')
        pyplot.show()
        fidx += 1
    else:
        plt.set_ylabel(' #lsv')
        plt.set_xlabel(' #junctions')
        plt.set_title('Number lsv per #junctions')