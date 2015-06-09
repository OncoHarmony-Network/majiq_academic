from matplotlib import use

use('Agg', warn=False)
from pylab import *
import numpy as np
from scipy.stats import pearsonr
import src.sample
import src.polyfitnb as polyfitnb
import os
import cPickle as pickle

DEBUG = True
TESTBREAK = 1500
LIM = 100
EPSILON = 1. / sys.maxint
BORDER = 5  # definition of what a border is


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
        savefig("%s/%s.png" % (plotpath, plotname.replace(" ", "_")), bbox_inches='tight', width=1000, height=2000,
                dpi=300)  # No spaces allowed, underscores!
        clf()
    else:
        show()


def calc_nonzeromeanvar(junctions):
    nonzeromean = []
    nonzerovar = []
    for junction in junctions:
        nonzerojunction = junction[junction <= 0]  # discard also -1 and -2, just in case
        if len(nonzerojunction) > 0:
            nonzeromean.append(mean(nonzerojunction))
            nonzerovar.append(var(nonzerojunction))
        else:
            nonzeromean.append(0)
            nonzerovar.append(0)

    return array(nonzeromean), array(nonzerovar)


def calc_weights(junction):
    # Not used yet. The weight_factors array is wrong
    weight_factors = array([0.0, 0.67221794713393235, 0.822110928783363, 1.0, 0.96555466736662077, 0.97034079788351413,
                            0.95888448040949026, 0.93916171702266071,
                            0.9586120159245054, 0.98886869384130682, 0.97031850460351132, 0.97919730108521363,
                            0.98845338964933505, 0.97762016912650862,
                            0.96000443463677787, 0.97795543996841916, 0.98483897464546255, 0.98430701054559211,
                            0.97621404631851538, 0.97557091482162561,
                            0.99783624218670419, 0.99420256804654417, 0.99996092005852, 1.0, 0.98891003681022327,
                            0.98408910298925234, 0.98588911669260937,
                            0.98944197552348001, 0.98861266787997559, 0.98334059809099128, 0.98616818121835459,
                            0.98356568445706294, 1.0, 1.0,
                            0.99415876734414588, 1.0, 0.99732413178991319, 0.98571657557568526, 0.98637294512249951,
                            0.98846297241187242, 1.0,
                            0.98857076368303576, 1.0, 0.98474007029306976, 0.98212050612598556, 0.99227062085183826,
                            0.98716235724225032, 0.98604617629365343,
                            0.96908030440229109, 0.97105918006649872, 0.97297718733803484, 0.98431591864639367,
                            0.98227616224387038, 0.9961571944449884,
                            0.97565056267585271, 0.96725772937340826, 0.95469906291036666, 0.94761567083759468,
                            0.96284719373281014, 1.0, 0])

    # TODO correction for the borders
    return weight_factors / len(weight_factors)


def _trimborders(junction):
    "Discard the borders of the junctions unless they have reads"
    # discard left side
    new_junction = []
    for i in range(0, BORDER):
        if junction[i] > 0:
            new_junction.append(junction[i])

    new_junction.extend(junction[BORDER:-BORDER])  # add the middle positions
    # discard right side
    for i in range(len(junction) - 1, len(junction) - BORDER - 1, -1):
        if junction[i] > 0:
            new_junction.append(junction[i])

    return array(new_junction)


def sample_from_junctions(junctions, m, k, discardzeros=False, nb=False, trimborder=False, fit_func=None,
                          poisson=False):
    if nb:
        return src.sample.sample_from_junctions(junctions, m, k, discardzeros=discardzeros, trimborder=trimborder,
                                                fitted_one_over_r=fit_func)

    # if parameters:
    # print "Loading parameters from %s..." % parameters
    # fitted_func = pickle.load(open('%sfitfunc.pickle' % parameters))
    #     a, b = fitted_func.c

    sampled_means = []
    sampled_var = []
    all_samples = []
    for i, junction in enumerate(junctions):
        # if i % 100 == 0:
        #     print "junction %s..."%i,
        #     sys.stdout.flush()

        junction = junction[junction > -EPSILON]  #discard the -1 (or lower) positions regardless of the dzero treatment

        if trimborder:
            junction = src.sample._trimborders(junction,
                                               trimborder)  #trim the zeroes from the borders regardless of the discardzeros flag
        if discardzeros:
            junction = junction[junction != 0]  #a junction array without the zeroes

        if poisson:
            mean_poisson = np.mean(junction, axis=0)
            sampled_means.append(mean_poisson)
            sampled_var.append(mean_poisson)
            all_samples.append(junction)
        else:
            import random

            if len(junction) == 0:
                sampled_means.append(0)
                sampled_var.append(0)
                all_samples.append([0] * (k * m))  #k*m zeroes
            else:
                samples = []
                for iternumber in xrange(m):
                    junction_samples = []
                    #using the matrix
                    #weights = calc_weights(junction)
                    #junction_samples = multinomial(k, weights, junction) #This function is very slow
                    for numsamples in xrange(k):
                        junction_samples.append(random.choice(junction))

                    samples.extend(junction_samples)

                #calculate the mean and the variance
                sampled_means.append(mean(samples))
                sampled_var.append(var(samples))
                all_samples.append(samples)

    return array(sampled_means), array(sampled_var), array(all_samples)


def plot_pearsoncorr(var1, var2, my_title, my_xlabel, my_ylabel, plotpath=None, max_value=None):
    if DEBUG:
        var1 = var1[:TESTBREAK]
        var2 = var2[:TESTBREAK]

    var1 = array(var1)
    var2 = array(var2)
    xlabel(my_xlabel)
    ylabel(my_ylabel)

    if not max_value:
        max_value = max(max(var1), max(var2))

    xlim(0, max_value)
    ylim(0, max_value)

    # plot([0, max_value], [0, max_value])
    print "PEARSON", var1, var2
    pear, pvalue = pearsonr(var1, var2)
    r_squared = pear ** 2

    a, b = polyfit(var1, var2, 1)
    fit_func = poly1d([a, b])
    plot(var1, fit_func(var1), '-r')
    # percentage_under = sum(var1 < var2)/float(len(var1))
    text(abs(max_value) * 0.1, max_value - abs(max_value) * 0.2, r'$R^2$: %.2f (p-value: %.2E)' % (r_squared, pvalue),
         fontsize=18, bbox={'facecolor': 'yellow', 'alpha': 0.3, 'pad': 10})
    title(my_title)
    print r"%s R^2: %.2f (p-value: %.2E)" % (my_title, r_squared, pvalue)
    plot(var1, var2, '.')
    if plotpath:
        _save_or_show(plotpath, my_title)


def filter_bulk(matrix_filter, *matrices):
    ret = []
    for m in matrices:
        ret.append(m[matrix_filter])

    return ret


def check_junctions_in_replicates(lsv_junc1, lsv_junc2, discard_empty_junctions=False):
    ids1 = set([x[1] for x in lsv_junc1[1]])
    ids2 = set([x[1] for x in lsv_junc2[1]])

    matched_names = ids1.intersection(ids2)
    print len(ids1), len(ids2)
    print len(matched_names)
    gc1 = []
    gc2 = []
    replica1 = []
    replica2 = []
    for ii in matched_names:
        for idx, nm in enumerate(lsv_junc1[1]):
            if nm[1] == ii:
                replica1.append(lsv_junc1[0][idx])
                gc1.append(lsv_junc1[2])
                dummy = lsv_junc1[0][idx].shape
                break
        for idx, nm in enumerate(lsv_junc2[1]):
            if nm[1] == ii:
                dummy2 = lsv_junc2[0][idx].shape
                if dummy != dummy2:
                    print "ERRRRORRRRRR", dummy, dummy2
                    replica1 = replica1[:-1]
                    gc1 = gc1[:-1]
                else:
                    replica2.append(lsv_junc2[0][idx])
                    gc2.append(lsv_junc2[2])

                break

    replica1 = np.concatenate(replica1)
    replica1 = replica1.astype(np.float64)
    replica2 = np.concatenate(replica2)
    replica2 = replica2.astype(np.float64)

    gc1 = np.concatenate(gc1)
    gc1 = gc1.astype(np.float64)
    gc2 = np.concatenate(gc2)
    gc2 = gc2.astype(np.float64)

    if discard_empty_junctions:
        idx_list = []
        for idx in range(replica1.shape[0]):
            if np.count_nonzero(replica1[idx]) == 0 or np.count_nonzero(replica2[idx]) == 0: idx_list.append(idx)
        replica1 = np.delete(replica1, idx_list, axis=0)
        replica2 = np.delete(replica2, idx_list, axis=0)
        gc1 = np.delete(gc1, idx_list, axis=0)
        gc2 = np.delete(gc2, idx_list, axis=0)

    return replica1, replica2, gc1, gc2


def discard_empty_junctions(junc_list1, junc_list2):
    ids1 = set(junc_list1[1])
    ids2 = set(junc_list2[1])

    matched_names = ids1.intersection(ids2)

    replica1 = []
    replica2 = []
    for ii in matched_names:
        for idx, nm in enumerate(junc_list1[1]):
            if nm == ii:
                replica1.append(junc_list1[0][idx])
                break
        for idx, nm in enumerate(junc_list2[1]):
            if nm == ii:
                replica2.append(junc_list2[0][idx])
                break

    # print replica1
    replica1 = array(replica1)
    replica1 = replica1.astype(np.float64)
    replica2 = array(replica2)
    replica2 = replica2.astype(np.float64)

    if discard_empty_junctions:
        idx_list = []
        for idx in range(replica1.shape[0]):
            if np.count_nonzero(replica1[idx]) == 0 or np.count_nonzero(replica2[idx]) == 0: idx_list.append(idx)
        replica1 = np.delete(replica1, idx_list, axis=0)
        replica2 = np.delete(replica2, idx_list, axis=0)

    return replica1, replica2


def load_data_lsv(path, group_name, logger=None):
    """Load data from the preprocess step. Could change to a DDBB someday"""
    data = pickle.load(open(path))
    lsv_cov_list = []
    lsv_info = []
    const_info = []
    num_pos = data[1][0].junction_list.shape[1]
    # num_pos = data[0][0].junction_list.shape[1]

    for lsv in data[1]:
        # for lsv in data[0]:
        try:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0, lsv.visual])
        except AttributeError, e:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0])

        cov = lsv.junction_list.toarray()
        lsv_cov_list.append(cov)

    import random

    clist = random.sample(data[2], min(5000, len(data[2])))
    # clist = random.sample(data[1], min(5000, len(data[1])))
    const_list = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('int'))

    for cidx, const in enumerate(clist):
        const_info.append(const.id)
        const_list[cidx, :] = const.coverage.toarray()

    return (lsv_cov_list, lsv_info), (const_list, const_info)


def load_junctions(filename1, filename2, args, fromlsv=False):
    # Parse LSV files
    # meta1, lsv_junc1, const1 = majiqio.load_data_lsv(filename1, os.path.basename(filename1).split('.')[0])
    # meta2, lsv_junc2, const2 = majiqio.load_data_lsv(filename2, os.path.basename(filename2).split('.')[0])
    lsv_junc1, const1 = load_data_lsv(filename1, os.path.basename(filename1).split('.')[0])
    lsv_junc2, const2 = load_data_lsv(filename2, os.path.basename(filename2).split('.')[0])

    # print const1[0].shape
    # print const2[0].shape

    fit_func1 = polyfitnb.fit_nb(const1[0], "%s_nbfit" % args.output, args.plotpath, nbdisp=args.dispersion,
                                 logger=None)
    fit_func2 = polyfitnb.fit_nb(const2[0], "%s_nbfit" % args.output, args.plotpath, nbdisp=args.dispersion,
                                 logger=None)

    return lsv_junc1, lsv_junc2, fit_func1, fit_func2, [], [], const1, const2


def split_junction_pool(replica1, replica2):
    rep1 = [[], [], []]
    rep2 = [[], [], []]

    for idx in range(replica1.shape[0]):
        if np.count_nonzero(replica1[idx]) in range(1, 6):
            rep1[0].append(replica1[idx])
            rep2[0].append(replica2[idx])
        elif np.count_nonzero(replica1[idx]) in range(6, 16):
            rep1[1].append(replica1[idx])
            rep2[1].append(replica2[idx])
        elif np.count_nonzero(replica1[idx]) > 15:
            rep1[2].append(replica1[idx])
            rep2[2].append(replica2[idx])

    for idx in range(3):
        rep1[idx] = array(rep1[idx])
        rep2[idx] = array(rep2[idx])
    return rep1, rep2


def main():
    print "Deprecated"
    import sys

    sys.exit(1)


if __name__ == '__main__':
    main()
