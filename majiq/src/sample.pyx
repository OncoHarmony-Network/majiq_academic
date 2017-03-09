import sys
from random import choice

import numpy as np
import majiq.src.checkpolyfitnb as majiq_fit

"""
Sampling from junctions using a Negative Binomial model.
"""
LIM = 100
EPSILON = 1. / sys.maxint
PSEUDO = 0.0000000001  # EPSILON is too small for some calculations


def _trimborders(junction, border):
    """Discard the borders of the junctions unless they have reads
    :param junction:
    :param border:
    :return:
    """


    # TODO use a masker strategy instead of generating a new array (warning: Assigning values in a numpy
    # array creates a new one, so it is more inneficient than this)
    # discard left side
    new_junction = []
    for i in range(0, border):
        if junction[i] > 0:
            new_junction.append(junction[i])

    new_junction.extend(junction[border:-border])  #add the middle positions, disregard masked positions
    #discard right side
    for i in range(len(junction) - 1, len(junction) - border - 1, -1):
        if junction[i] > 0:
            new_junction.append(junction[i])

    return np.array(new_junction)


def mean_junction(junctions, discardzeros=True):
    """Simple mean of junctions without bootstrapping, but discarding zeroes and flagged stuff
    :param junctions:
    :param discardzeros:
    :return:
    """
    ret = []
    for junc in junctions:
        junc = junc[junc > -EPSILON]
        # mask the -1 (or lower) positions regardless of the discardzero treatment
        if discardzeros:
            junc = junc[junc != 0]
            # a junc array without the zeroes

        if len(junc) == 0:
            ret.append(0)
        else:
            ret.append(junc.mean())

    return np.array(ret)


def sample_from_junctions(junction_list, m, k, discardzeros=5, trimborder=True, fitted_one_over_r=None,
                          debug=False):
    """Given the filtered reads, bootstrap samples from every junction
    :param m:
    :param k:
    :param discardzeros:
    :param trimborder:
    :param fitted_one_over_r:
    :param debug:
    :return:
    :param junction_list:
    """
    sampled_means = []
    sampled_var = []
    all_samples = []

    for i, junction in enumerate(junction_list):
        if 0 < debug == i:
            break
        if i % 100 == 0 and debug > 0:
            print "junction %s..." % i,
            sys.stdout.flush()

        if trimborder:
            junction = _trimborders(junction, trimborder)
            # trim the zeroes from the borders regardless of the discardzeros flag

        junction = junction[junction > -EPSILON]
        # mask the -1 (or lower) positions regardless of the discardzero treatment

        if discardzeros > 0:
            junction = junction[junction != 0]
            # a junction array without the zeroes
            sys.stdout.flush()
            if junction.shape[0] < discardzeros:
                z = np.zeros(shape=(discardzeros - junction.shape[0]), dtype=int)
                junction = np.concatenate((junction, z))
                #a junction array without the zeroes

        if np.count_nonzero(junction) == 0:
            sampled_means.append(0)
            sampled_var.append(0)
            all_samples.append([0] * m)
            # k*m zeroes
        else:

            npos_mult = np.count_nonzero(junction)
            samples = []
            for iternumber in xrange(m):
                junction_samples = []
                for numsamples in xrange(k):
                    junction_samples.append(choice(junction))

                sampled_mean = np.mean(junction_samples)
                if sampled_mean == 0.0:
                    samples.append(0.0)
                    continue

                # recalculating
                nb50 = majiq_fit.sample_over_nb(one_over_r=fitted_one_over_r, mu=sampled_mean, num_samples=k)

                smpl = np.mean(nb50)
                samples.append(smpl)
            # calculate the mean and the variance

            lsv_mean = np.mean(samples)
            var_nb = lsv_mean + fitted_one_over_r * (lsv_mean ** 2)
            sampled_means.append(lsv_mean)
            sampled_var.append(var_nb)

            #            samples = [ npos_mult* (x+1) for x in samples]
            samples = npos_mult * (np.array(samples) + 1)
            #            print i, Nz, samples
            all_samples.append(samples)

    return np.array(sampled_means), np.array(sampled_var), np.array(all_samples)
