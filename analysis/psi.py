import argparse
from collections import defaultdict
from scipy.special import gamma, gammaln
#deprecation WARNING comes from this import!!!
from scipy.stats import binom_test
from numpy.random import dirichlet
import numpy as np
import sys
import cPickle as pickle
import matplotlib.pyplot as plt
import analysis.filter as majiq_filter
import analysis.adjustdelta as majiq_delta
import analysis.sample as majiq_sample
import operator
import os

"""
Calculate and manipulate PSI and Delta PSI values
"""
BSIZE = 0.025 #TODO To parameters
BINS = np.arange(0, 1, BSIZE)
# The bins for PSI values. With a BSIZE of 0.025, we have 40 BINS
BINS_CENTER = np.arange(0+BSIZE/2, 1, BSIZE)
#The center of the previous BINS. This is used to calculate the mean value of each bin.


def plot_matrix(matrix, my_title, plotname, plotpath):
    plt.clf()
    ax = plt.subplot(1, 1, 1)
    plt.title(my_title)
    #plt.imshow(matrix)
    plt.xlabel(u"PSI i")
    plt.ylabel(u"PSI j")
    ax.set_xticklabels([0, 0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0, 0.25, 0.5, 0.75, 1])

    _save_or_show(plotpath, plotname=plotname)


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        plt.savefig("%s%s.png" % (plotpath, plotname), bbox_inches='tight')
        plt.clf()
    else:
        plt.show()


def median_psi(junctions, discardzeros=True):
    "Calculates the median PSI for all events"
    medians = []
    for junction in junctions:
        if discardzeros:
            junction = junction[junction != 0]
            #a junction array without the zeroes

        medians.append(np.median(junction))

    return np.array(medians)


def empirical_delta_psi(lsv_list1, lsv_list2, logger=None):
    """Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions"""

#    if  not logger: logger.info("Calculating PSI for 'best set'...")

    delta_psi = []
    delta_psi_ir = []
    for idx, lsv in enumerate(lsv_list1):
        psi1 = np.zeros(shape=len(lsv), dtype=np.dtype('float'))
        psi2 = np.zeros(shape=len(lsv), dtype=np.dtype('float'))
        for ii, rate in enumerate(lsv):
            val = float(rate) / float(np.sum(lsv))
            if np.isnan(val):
                val = 0.5
            psi1[ii] = val

        for ii, rate in enumerate(lsv_list2[idx]):
            val = float(rate) / float(np.sum(lsv_list2[idx]))
            if np.isnan(val):
                val = 0.5
            psi2[ii] = val

        sys.stdout.flush()

        delta_psi.append(psi1 - psi2)
 #   if logger: logger.info("Calculating delta PSI for 'best set'...")
    return delta_psi 


def simple_psi(inc, exc):
    """Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions"""
    psi = inc/(exc+inc)
    psi[np.isnan(psi)] = 0.5
    #if NaN, is because exc+inc = 0. If we know nothing, then we don't know if its 0 (exclusion) or 1 (inclusion)
    return psi 


def reads_given_psi_lsv(lsv, psi_space):
    #P(vector_i | PSI_i)
    "We do a simple binomial test to evaluate how probable is the data given a PSI range"
    psi = np.zeros(shape=(lsv.shape[0], psi_space.sh8ape[0]), dtype=np.float)
    for idx, junc in enumerate(lsv):
        total_psi = np.zeros(shape=(100, psi_space.shape[0]), dtype=np.float)
        for pidx, smpl in enumerate(junc):
            bin_test = [binom_test(smpl, lsv[:, pidx].sum(), p=x) for x in psi_space]
            bin_test = np.array(bin_test) + 1e-10
            total_psi[pidx] = (bin_test/bin_test.sum())
        total_psi = np.mean(total_psi, axis=0)
        psi[idx] = total_psi/total_psi.sum()

    return psi


def reads_given_psi(inc_samples, exc_samples, psi_space):
    #P(vector_i | PSI_i)
    "We do a simple binomial test to evaluate how probable is the data given a PSI range"
    ret = []
    inc = inc_samples.sum(axis=1)
    exc = exc_samples.sum(axis=1)
    for i in xrange(inc.shape[0]):
        event = []
        for psi_val in psi_space:
            event.append(binom_test(inc[i], exc[i]+inc[i], p=psi_val))

        ret.append(np.array(event) / sum(event))
    
    return np.array(ret).reshape(-1, len(psi_space))


class DirichletCalc:
    def __init__(self):
        self.cache = defaultdict(float)

    def pdf(self, x, alpha):
        k = x
        k.extend(alpha)
        key = " ".join(map(str, k))
        if self.cache[key]:
        #try to figure out if we already calculated this pdf
            return self.cache[key] 
        else: 
            #formula taken from stackoverflow "How to calculate dirichlet PDF", author based it on
            #  Wikipedia's PDF definition
            ret = gamma(sum(alpha)) / reduce(operator.mul, [gamma(a) for a in alpha]) \
                * reduce(operator.mul, [x[i] ** (alpha[i]-1.0) for i in xrange(len(alpha))])
            self.cache[key] = ret
            return ret


def dirichlet_pdf(x, alpha):
    """Returns a Dirichlet PDF function"""
    alphap = alpha - 1
    c = np.exp(gammaln(alpha.sum()) - gammaln(alpha).sum())
    dir_res = c * (x**alphap).prod(axis=1)
    dir_res = np.array(dir_res)
    psi = dir_res / float(dir_res.sum())

    return psi


def recalibrate_delta(deltapsi):
    #TODO make deltaPSI follow the following binning system
    np.arange(-98.75, 100, 2.5)


def lsv_psi(samples_events, alpha, n, debug):
    """
    Given a set of matching inclusion and exclusion samples, calculate psi, save it in disk, and
    return the psi-per-juntion matrix
    """

    psi_scores = []
    dircalc = DirichletCalc() 
    for i, lsv in enumerate(samples_events):

        if i % 50 == 0:
            print "event %s..." % i,
            sys.stdout.flush()
        if 0 < debug == i: break
        psi = np.zeros(shape=(lsv.shape[0],BINS.shape[0]), dtype=np.float)
        #if debug: print "Paired samples to dirichlet..."
        #sampling PSI by pairing the samples of the previous step sequentially
        for idx, junc in enumerate(lsv):
            total_acum = 0.
            acum_samples = np.zeros(shape=(BINS.shape[0]))
            aggr = np.zeros(shape=(junc.shape[0]))
            for xidx, xx in enumerate(lsv):
                if idx == xidx : continue
                aggr += xx+alpha

            samples = np.ndarray(shape=(2, junc.shape[0]))
            samples[0, :] = junc + alpha
            samples[1, :] = aggr

            total_psi = np.zeros(shape=(n, BINS.shape[0]),dtype=np.float)
            for pidx, paired_samples in enumerate(samples.T):
                total_psi[pidx] = dirichlet_pdf(np.array([BINS_CENTER, 1-BINS_CENTER]).T, paired_samples)

            total_psi = np.median(total_psi, axis=0)
            psi[idx] = total_psi/total_psi.sum()

        psi_scores.append(psi)

    return psi_scores


def calc_psi(inc_samples, exc_samples, name, alpha, n, debug, psiparam):
    """
    Given a set of matching inclusion and exclusion samples, calculate
    psi, save it in disk, and return the psi-per-juntion matrix
    """
    
    print inc_samples.shape
    samples = np.vstack([inc_samples, exc_samples]).reshape(2, inc_samples.shape[0], inc_samples.shape[1])
    psi_scores = calc_dirichlet(alpha, n, samples, debug=debug, psiparam=psiparam)
    if psiparam:
        return psi_scores
    else:
        return psi_scores[:, 0]
        #psi_scores[:,1] is PSE


def mean_psi(psi_events):
    "Calculate the mean for every junction. Used for delta PSI calculation."
    ret = []
    for psi_dist in psi_events:
        #print "PSI_DIST", psi_dist
        #print "sum(PSI_DIST)", sum(psi_dist)
        ret.append(sum(psi_dist*BINS_CENTER))

    return np.array(ret)


def calc_dirichlet(alpha, n, samples_events, debug=False, psiparam=False):
    "Expects 3 dimensional matrix in samples_events"
    psi_matrix = []
    dircalc = DirichletCalc() 
    if psiparam:
        for i, event_samples in enumerate(np.rollaxis(samples_events, 1)):
        #The second dimension of the matrix corresponds to the paired samples per event (3rd dimension) for different
        # experiments (1st dimension)
            if i % 5 == 0:
                print "event %s..."%i,
                sys.stdout.flush()

            if debug > 0 and i == debug:
                break
            #if debug: print "Paired samples to dirichlet..."
            #sampling PSI by pairing the samples of the previous step sequentially
            total_acum = 0.
            acum_samples = np.array([0]*BINS)
            for h, paired_samples in enumerate(event_samples.T):
                dir_pdf = [dircalc.pdf([x, 1-x], alpha+paired_samples) for x in BINS_CENTER]
                acum_samples += dir_pdf
                total_acum += sum(dir_pdf) 

            #if debug: print "Dividing by total acum..."
            psi_matrix.append(acum_samples/total_acum)
            #print "Junction %s PSI distribution: %s sum_N: %s"%(i, psi_matrix[-1], sum(psi_matrix[-1]))
    else:
        for i, event_samples in enumerate(np.rollaxis(samples_events, 1)):
        #we iterate through the second dimension of the matrix, which corresponds to the paired samples per event
        # for different experiments
            #This is sampling the PSI. Instead of doing this, we want to fit a parametric form.
            if i % 50 == 0:
                print "event %s..." % i,
                sys.stdout.flush()

            if 0 < debug == i:
                break

            if len(event_samples.shape) == 1:
            #only one sample (this is only used if bootstrapping of reads is deactivated)
                event_psi_samples = dirichlet(event_samples+alpha, n)
            else:
                event_psi_samples = []
                #sampling PSI by pairing the samples of the previous step sequentially (to gain execution time)
                for paired_samples in event_samples.T:
                    event_psi_samples.extend(dirichlet(paired_samples+alpha, n))

            #discretization step. Get the psi samples and transform them into a histogram-like distribution
            event_psi_discrete = []
            for psi_dist in np.array(event_psi_samples).transpose():
                my_bins = list(BINS)
                my_bins.extend([1])
                #extend because:  If `bins` is a sequence,it defines the bin edges, including the rightmost edge,
                # allowing for non-uniform bin widths (form histogram docs)
                counts, limits = np.histogram(psi_dist, bins=my_bins)
                event_psi_discrete.append(counts/float(sum(counts)))

            #print event_psi_discrete[-1], len(event_psi_discrete), len(event_psi_discrete[-1])
            psi_matrix.append(event_psi_discrete)
            #print "Junction %s PSI distribution:"%i, psi_matrix[-1]

    psi_matrix = np.array(psi_matrix)
    return psi_matrix


def gen_prior_matrix(pip, lsv_exp1, lsv_exp2, output, numbins=20, defaultprior=False):

    #Start prior matrix
    pip.logger.info("Calculating prior matrix...")
    psi_space = np.linspace(0, 1-pip.binsize, num=numbins) + pip.binsize/2
    if defaultprior:
        direc = "%s/../data" % os.path.dirname(os.path.realpath(__file__))
        def_mat = pickle.load(open('%s/defaultprior.pickle' % direc, 'r'))
        prior_matrix = [def_mat, def_mat]
        return psi_space, prior_matrix

    pip.logger.info('Filtering to obtain "best set"...')

    filtered_lsv1 = majiq_filter.lsv_quantifiable(lsv_exp1, minnonzero=10, min_reads=20, logger=pip.logger)
    filtered_lsv2 = majiq_filter.lsv_quantifiable(lsv_exp2, minnonzero=10, min_reads=20, logger=pip.logger)

    ids1 = set([(xx[1], xx[2]) for xx in filtered_lsv1[1]])
    ids2 = set([(xx[1], xx[2]) for xx in filtered_lsv2[1]])
    matched_names = ids1.intersection(ids2)
    best_set_mean1 = [[], []]
    best_set_mean2 = [[], []]
    best_set_mean_ir1 = [[], []]
    best_set_mean_ir2 = [[], []]

    for ii, tt in matched_names:
        if 'i' in tt:
            continue
        for idx, nm in enumerate(filtered_lsv1[1]):
            if nm[1] == ii:
                nz = np.count_nonzero(filtered_lsv1[0][idx])
                if 'i' in nm[2]:
                    best_set_mean_ir1[0].append(nz * majiq_sample.mean_junction(filtered_lsv1[0][idx]))
                    best_set_mean_ir1[1].append(filtered_lsv1[1][idx])
                else:
                    best_set_mean1[0].append(nz * majiq_sample.mean_junction(filtered_lsv1[0][idx]))
                    best_set_mean1[1].append(filtered_lsv1[1][idx])
                break
        for idx, nm in enumerate(filtered_lsv2[1]):
            if nm[1] == ii:
                nz = np.count_nonzero(filtered_lsv2[0][idx])
                if 'i' in nm[2]:
                    best_set_mean_ir2[0].append(nz * majiq_sample.mean_junction(filtered_lsv2[0][idx]))
                    best_set_mean_ir2[1].append(filtered_lsv2[1][idx])
                else:
                    best_set_mean2[0].append(nz * majiq_sample.mean_junction(filtered_lsv2[0][idx]))
                    best_set_mean2[1].append(filtered_lsv2[1][idx])
                break

    pip.logger.info("'Best set' is %s events (out of %s)" % (len(best_set_mean1[0]), len(lsv_exp1[0])))
    best_dpsi = empirical_delta_psi(best_set_mean1[0], best_set_mean2[0])
    pip.logger.info("'Best set IR' is %s events (out of %s)" % (len(best_set_mean_ir1[0]), len(lsv_exp1[0])))
    best_dpsi_ir = empirical_delta_psi(best_set_mean_ir1[0], best_set_mean_ir2[0])

    num_priors = 1

    prior_matrix = [[], []]

    for prior_idx, best_delta_psi in enumerate((best_dpsi, best_dpsi_ir)):
        njun_prior = [[]]

        for lsv in best_delta_psi:
            if lsv.shape[0] != 2:
                continue
            njun_prior[0].append(lsv[0])

        for nj in range(len(njun_prior)):
            best_delta_psi = np.array(njun_prior[nj])

            pip.logger.info("Parametrizing 'best set'...%s",  prior_idx)
            mixture_pdf = majiq_delta.adjustdelta_lsv(best_delta_psi, output, plotpath=pip.plotpath,
                                                      title=" ".join(pip.names), numiter=pip.iter,
                                                      breakiter=pip.breakiter, njunc=nj, logger=pip.logger)
            #pickle.dump(mixture_pdf, open("%s%s_%s_bestset_junc_%s.pickle"%(output, pip.names[0], pip.names[1], nj), 'w'))

            for i in xrange(numbins):
                prior_matrix[prior_idx].extend(mixture_pdf[numbins-i:(numbins*2)-i])

            prior_matrix[prior_idx] = np.array(prior_matrix[prior_idx]).reshape(numbins, -1)
            if np.isnan(prior_matrix[prior_idx]).any():
                if prior_idx == 1:
                    pip.logger.info("Not enought statistic power to calculate the intron retention specific prior, "
                                    "in that case we will use the global prior")
                    prior_matrix[prior_idx] = prior_matrix[0]
                else:
                    raise ValueError(" The input data does not have enought statistic power in order to calculate "
                                     "the prior. Check if the input is correct or use the --default_prior option in "
                                     " order to use a precomputed prior")
            else:
                prior_matrix[prior_idx] /= sum(prior_matrix[prior_idx])
             #renormalize so it sums 1

        plot_matrix(prior_matrix[prior_idx], "Prior Matrix nj%s, version %s" % (nj, prior_idx),
                    "prior_matrix_jun_%s" % nj, pip.plotpath)
    pip.logger.info("Saving prior matrix for %s..." % pip.names)
    pickle.dump(prior_matrix, open("%s/%s_%s_priormatrix_jun_%s.pickle" % (output, pip.names[0], pip.names[1], nj),
                                   'w'))

    return psi_space, prior_matrix


def __extract_cassette(delta_psi, info, psi1, psi2):
    cas = "|1e1.1|1e2.1"
#    print "EXTRACT"
    listd = []
    for idx, linfo in enumerate(info):
#        print "linfo",linfo[2]
#        if len(linfo[2].split('|')) == 3 and linfo[2].find('e0')==-1:
        if linfo[2][1:] == cas:
            listd.append(delta_psi[idx][0])

    fp = open('./toyoseph.pickle', 'wb')
    pickle.dump(listd, fp)
    fp.close()

    out = open('./psi_delta0.tab', 'w+')
    for idx, linfo in enumerate(info):
        if linfo[2][1:] == cas:
#        if len(linfo[2].split('|')) == 3 and linfo[2].find('e0')==-1:
            lsv = delta_psi[idx]
            if 0.0125 > lsv[0] >= -0.0125:
#                print psi1, idx, np.sum(psi1[idx][0]), np.sum(psi1[0][idx])
                val1 = float(np.sum(psi1[0][idx][0])) / float(np.sum(psi1[0][idx]))
                if np.isnan(val1):
                    val1 = 0.5
                val2 = float(np.sum(psi2[0][idx][0])) / float(np.sum(psi2[0][idx]))
                if np.isnan(val2):
                    val2 = 0.5
                out.write("%d\t%d\n" % (val1, val2))

    out.close()


#deprecated
def sample_psi(psi_scores):
    """
    Input is all junctions PSI distributions for 1 replica
    """
    samples = []
    for pval, limits in psi_scores:
        event_samples = []
        sample_pos = np.multinomial(100, pval)
        for p in sample_pos:
            event_samples.append(limits[p])

        samples.append(np.mean(event_samples))

    return np.array(samples)


def main():
    """
    Script for initial testing of the MAJIQ algorithms for sampling and initial PSI values generator. 

    TODO: Many functions from this script will be extracted for general usage in the pipeline. 
    """
    parser = argparse.ArgumentParser() 
    parser.add_argument('samples', nargs='+', help='Path for samples of conditions 1 to N')
    parser.add_argument('--n', default=1, type=int, help='Number of PSI samples per sample paired') 
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution') 
    parser.add_argument('--output', required=True, help="Path to save the results to.")
    parser.add_argument('--name1', default='Inc')
    parser.add_argument('--name2', default='Exc')
    args = parser.parse_args()

    print "Loading samples..."
    samples = []
    for sample in args.samples:
        samples.append(pickle.load(open(sample)))

    samples = np.vstack(samples)
    print "Calculating PSI for %s and %s..."%(args.name1, args.name2)
    psi_scores = calc_dirichlet(args.alpha, args.n, samples)  
    pickle.dump(psi_scores, open("%s%s_vs_%s_psivalues.pickle"%(args.output, args.name1, args.name2), 'w'))
    print "Done."


if __name__ == '__main__':
    main()

