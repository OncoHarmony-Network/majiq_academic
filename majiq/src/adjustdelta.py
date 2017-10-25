from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta
from scipy.misc import logsumexp
import majiq.src.plotting as mplot
PSEUDO = 1e-300


def calc_mixture_pdf(beta_dists):
    mixture_pdf = []
    x_pos = np.arange(0, 1, 0.025)  #TODO Should come from parameter
    for x in x_pos:
        local_sum = 0
        for a, b, pi in beta_dists:
            local_sum += beta.pdf(x=x, a=a, b=b) * pi

        mixture_pdf.append(local_sum)

    return x_pos, np.array(mixture_pdf)



def calc_beta_pdf(a, b, binsize=0.025):
    x_pos = np.arange(0, 1, binsize)
    beta_pdfs = []
    for x in x_pos:
        beta_pdfs.append(beta.pdf(x=x, a=a, b=b))

    return np.array(beta_pdfs), x_pos

def check_valid(a_b):
    "Make sure that init parameters are not Inf or NaN"
    ALT_INIT = 2000
    a, b = a_b
    if np.isinf(a) or np.isnan(a) or np.isinf(a) or np.isnan(b):
        return ALT_INIT, ALT_INIT
    else:
        return a, b


def calc_mixture_pdf_lsv(beta_param, pmix):
    mixture_pdf = []
    x_pos = np.arange(0, 1, 0.025)  #TODO Should come from parameter
    for x in x_pos:
        local_sum = 0
        for ii, bt in enumerate(beta_param):
            local_sum += beta.pdf(x=x, a=bt[0], b=bt[1]) * pmix[ii]
        mixture_pdf.append(local_sum)

    return x_pos, np.array(mixture_pdf)


def calculate_beta_params(mean, vari, sample_size):
    b = (1 - mean) * sample_size
    a = mean * sample_size

    return [a, b]


def adjustdelta_lsv(deltapsi, output, plotpath=None, title=None, numiter=10, breakiter=0.01, V=0.1, njunc=1,
                    logger=False):
    D = np.zeros(shape=(79, 2), dtype=np.float)
    xpos = np.arange(-1 + (0.025 / 2), 1, 0.025)

    for idx, ii in enumerate(xpos[:-1]):
        D[idx, 0] = round((ii + xpos[idx + 1]) / 2, 5)
        if D[idx, 0] == 0:
            zero_idx = idx

    for ppv in deltapsi:
        for idx, ii in enumerate(xpos[:-1]):
            if ii <= ppv < xpos[idx + 1]:
                D[idx, 1] += 1
                break

    total = D[:, 1].sum()
    num_spike = D[zero_idx, 1]
    num_lowcenter = D[zero_idx - 3:zero_idx + 3, 1].sum() - num_spike

    p_mixture = np.zeros(shape=3, dtype=np.float)

    spike = calculate_beta_params(0.5, 0, num_spike)
    p_mixture[2] = num_spike / total
    #    pdb.set_trace()
    center = calculate_beta_params(0.5, 0, num_lowcenter)
    p_mixture[1] = num_lowcenter / total
    uniform = [1, 1]
    p_mixture[0] = 1 - ((num_lowcenter + num_spike) / total)
    beta_params = np.array([uniform, center, spike])

    labels = ['Uniform', 'center', 'spike']

    beta_params, pmix = EMBetaMixture(D, p_mixture, beta_params, 0, logger=logger, plotpath=plotpath, nj=njunc,
                                      labels=labels)

    x_pos, z_mixture_pdf = calc_mixture_pdf_lsv(beta_params, pmix)
    return z_mixture_pdf



def loglikelihood(D, beta_mix, logp_mix, logger=False):

    N = D.shape[0]
    K = beta_mix.shape[0]
    ''' logp_DgK = log P (D | model K ) for each data point without the weight '''
    logp_DgK = np.zeros(shape=(N, K), dtype=np.float)

    for k in range(K):
        logp_DgK[:, k] = np.log(beta.pdf(D[:, 0], beta_mix[k, 0], beta_mix[k, 1]) + PSEUDO)

    logp_D = logp_DgK + logp_mix * np.ones(shape=(N, 1), dtype=np.float)

    dm = np.sum(logp_DgK, axis=1)
    zrow = dm.astype(np.bool)
    no_zrow = np.logical_not(zrow)

    logp_Dsum = np.zeros(shape=(N), dtype=np.float)
    logp_Dsum[no_zrow] = logsumexp(logp_D[no_zrow, :], axis=1)
    logp_Dsum[zrow] = logsumexp(logp_D[zrow, :-1], axis=1)

    LL = np.sum(logp_Dsum * D[:, 1], axis=0)

    return logp_D, logp_Dsum, LL, zrow


def EMBetaMixture(D, p0_mix, beta0_mix, num_iter, min_ratio=1e-5, logger=False, plotpath=None, nj=0, labels=None):
    D0 = D.copy()
    N = D.shape[0]
    K = beta0_mix.shape[0]

    c = 1
    a = -1

    if min(D[:, 0]) < 0.0: D[:, 0] = (D[:, 0] + 1) / (c - a)
    pmix = p0_mix
    beta_mix = beta0_mix
    logp_mix = np.log(pmix)

    logp_D, logp_Dsum, LL, zrow = loglikelihood(D, beta_mix, logp_mix)
    mplot.plot_all_lsv(D0, beta_mix, pmix, labels, 'iteration 0')
    mplot.save_or_show(plotpath, "iter_0.jun_%s" % nj)
    if logger:
        logger.debug("[NJ:%s] Initial Log_Likelihood %.3f \n" % (nj, LL))

    ones_1k = np.ones(shape=(1, K), dtype=np.float)
    for mm in range(num_iter):
        new_beta_mix = beta_mix
        new_pmix = pmix

        ''' E STEP: '''
        p_KgD = np.exp(logp_D - (logp_Dsum * ones_1k.T).T)
        p_KgD[zrow, K - 1] = 0
        #        pdb.set_trace()

        avgxPerK = np.sum(p_KgD * (D[:, 0] * ones_1k.T).T * (D[:, 1] * ones_1k.T).T, axis=0) / np.sum(
            p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        avgx2PerK = np.sum(p_KgD * np.square((D[:, 0] * ones_1k.T).T) * (D[:, 1] * ones_1k.T).T, axis=0) / np.sum(
            p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        varxPerK = avgx2PerK - (np.square(avgxPerK))

        new_beta_mix = np.zeros(shape=beta0_mix.shape, dtype=np.float)
        new_beta_mix[:, 0] = avgxPerK * (((avgxPerK * (1 - avgxPerK)) / varxPerK) - 1)
        new_beta_mix[:, 1] = (1 - avgxPerK) * (((avgxPerK * (1 - avgxPerK)) / varxPerK) - 1)

        new_pmix = np.sum(p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        new_pmix = new_pmix / np.sum(new_pmix, axis=0)

        LLold = LL
        logp_D, logp_Dsum, LL, zrow = loglikelihood(D, new_beta_mix, np.log(new_pmix))
        if logger: logger.debug("[NJ:%s] EM Iteration %d:\t LL: %.3f\n" % (nj, mm, LL))
        mplot.plot_all_lsv(D0, beta_mix, pmix, labels, 'iteration %s' % str(mm + 1))
        mplot.save_or_show(plotpath, "iter_%05d.jun_%s" % (mm + 1, nj))

        if LL < LLold:
            if logger:
                logger.debug("Log_Likelihood DECREASE new %d old %d - Aborting ....\n" % (LL, LLold))
            break

        pmix = new_pmix
        beta_mix = new_beta_mix
        logp_mix = np.log(pmix)

        if np.exp(LL - LLold) < (1.0 + min_ratio):
            if logger:
                logger.debug("Ratio = %.3f < 1+R(%.3f) - Aborting ... \n" % (LL - LLold, min_ratio))
            break

    return beta_mix, np.array(pmix)


