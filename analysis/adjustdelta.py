import pickle
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from scipy.stats import beta
from scipy.misc import logsumexp

PSEUDO = 0.0001

import pdb


def likelihood(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata):
    first = True
    pi_change = np.log(pi_change)
    pi_center = np.log(pi_center)
    for value in deltadata:
        if value > 0:
            change_prob = np.logaddexp(pi_change, beta.pdf(x=value, a=a_change, b=b_change))
            center_prob = np.logaddexp(pi_center, beta.pdf(x=value, a=a_center, b=b_center))
            all_prob = np.logaddexp(center_prob, change_prob)
            if first:
                result = all_prob
                first = False
            else:
                result = np.logaddexp(result, all_prob)

    return result


def a_from_meanvar(mean, var):
    return mean * (((mean * (1 - mean)) / var + PSEUDO) - 1)


def b_from_meanvar(mean, var):
    return (1 - mean) * (((mean * (1 - mean)) / var + PSEUDO) - 1)


def estimate_parameters(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata):
    # calculate the number of values that fall in one of the distributions
    N_change = 0
    N_center = 0
    mu_change_acum = 0
    mu_center_acum = 0
    respons = []
    #calculate responsability functions, means and N (total responsability) 
    for value in deltadata:
        change_prob = beta.pdf(x=value, a=a_change, b=b_change)
        center_prob = beta.pdf(x=value, a=a_center, b=b_center)
        #        change_prob = pi_change*beta.pdf(x=value, a=a_change, b=b_change)
        #        center_prob = pi_center*beta.pdf(x=value, a=a_center, b=b_center)

        #total probability
        total_prob = change_prob + center_prob
        #        respons_change = change_prob / total_prob
        #        respons_center = center_prob / total_prob
        respons_change = pi_change * change_prob / total_prob
        respons_center = pi_center * center_prob / total_prob
        respons.append([respons_change, respons_center])
        N_change += respons_change
        N_center += respons_center
        mu_change_acum += respons_change * value
        mu_center_acum += respons_center * value

    N_total = N_change + N_center
    pi_change = N_change / N_total
    pi_center = N_center / N_total

    mu_change = mu_change_acum / N_change
    mu_center = mu_center_acum / N_center

    #calculate variance
    acum_var_change = 0
    acum_var_center = 0
    for i, value in enumerate(deltadata):
        acum_var_change += respons[i][0] * (value - mu_change) ** 2
        acum_var_center += respons[i][1] * (value - mu_center) ** 2

    #calculate variance
    var_change = acum_var_change / N_change + PSEUDO
    var_center = acum_var_center / N_center + PSEUDO

    #calculate beta(a, b) from norm(mean, variance)
    #a_change = a_from_meanvar(mu_change, var_change)
    a_change = 1
    a_center = a_from_meanvar(mu_center, var_center)
    b_change = 1
    #b_change = b_from_meanvar(mu_change, var_change)
    b_center = b_from_meanvar(mu_center, var_center)

    return a_change, b_change, a_center, b_center, pi_change, pi_center


def label_beta(a, b, pi):
    return "(a=%.2f b=%.2f pi=%.4f)" % (a, b, pi)


def plot_all(a_center, b_center, pi_center, label_center, a_change, b_change, pi_change, label_change, figure_title,
             deltadata):
    ax = plt.subplot(2, 2, 1)
    plot_densities(deltadata, ax)
    plt.subplot(2, 2, 2)
    plot_mixture(a_center, b_center, pi_center, label_center)
    plot_mixture(a_change, b_change, pi_change, label_change)
    plt.subplot(2, 2, 3)
    plot_combpdf([[a_center, b_center, pi_center], [a_change, b_change, pi_change]])
    plt.subplot(2, 2, 4)
    plot_pi(pi_center, pi_change)
    plt.suptitle(figure_title, fontsize=24)


def plot_densities(deltadata, ax=None, my_title="Empirical Data"):
    # deltadata = nan_to_num(deltadata) #substitute nan with zero, because histogram is a shitty function that cant take nans. Shame, shame on histogram. You should be a more manly function and take NaNs without crying, you are part of matplotlib.
    ax.bar(deltadata[:, 0] - 0.0125, deltadata[:, 1], width=0.025)
    ax.set_title(my_title)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Delta PSI")
    ax.set_ylabel("Density")


# if len(deltadata[deltadata > 1]):
#        print "DELTADATA BAD", deltadata[deltadata > 1]
#        sys.exit(1)

#    values, bins = histogram(deltadata, bins = 100, range=(-1, 1))

#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    if ax:
#        ax.set_xticks([0, 0.25, 0.5, 0.75, 1]) #for cosmetics because of the z-space
#        ax.set_xticklabels([-1, -0.5, 0, 0.5, 1]) #for cosmetics because of the z-space


def truncate_betadists(beta_dists):
    "Calculate a new mean, assume variance doesn't change much"
    new_betadist = []
    x_pos = np.arange(0, 1, 0.025)  #TODO Should come from parameter
    for a, b, pi in beta_dists:
        variance = beta.stats(a, b, moments='v')
        mean_acum = 0
        for x in x_pos:
            mean_acum += beta.pdf(x=x, a=a, b=b) * x

        mean = mean_acum / len(x_pos)
        new_a, new_b = ab_from_meanvar(mean, variance)
        new_betadist.append([new_a, new_b, pi])

    return new_betadist


def calc_mixture_pdf(beta_dists):
    mixture_pdf = []
    x_pos = np.arange(0, 1, 0.025)  #TODO Should come from parameter
    for x in x_pos:
        local_sum = 0
        for a, b, pi in beta_dists:
            local_sum += beta.pdf(x=x, a=a, b=b) * pi

        mixture_pdf.append(local_sum)

    return x_pos, np.array(mixture_pdf)


def plot_combpdf(beta_dists, fig):
    fig.set_xlim(-1, 1)
    fig.set_title("Mixture PDF")
    fig.set_xlabel("PDF")
    fig.set_ylabel("Density")
    x_pos, mixture_pdf = calc_mixture_pdf(beta_dists)
    fig.plot(np.linspace(-1, 1, num=len(mixture_pdf)), mixture_pdf / 2)


def plot_mixture(a, b, pi, label_name, fig):
    fig.set_xlim(-1, 1)
    points, x_pos = calc_beta_pdf(a, b)
    fig.set_title("Beta mixtures")
    fig.set_xlabel("Delta PSI")
    fig.set_ylabel("Density")
    fig.plot(np.linspace(-1, 1, num=len(points)), points / 2, label="%s %s" % (label_name, label_beta(a, b, pi)))
    # legend()


def plot_pi(p_mixture, fig):
    fig.set_title("Pi distributions")
    fig.set_ylim(0, 1)
    fig.set_ylabel("Weight")
    fig.bar(np.arange(len(p_mixture)), p_mixture)
    #fig.set_xticks(arange(2)+0.3, ["Center", "change"], rotation=50)
    fig.set_xticks(np.arange(2) + 0.3, ["Center", "change"])


#def plot_pi(pi_center, pi_change):
#    title("Pi distributions")
#    ylim(0, 1)
#    ylabel("Weight")
#    bar(arange(2), [pi_center, pi_change])
#    xticks(arange(2)+0.3, ["Center", "change"], rotation=50)

def _save_or_show(plotpath, name):
    if plotpath:
        plt.savefig("%s_%s.png" % (plotpath, name), width=200, height=400, dpi=100)
        plt.clf()
    else:
        plt.show()


def EM(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata, num_iter, plotpath, logger=False):
    fig = plt.figure(figsize=[15, 10])
    prev_likelihood = likelihood(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata)
    if logger: logger.info("INIT: Center %s Change %s Likelihood: %.5f" % (
        label_beta(a_center, b_center, pi_center), label_beta(a_change, b_change, pi_change), prev_likelihood))

    plot_all(a_center, b_center, pi_center, "Center", a_change, b_change, pi_change, "change", "Initial Parameters",
             deltadata)
    _save_or_show(plotpath, "init")
    for iteration in xrange(num_iter):
        a_change, b_change, a_center, b_center, pi_change, pi_center = estimate_parameters(a_change, b_change, a_center,
                                                                                           b_center, pi_change,
                                                                                           pi_center, deltadata)
        current_likelihood = likelihood(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata)
        if logger: logger.info("EM iteration %s: Center %s Change %s Likelihood: %.5f" % (
            iteration, label_beta(a_center, b_center, pi_center), label_beta(a_change, b_change, pi_change),
            current_likelihood))
        plot_all(a_center, b_center, pi_center, "Center", a_change, b_change, pi_change, "Change",
                 "iteration %s (likelihood: %.5f)" % (iteration, current_likelihood), deltadata)
        _save_or_show(plotpath, "iter%05d" % iteration)
        prev_likelihood = current_likelihood

    return [[a_center, b_center, pi_center], [a_change, b_change, pi_change]]


def calc_beta_pdf(a, b, binsize=0.025):
    x_pos = np.arange(0, 1, binsize)
    beta_pdfs = []
    for x in x_pos:
        beta_pdfs.append(beta.pdf(x=x, a=a, b=b))

    return np.array(beta_pdfs), x_pos


def ab_from_meanvar(mean, var):
    return a_from_meanvar(mean, var), b_from_meanvar(mean, var)


def check_valid(a_b):
    "Make sure that init parameters are not Inf or NaN"
    ALT_INIT = 2000
    a, b = a_b
    if np.isinf(a) or np.isnan(a) or np.isinf(a) or np.isnan(b):
        return ALT_INIT, ALT_INIT
    else:
        return a, b


def adjustdelta(deltapsi, output, plotpath=None, title=None, numiter=10, breakiter=0.01, V=0.1, logger=False):
    #TODO make breakiter work
    #transform to z-space
    z_deltapsi = 0.5 * (deltapsi + 1)

    #calculate init values
    a_change = b_change = 1
    a_center = b_center = 2000

    #pi_change = len(change_delta)/float(len(z_deltapsi))
    #pi_center = len(center_delta)/float(len(z_deltapsi))
    pi_change = 0.05
    pi_center = 0.95

    if not pi_change or not pi_center:  #if any of the 'pi' parameters are 0, one of the distributions will never be considered, so reboot all pi to equal.
        pi_change = 0.05
        pi_center = 0.95

    beta_dists = EM(a_change, b_change, a_center, b_center, pi_change, pi_center, z_deltapsi, numiter, plotpath, logger)
    #truncate the beta distributions limiting them to the 0 to 1 space
    beta_dists = truncate_betadists(beta_dists)
    x_pos, z_mixture_pdf = calc_mixture_pdf(beta_dists)
    #No need to convert back from the z space, it is a distribution
    return z_mixture_pdf


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


    #    p_mixture = np.array([0.05, 0.95])
    #    beta_params = np.array([ [1 , 1], [2000, 2000]])
    #    labels = ['Uniform','Center']

    #    p_mixture = np.array([0.03, 0.03, 0.03, 0.91])
    #    beta_params = np.array([[2.5, 3],[3, 2.5],[2, 2],[2000, 2000]])
    #    labels = ['beta left', 'beta center', 'beta right','delta']



    temp = open('./temp.pickle', 'wb')
    pickle.dump((D, p_mixture, beta_params), temp)
    temp.close()

    beta_params, pmix = EMBetaMixture(D, p_mixture, beta_params, 0, logger=logger, plotpath=plotpath, nj=njunc,
                                      labels=labels)

    x_pos, z_mixture_pdf = calc_mixture_pdf_lsv(beta_params, pmix)
    return z_mixture_pdf


def plot_all_lsv(deltadata, beta_params, pmix, labels, figure_title):
    f, sp = plt.subplots(2, 2)
    plt.subplots_adjust(hspace=.4)
    #    print deltadata
    plot_densities(deltadata, sp[0, 0])
    cmb = []

    for pl in xrange(beta_params.shape[0]):
        plot_mixture(beta_params[pl, 0], beta_params[pl, 1], pmix[pl], labels[pl], sp[0, 1])
        cmb.append([beta_params[pl, 0], beta_params[pl, 1], pmix[pl]])

    plot_combpdf(cmb, sp[1, 0])
    plot_pi(pmix, sp[1, 1])
    plt.suptitle(figure_title, fontsize=24)


def loglikelihood(D, beta_mix, logp_mix, logger=False):
    N = D.shape[0]
    K = beta_mix.shape[0]
    ''' logp_DgK = log P (D | model K ) for each data point without the weight '''
    logp_DgK = np.zeros(shape=( N, K), dtype=np.float)
    logp_DgK = ma.asarray(logp_DgK)
    #print "logp_mix=%s"%(logp_mix)
    for k in xrange(K):
        #print "beta_mix[%d] = %s"%(k, beta_mix[k])
        logp_DgK[:, k] = ma.log(beta.pdf(D[:, 0], beta_mix[k, 0], beta_mix[k, 1]))

    logp_D = logp_DgK + logp_mix * np.ones(shape=(N, 1), dtype=np.float)

    dm = np.sum(logp_DgK.mask, axis=1)
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
    #transform the data to a z-space of 0-1
    if min(D[:, 0]) < 0.0: D[:, 0] = (D[:, 0] + 1) / ( c - a)
    pmix = p0_mix
    beta_mix = beta0_mix
    logp_mix = np.log(pmix)

    logp_D, logp_Dsum, LL, zrow = loglikelihood(D, beta_mix, logp_mix)
    plot_all_lsv(D0, beta_mix, pmix, labels, 'iteration 0')
    _save_or_show(plotpath, "iter_0.jun_%s" % nj)
    if logger: logger.info("[NJ:%s] Initial Log_Likelihood %.3f \n" % (nj, LL))
    #    pdb.set_trace()
    ones_1k = np.ones(shape=(1, K), dtype=np.float)
    for mm in xrange(num_iter):
        new_beta_mix = beta_mix
        new_pmix = pmix

        ''' E STEP: '''
        p_KgD = np.exp(logp_D - (logp_Dsum * ones_1k.T).T)
        p_KgD[zrow, K - 1] = 0
        #        pdb.set_trace()

        avgxPerK = np.sum(p_KgD * (D[:, 0] * ones_1k.T).T * (D[:, 1] * ones_1k.T).T, axis=0) / np.sum(p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        avgx2PerK = np.sum(p_KgD * np.square((D[:, 0] * ones_1k.T).T) * (D[:, 1] * ones_1k.T).T, axis=0) / np.sum(p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        varxPerK = avgx2PerK - (np.square(avgxPerK))

        new_beta_mix = np.zeros(shape=beta0_mix.shape, dtype=np.float)
        new_beta_mix[:, 0] = avgxPerK * (((avgxPerK * (1 - avgxPerK)) / varxPerK) - 1)
        new_beta_mix[:, 1] = (1 - avgxPerK) * (((avgxPerK * (1 - avgxPerK)) / varxPerK) - 1)

        new_pmix = np.sum(p_KgD * (D[:, 1] * ones_1k.T).T, axis=0)
        new_pmix = new_pmix / np.sum(new_pmix, axis=0)

        LLold = LL
        logp_D, logp_Dsum, LL, zrow = loglikelihood(D, new_beta_mix, np.log(new_pmix))
        if logger: logger.info("[NJ:%s] EM Iteration %d:\t LL: %.3f\n" % (nj, mm, LL))
        plot_all_lsv(D0, beta_mix, pmix, labels, 'iteration %s' % str(mm + 1))
        _save_or_show(plotpath, "iter_%05d.jun_%s" % (mm + 1, nj))

        if LL < LLold:
            if logger: logger.info("Log_Likelihood DECREASE new %d old %d - Aborting ....\n" % (LL, LLold))
            break

        pmix = new_pmix
        beta_mix = new_beta_mix
        logp_mix = np.log(pmix)

        if np.exp(LL - LLold) < (1.0 + min_ratio):
            if logger:
                logger.info("Ratio = %.3f < 1+R(%.3f) - Aborting ... \n" % (LL - LLold, min_ratio))
            break

    return beta_mix, np.array(pmix)


