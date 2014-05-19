import sys
import argparse
import pickle
from pylab import *

from scipy.stats import beta

PSEUDO = 0.0001

def likelihood(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata):
    first = True
    pi_change = log(pi_change)
    pi_center = log(pi_center)
    for value in deltadata:
        if value > 0:
            change_prob = logaddexp(pi_change, beta.pdf(x=value, a=a_change, b=b_change))
            center_prob = logaddexp(pi_center, beta.pdf(x=value, a=a_center, b=b_center))
            all_prob = logaddexp(center_prob, change_prob)
            if first:
                result = all_prob
                first = False
            else:
                result = logaddexp(result, all_prob)

    return result

def a_from_meanvar(mean, var):
    return mean*(((mean*(1-mean))/var+PSEUDO)-1)

def b_from_meanvar(mean, var):
    return (1-mean)*(((mean*(1-mean))/var+PSEUDO)-1)

def estimate_parameters(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata):
    #calculate the number of values that fall in one of the distributions
    N_change = 0
    N_center = 0
    mu_change_acum = 0
    mu_center_acum = 0
    respons = []
    #calculate responsability functions, means and N (total responsability) 
    for value in deltadata:
        change_prob = pi_change*beta.pdf(x=value, a=a_change, b=b_change)
        center_prob = pi_center*beta.pdf(x=value, a=a_center, b=b_center)
        #total probability
        total_prob = change_prob + center_prob
        respons_change = change_prob / total_prob
        respons_center = center_prob / total_prob
        respons.append([respons_change, respons_center])
        N_change += respons_change
        N_center += respons_center
        mu_change_acum += respons_change*value
        mu_center_acum += respons_center*value

    N_total = N_change + N_center
    pi_change = N_change/N_total
    pi_center = N_center/N_total

    mu_change = mu_change_acum/N_change
    mu_center = mu_center_acum/N_center

    #calculate variance
    acum_var_change = 0
    acum_var_center = 0
    for i, value in enumerate(deltadata):
        acum_var_change += respons[i][0]*(value - mu_change)**2
        acum_var_center += respons[i][1]*(value - mu_center)**2

    #calculate variance
    var_change = acum_var_change/N_change+PSEUDO
    var_center = acum_var_center/N_center+PSEUDO   

    #calculate beta(a, b) from norm(mean, variance)
    #a_change = a_from_meanvar(mu_change, var_change)
    a_change = 1
    a_center = a_from_meanvar(mu_center, var_center) 
    b_change = 1
    #b_change = b_from_meanvar(mu_change, var_change)
    b_center = b_from_meanvar(mu_center, var_center)  

    return a_change, b_change, a_center, b_center, pi_change, pi_center

def label_beta(a, b, pi):
    return "(a=%.2f b=%.2f pi=%.4f)"%(a, b, pi)


def plot_all(a_center, b_center, pi_center, label_center, a_change, b_change, pi_change, label_change, figure_title, deltadata):
    ax = subplot(2,2,1)
    plot_densities(deltadata, ax)
    subplot(2,2,2)
    plot_mixture(a_center, b_center, pi_center, label_center)
    plot_mixture(a_change, b_change, pi_change, label_change)
    subplot(2,2,3)
    plot_combpdf([[a_center, b_center, pi_center], [a_change, b_change, pi_change]])   
    subplot(2,2,4)
    plot_pi(pi_center, pi_change)
    suptitle(figure_title, fontsize=24)


def plot_densities(deltadata, ax = None, my_title="Empirical Data"):
    deltadata = nan_to_num(deltadata) #substitute nan with zero, because histogram is a shitty function that cant take nans. Shame, shame on histogram. You should be a more manly function and take NaNs without crying, you are part of matplotlib.
    title(my_title)
    xlim(0, 1)
    xlabel("Delta PSI")
    ylabel("Density")
    if len(deltadata[deltadata > 1]):
        print "DELTADATA BAD", deltadata[deltadata > 1]
        sys.exit(1)

    values, bins = histogram(deltadata, bins = 100, range=(-1, 1))
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    bar(center, values, align='center', width=width)
    if ax:
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1]) #for cosmetics because of the z-space
        ax.set_xticklabels([-1, -0.5, 0, 0.5, 1]) #for cosmetics because of the z-space


def truncate_betadists(beta_dists):
    "Calculate a new mean, assume variance doesn't change much"
    new_betadist = []
    x_pos = arange(0, 1, 0.025) #TODO Should come from parameter
    for a, b, pi in beta_dists:
        variance = beta.stats(a, b, moments='v')
        mean_acum = 0
        for x in x_pos:
            mean_acum += beta.pdf(x=x, a=a, b=b)*x

        mean = mean_acum / len(x_pos)
        new_a, new_b = ab_from_meanvar(mean, variance)
        new_betadist.append([new_a, new_b, pi])

    return new_betadist


def calc_mixture_pdf(beta_dists):
    mixture_pdf = []
    x_pos = arange(0, 1, 0.025) #TODO Should come from parameter
    for x in x_pos:
        local_sum = 0
        for a, b, pi in beta_dists:
            local_sum += beta.pdf(x=x, a=a, b=b)*pi

        mixture_pdf.append(local_sum)    

    return x_pos, array(mixture_pdf)

def plot_combpdf(beta_dists):
    xlim(-1, 1)
    title("Mixture PDF")
    xlabel("PDF")
    ylabel("Density")
    x_pos, mixture_pdf = calc_mixture_pdf(beta_dists)
    plot(linspace(-1, 1, num=len(mixture_pdf)), mixture_pdf / 2)

def plot_mixture(a, b, pi, label_name):
    xlim(-1, 1)
    points, x_pos = calc_beta_pdf(a, b)
    title("Beta mixtures")
    xlabel("Delta PSI")
    ylabel("Density")
    plot(linspace(-1, 1, num=len(points)), points / 2, label="%s %s"%(label_name, label_beta(a, b, pi)))
    legend()

def plot_pi(pi_center, pi_change):
    title("Pi distributions")
    ylim(0, 1)
    ylabel("Weight")
    bar(arange(2), [pi_center, pi_change])
    xticks(arange(2)+0.3, ["Center", "change"], rotation=50)

def _save_or_show(plotpath, name):
    if plotpath:
        savefig("%s_%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  


def EM(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata, num_iter, plotpath, logger=False):
    fig = figure(figsize=[15, 10])
    prev_likelihood = likelihood(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata)
    if logger: logger.info("INIT: Center %s Change %s Likelihood: %.5f"%(label_beta(a_center, b_center, pi_center), label_beta(a_change, b_change, pi_change), prev_likelihood))
    
    plot_all(a_center, b_center, pi_center, "Center", a_change, b_change, pi_change, "change", "Initial Parameters", deltadata)
    _save_or_show(plotpath, "init")    
    for iteration in xrange(num_iter):
        a_change, b_change, a_center, b_center, pi_change, pi_center = estimate_parameters(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata)
        current_likelihood = likelihood(a_change, b_change, a_center, b_center, pi_change, pi_center, deltadata)
        if logger: logger.info("EM iteration %s: Center %s Change %s Likelihood: %.5f"%(iteration, label_beta(a_center, b_center, pi_center), label_beta(a_change, b_change, pi_change), current_likelihood))        
        plot_all(a_center, b_center, pi_center, "Center", a_change, b_change, pi_change, "Change", "iteration %s (likelihood: %.5f)"%(iteration, current_likelihood), deltadata)
        _save_or_show(plotpath, "iter%05d"%iteration)
        prev_likelihood = current_likelihood

    return [[a_center, b_center, pi_center], [a_change, b_change, pi_change]]


def calc_beta_pdf(a, b, binsize=0.025):
    x_pos = arange(0, 1, binsize)
    beta_pdfs = []
    for x in x_pos:
        beta_pdfs.append(beta.pdf(x=x, a=a, b=b))

    return array(beta_pdfs), x_pos

def ab_from_meanvar(mean, var):
    return a_from_meanvar(mean, var), b_from_meanvar(mean, var)

def check_valid(a_b):
    "Make sure that init parameters are not Inf or NaN"
    ALT_INIT = 2000
    a, b = a_b
    if isinf(a) or isnan(a) or isinf(a) or isnan(b):
        return ALT_INIT, ALT_INIT
    else:
        return a, b


def adjustdelta(deltapsi, output, plotpath=None, title=None, numiter=10, breakiter=0.01, V=0.1, logger=False):
    #TODO make breakiter work
    #transform to z-space
    z_deltapsi = 0.5*(deltapsi+1)

    #calculate init values
    a_change = b_change = 1
    a_center = b_center = 2000

    #pi_change = len(change_delta)/float(len(z_deltapsi))
    #pi_center = len(center_delta)/float(len(z_deltapsi))
    pi_change = 0.05
    pi_center = 0.95

    if not pi_change or not pi_center: #if any of the 'pi' parameters are 0, one of the distributions will never be considered, so reboot all pi to equal.
        pi_change = 0.05
        pi_center = 0.95

    beta_dists = EM(a_change, b_change, a_center, b_center, pi_change, pi_center, z_deltapsi, numiter, plotpath, logger)
    #truncate the beta distributions limiting them to the 0 to 1 space
    beta_dists = truncate_betadists(beta_dists) 
    x_pos, z_mixture_pdf = calc_mixture_pdf(beta_dists)
    #No need to convert back from the z space, it is a distribution
    return z_mixture_pdf 


def adjustdelta_lsv( deltapsi, output, plotpath=None, title=None, numiter=10, breakiter=0.01, V=0.1, logger=False ) :
    #TODO make breakiter work
    #transform to z-space
    z_deltapsi = 0.5*(deltapsi+1)

    #calculate init values
    a_change = b_change = 1
    a_center = b_center = 2000

    #pi_change = len(change_delta)/float(len(z_deltapsi))
    #pi_center = len(center_delta)/float(len(z_deltapsi))
    pi_change = 0.05
    pi_center = 0.95

    if not pi_change or not pi_center: #if any of the 'pi' parameters are 0, one of the distributions will never be considered, so reboot all pi to equal.
        pi_change = 0.05
        pi_center = 0.95

    beta_dists = EM(a_change, b_change, a_center, b_center, pi_change, pi_center, z_deltapsi, numiter, plotpath, logger)
    #truncate the beta distributions limiting them to the 0 to 1 space
    beta_dists = truncate_betadists(beta_dists) 
    x_pos, z_mixture_pdf = calc_mixture_pdf(beta_dists)
    #No need to convert back from the z space, it is a distribution
    return z_mixture_pdf 

