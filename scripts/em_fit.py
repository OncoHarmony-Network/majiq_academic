import sys
import argparse
import pickle
from pylab import *
#import warnings
#warnings.filterwarnings('error')

from scipy.stats import beta

"""
TODO: Use the Beta distribution with 2 parameters and then fit them in the 3 mixed distributions.

Do the following:

1. Transform your original values from the range X \in [-1,1] to Z \in [0,1]. with the simple 
transformation Z = 0.5*(X+1). You can now work only in the transformed space with all the mixture components 
and that will simplify the code and make it faster.

2. Initialize the 3 mixtures in the following way: take the %V right most values 
(where V is a global param, say 10 or 20), compute their mean and variance and use these 
to set the a,b for the right component. Do the same for the left component. 
Take the remaining 1-2V% points that are more in the center and do the same to initialize 
the middle component. You now also have an initial guess for the proportion of Pi i.e. 
the prior for each mixture component (V% for each side one and 100-2V for the middle one). 
Alternatively, you could use the same technique you used here where the variances are 
the same but the means are shifted for each component and you start of by having Pi = 1/3 for each of them.
 This is simpler. Probably will give very similar/same results, but may take longer to converge. 

3. Run the EM now on the modified Z space, do the plots and everything in this space.

4. At the end, you will need to go back to the X space between -1 and 1. The transformation 
between your mixture PDF in the Z space and the X space is simple: P(X) = P(Z)/2. 
So, when you get a value in X space, transform it to Z space, compute the PDF 
from the mixture of beta distributions, and divide the result by 2. That will be the returned value.

Other comments:
The mixture plots are nice, but these only give the shape of each distribution separately. 
You should plot (a) the combined PDF as we discussed, where you plot P(x) for  x = -1:0.1:+1 
by the mixture distribution and (b) plot the Pi i.e. the probabilitiy for each mixture 
component at each step e.g. as a bar chart.

"""


def likelihood(a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center, deltadata):
    first = True
    pi_left = log(pi_left)
    pi_right = log(pi_right)
    pi_center = log(pi_center)
    for value in deltadata:
        if value > 0:
            left_prob = logaddexp(pi_left, beta.pdf(x=value, a=a_left, b=b_left))
            right_prob = logaddexp(pi_right, beta.pdf(x=value, a=a_right, b=b_right))
            center_prob = logaddexp(pi_center, beta.pdf(x=value, a=a_center, b=b_center))
            all_prob = logaddexp(left_prob, right_prob)
            all_prob = logaddexp(all_prob, center_prob)
            if first:
                result = all_prob
                first = False
            else:
                result = logaddexp(result, all_prob)


    return result

def a_from_meanvar(mean, var):
    return mean*(((mean*(1-mean))/var)-1)

def b_from_meanvar(mean, var):
    return (1-mean)*(((mean*(1-mean))/var)-1)

def estimate_parameters(a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center, deltadata):
    #calculate the number of values that fall in one of the distributions
    N_left = 0
    N_right = 0
    N_center = 0
    mu_left_acum = 0
    mu_right_acum = 0
    mu_center_acum = 0
    respons = []
    #calculate responsability functions, means and N 
    for value in deltadata:
        left_prob = pi_left*beta.pdf(x=value, a=a_left, b=b_left)
        right_prob = pi_right*beta.pdf(x=value, a=a_right, b=b_right)
        center_prob = pi_center*beta.pdf(x=value, a=a_center, b=b_center)
        #total probability
        total_prob = left_prob + right_prob + center_prob

        respons_left = left_prob / total_prob
        respons_right = right_prob / total_prob
        respons_center = center_prob / total_prob

        respons.append([respons_left, respons_right, respons_center])

        N_left += respons_left
        N_right += respons_right
        N_center += respons_center

        mu_left_acum += respons_left*value
        mu_right_acum += respons_right*value
        mu_center_acum += respons_center*value

    N_total = N_left + N_right + N_center
    pi_left = N_left/N_total
    pi_right = N_right/N_total
    pi_center = N_center/N_total

    mu_left = mu_left_acum/N_left
    mu_right = mu_right_acum/N_right
    mu_center = mu_center_acum/N_center

    #calculate var
    acum_var_left = 0
    acum_var_right = 0
    acum_var_center = 0
    for i, value in enumerate(deltadata):
        acum_var_left += respons[i][0]*(value - mu_left)**2
        acum_var_right += respons[i][1]*(value - mu_right)**2
        acum_var_center += respons[i][2]*(value - mu_center)**2

    var_left = acum_var_left/N_left
    var_right = acum_var_right/N_right
    var_center = acum_var_center/N_center    

    a_left = a_from_meanvar(mu_left, var_left)
    a_right = a_from_meanvar(mu_right, var_right)
    a_center = a_from_meanvar(mu_center, var_center) 
    b_left = b_from_meanvar(mu_left, var_left)
    b_right = b_from_meanvar(mu_right, var_right)
    b_center = b_from_meanvar(mu_center, var_center)  

    return a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center

def label_beta(a, b, pi):
    return "(a=%.2f b=%.2f pi=%.4f)"%(a, b, pi)


def plot_all(a_left, b_left, pi_left, label_left, a_center, b_center, pi_center, label_center, a_right, b_right, pi_right, label_right, figure_title, deltadata):
    plot_densities(deltadata)
    plot_mixture(a_left, b_left, pi_left, label_left)
    plot_mixture(a_center, b_center, pi_center, label_center)
    plot_mixture(a_right, b_right, pi_right, label_right)
    plot_combpdf([[a_left, b_left, pi_left], [a_center, b_center, pi_center], [a_right, b_right, pi_right]])   
    plot_pi(pi_left, pi_center, pi_right)
    suptitle(figure_title, fontsize=24)


def plot_densities(deltadata):
    subplot(2,2,1)
    title("Empirical Data")
    xlim(0, 1)
    xlabel("Delta PSI")
    ylabel("Density")
    hist(deltadata, bins = 40, histtype='step')      


def plot_combpdf(beta_dists):
    subplot(2,2,3)
    xlim(0, 1)
    mixture_pdf = []
    title("Mixture PDF")
    xlabel("PDF")
    ylabel("Density")
    x_pos = arange(0, 1, 0.01)
    for x in x_pos:
        local_sum = 0
        for a, b, pi in beta_dists:
            local_sum += beta.pdf(x=x, a=a, b=b)*pi

        mixture_pdf.append(local_sum)
    plot(x_pos, mixture_pdf)

def plot_mixture(a, b, pi, label_name):
    subplot(2,2,2)
    xlim(0, 1)
    points, x_pos = calc_beta_pdf(a, b)
    title("Beta mixtures")
    xlabel("Delta PSI")
    ylabel("Density")
    plot(x_pos, points, label="%s %s"%(label_name, label_beta(a, b, pi)))
    legend()

def plot_pi(pi_left, pi_center, pi_right):
    subplot(2,2,4)
    title("Pi distributions")
    ylim(0, 1)
    ylabel("Weight")
    bar(arange(3), [pi_left, pi_center, pi_right])
    xticks(arange(3)+0.3, ["Left", "Center", "Right"], rotation=50)

def _save_or_show(plotpath, name):
    if plotpath:
        savefig("%s_%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  


def EM(a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center, deltadata, num_iter, plotpath):
    fig = figure(figsize=[15, 10])
    prev_likelihood = likelihood(a_left, b_left, 
                                 a_right, b_right, a_center, b_center, 
                                 pi_left, pi_right, pi_center, deltadata)
    print "INIT: Left %s Center %s Right %s Likelihood: %.5f"%(label_beta(a_left, b_left, pi_left), label_beta(a_center, b_center, pi_center), label_beta(a_right, b_right, pi_right), prev_likelihood)
    
    plot_all(a_left, b_left, pi_left, "Left", a_center, b_center, pi_center, "Center", a_right, b_right, pi_right, "Right", "Initial Parameters", deltadata)
    _save_or_show(plotpath, "init")    
    for iteration in xrange(num_iter):
        a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center = estimate_parameters(a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center, deltadata)
        current_likelihood = likelihood(a_left, b_left, a_right, b_right, a_center, b_center, 
                                        pi_left, pi_right, pi_center, deltadata)
        print "EM iteration %s: Left %s Center %s Right %s Likelihood: %.5f"%(iteration, label_beta(a_left, b_left, pi_left), label_beta(a_center, b_center, pi_center), label_beta(a_right, b_right, pi_right), current_likelihood)        

        plot_all(a_left, b_left, pi_left, "Left", a_center, b_center, pi_center, "Center", a_right, b_right, pi_right, "Right", "iteration %s (likelihood: %.5f)"%(iteration, current_likelihood), deltadata)
        _save_or_show(plotpath, "iter%05d"%iteration)   

        prev_likelihood = current_likelihood

    return a_left, b_left, a_right, b_right, a_center, b_center


def calc_beta_pdf(a, b):
    x_pos = arange(0, 1, 0.01)
    beta_pdfs = []
    for x in x_pos:
        beta_pdfs.append(beta.pdf(x=x, a=a, b=b))    

    return array(beta_pdfs), x_pos

def ab_from_meanvar(mean, var):
    return a_from_meanvar(mean, var), b_from_meanvar(mean, var)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('deltapsi', help='Path for pickle with the deltapsi values')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--title', default=None, help='') 
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    parser.add_argument('--iter', default=10, type=int, help='Max number of iterations of the EM')
    parser.add_argument('--breakiter', default=0.01, type=float, help='If the log likelihood increases less that this flag, do not do another EM step')
    parser.add_argument('--V', default=0.1, type=float, help='Value of DeltaPSI used for initialization of the EM model')
    args = parser.parse_args()
    deltapsi = pickle.load(open(args.deltapsi))

    #transform to z-space
    z_deltapsi = 0.5*(deltapsi+1)
    z_right_value = 0.5*(args.V+1)
    z_left_value = 1-z_right_value

    #calculate init values
    right_delta = z_deltapsi[z_deltapsi > z_right_value]
    left_delta = z_deltapsi[z_deltapsi < z_left_value]
    center_delta = z_deltapsi[z_deltapsi > z_left_value]
    center_delta = center_delta[center_delta < z_right_value]

    a_left, b_left = ab_from_meanvar(mean(left_delta), var(left_delta))
    a_right, b_right = ab_from_meanvar(mean(right_delta), var(right_delta))
    a_center, b_center = ab_from_meanvar(mean(center_delta), var(center_delta))

    pi_left = len(left_delta)/float(len(z_deltapsi))
    pi_right = len(right_delta)/float(len(z_deltapsi))
    pi_center = len(center_delta)/float(len(z_deltapsi))

    a_left, b_left, a_right, b_right, a_center, b_center = EM(a_left, b_left, a_right, b_right, a_center, b_center, 
                                                              pi_left, pi_right, pi_center, 
                                                              z_deltapsi, args.iter, args.plotpath)

    #TODO truncate matrixes by calculating the mean and variance from 0 to 1 for every beta(a, b) and use a_from_meanvar() and b_from_meanvar()
    #TODO convert back to 1-:1 space and pickle save
    arange(-98.75, 100, 2.5)


if __name__ == '__main__':
    main()

