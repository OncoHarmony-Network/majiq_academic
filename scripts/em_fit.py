import sys
import argparse
import pickle
from pylab import *
import warnings
warnings.filterwarnings('error')

from scipy.stats import beta

"""
TODO: Use the Beta distribution with 2 parameters and then fit them in the 3 mixed distributions.
"""

def _add_pi_pdf(pi, pdf):
    if pdf > 0:
        return logaddexp(pi, log(pdf))
    else:
        return 0

def likelihood(loc_left, loc_right, loc_center, a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center, deltadata):
    first = True
    pi_left = log(pi_left)
    pi_right = log(pi_right)
    pi_center = log(pi_center)
    for value in deltadata:
        if value > 0:

            left_prob = _add_pi_pdf(pi_left, beta.pdf(x=value, a=a_left, b=b_left, loc=loc_left))
            right_prob = logaddexp(pi_right, beta.pdf(x=value, a=a_right, b=b_right, loc=loc_right))
            center_prob = logaddexp(pi_center, beta.pdf(x=value, a=a_center, b=b_center, loc=loc_center))
            all_prob = logaddexp(left_prob, right_prob)
            all_prob = logaddexp(all_prob, center_prob)
            if first:
                result = all_prob
                first = False
            else:
                result = logaddexp(result, all_prob)


    return result

def mu_transform(mu, a, c):
    return (mu - a)/(c - a)    

def var_transform(var, a, c):
    return var/(c - a)**2  

def a_from_meanvar(mean, var):
    return mean*(((mean*(1-mean))/var)-1)

def b_from_meanvar(mean, var):
    return (1-mean)*(((mean*(1-mean))/var)-1)

def estimate_parameters(a_left, b_left, a_right, b_right, a_center, b_center, loc_left, loc_right, loc_center, pi_left, pi_right, pi_center, deltadata):
    #calculate the number of values that fall in one of the distributions
    N_left = 0
    N_right = 0
    N_center = 0
    mu_left_acum = 0
    mu_right_acum = 0
    mu_center_acum = 0
    respons = []
    for value in deltadata:
        left_prob = pi_left*beta.pdf(x=value, a=a_left, b=b_left, loc=loc_left)
        right_prob = pi_right*beta.pdf(x=value, a=a_right, b=b_right, loc=loc_right)
        center_prob = pi_center*beta.pdf(x=value, a=a_center, b=b_center, loc=loc_center)
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
    #end calculate var

    #mu and var on the left and center are not in the 0,1 interval, replace
    mu_left = mu_transform(mu_left, loc_left, loc_left+1)
    mu_center = mu_transform(mu_center, loc_center, loc_center+1)
    var_left = var_transform(var_left, loc_left, loc_left+1)
    var_center = var_transform(var_center, loc_center, loc_center+1)

    a_left = a_from_meanvar(mu_left, var_left)
    a_right = a_from_meanvar(mu_right, var_right)
    a_center = a_from_meanvar(mu_center, var_center) 
    b_left = b_from_meanvar(mu_left, var_left)
    b_right = b_from_meanvar(mu_right, var_right)
    b_center = b_from_meanvar(mu_center, var_center)  

    return a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center

def EM(a_left, b_left, a_right, b_right, a_center, b_center, loc_left, loc_right, loc_center, pi_left, pi_right, pi_center, deltadata, num_iter):
    #direction = True, bigger values, false, smaller
    prev_likelihood = likelihood(loc_left, loc_right, loc_center, a_left, b_left, 
                                 a_right, b_right, a_center, b_center, 
                                 pi_left, pi_right, pi_center, deltadata)
    print "INIT: Left (a=%.2f b=%.2f pi=%.4f) Center (a=%.2f b=%.2f pi=%.4f) Right (a=%.2f b=%.2f pi=%.4f) Likelihood: %.5f"%(a_left, b_left, pi_left, a_center, b_center, pi_center, a_right, b_right, pi_right, prev_likelihood)
    for iteration in xrange(num_iter):
        print "EM iteration %s..."%(iteration)
        a_left, b_left, a_right, b_right, a_center, b_center, pi_left, pi_right, pi_center = estimate_parameters(a_left, b_left, a_right, b_right, a_center, b_center, loc_left, loc_right, loc_center, pi_left, pi_right, pi_center, deltadata)
        current_likelihood = likelihood(loc_left, loc_right, loc_center, 
                                        a_left, b_left, a_right, b_right, a_center, b_center, 
                                        pi_left, pi_right, pi_center, deltadata)
        print "Left (a=%.2f b=%.2f pi=%.4f) Center (a=%.2f b=%.2f pi=%.4f) Right (a=%.2f b=%.2f pi=%.4f) Likelihood: %.5f"%(a_left, b_left, pi_left, a_center, b_center, pi_center, a_right, b_right, pi_right, current_likelihood)
        

        prev_likelihood = current_likelihood

    return a_left, b_left, a_right, b_right, a_center, b_center


def calc_beta_pdf(a, b, loc):
    x_pos = arange(loc, loc+1, 0.01)
    beta_pdfs = []
    for x in x_pos:
        beta_pdfs.append(beta.pdf(x=x, a=a, b=b, loc=loc))    

    return array(beta_pdfs), x_pos

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('deltapsi', help='Path for pickle with the deltapsi values')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--title', default=None, help='') 
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    parser.add_argument('--iter', default=20, type=int, help='Max number of iterations of the EM')
    parser.add_argument('--breakiter', default=0.01, type=float, help='If the log likelihood increases')
    args = parser.parse_args()
    deltapsi = pickle.load(open(args.deltapsi))

    subplot(2,1,1)
    xlim(-1, 1)
    hist(deltapsi, bins = 60, histtype='step')    

    ######  PARAMETERS    #####
    #loc values are the start of the interval of the interval 0,1 (default 0) in a beta distribution for the numpy function
    loc_left   = -1 
    loc_center = -0.5
    loc_right  = 0
    a_left = b_left = a_center = b_center = a_right = b_right = 2 #A and B of the beta distribution
    pi_left = pi_right = pi_center = 1./3

    a_left, b_left, a_right, b_right, a_center, b_center = EM(a_left, b_left, a_right, b_right, a_center, b_center, loc_left, loc_right, loc_center, pi_left, pi_right, pi_center, deltapsi, args.iter)

    center_points, x_pos_center = calc_beta_pdf(a_center, b_center, loc_center)
    right_points, x_pos_right  = calc_beta_pdf(a_right, b_right, loc_right)
    left_points, x_pos_left = calc_beta_pdf(a_left, b_left, loc_left)

    subplot(2,1,2)

    plot(x_pos_center, center_points, label="center")
    plot(x_pos_right, right_points, label="right")
    plot(x_pos_left, left_points, label="left")
    xlim(-1, 1)
    legend()
    show()

if __name__ == '__main__':
    main()

