import sys
import argparse
import pickle
from pylab import *

from scipy.stats import beta

"""
TODO: Use the Beta distribution with 2 parameters and then fit them in the 3 mixed distributions.
"""

def EM(loc, a, b, deltapsi_slice, num_iter):
    #direction = True, bigger values, false, smaller
    names = ["A", "B"]
    parameters = [a, b]
    precision = [20., 20.]

    for param in xrange(len(parameters)): #adjust all parameters one by one
        lower = [False, False] #if its lower to the left and to the right respectively  
        direction = False #start looking to the left
        prev_likelihood = likelihood(loc, parameters[0], parameters[1], deltapsi_slice)   
        prev_value = parameters[param]
        for iteration in xrange(num_iter):
            print "EM iteration %s..."%iteration
            print "Adjusting %s: value: %s Likelihood: %s precision: %s Direction: %s"%(names[param], parameters[param], prev_likelihood, precision[param], direction)
            if direction:
                parameters[param] += precision[param]
            else:
                parameters[param] -= precision[param]

            current_likelihood = likelihood(loc, parameters[0], parameters[1], deltapsi_slice)       

            if current_likelihood < prev_likelihood or isnan(current_likelihood): 
                #if we decrease the likelihood, rewind to previous value
                parameters[param] = prev_value
                current_likelihood = prev_likelihood
                #indicate that we failed in that direction
                if direction:
                    lower[1] = True
                else:
                    lower[0] = True

                if lower[0] and lower[1]: # if we have looked in both directions, lower precision, keep direction
                    precision[param] /= 2
                    lower = [False, False]
                elif lower[0]: #if we failed left, change to right
                    direction = True
                else: #if we failed right, change to left
                    direction = False

            else:
                lower = [False, False]

            prev_likelihood = current_likelihood
            prev_value = parameters[param]

        print "BEST EM: Loc: %s A: %s B: %s. Likelihood: %s\n"%(loc, a, b, current_likelihood)
    return loc, parameters[0], parameters[1]


def EM_center(loc, a, b, deltapsi_slice, num_iter):
    """
    This one is expected to be almost symetrical, so we increase and decrease a and b together 
    """
    lower = [False, False] #if its lower to the left and to the right respectively  
    direction = True #start looking to the right
    precision = 20.
    prev_likelihood = likelihood(loc, a, b, deltapsi_slice)   
    prev_a = a
    prev_b = b
    for iteration in xrange(num_iter):
        print "EM iteration %s..."%iteration
        print "Adjusting %s: A: %s B: %s Likelihood: %s precision: %s Direction: %s"%("CENTER", a, b, prev_likelihood, precision, direction)
        if direction:
            a += precision
            b += precision
        else:
            a -= precision
            b -= precision

        current_likelihood = likelihood(loc, a, b, deltapsi_slice)       

        if current_likelihood < prev_likelihood or isnan(current_likelihood): 
            #if we decrease the likelihood, rewind to previous value
            a = prev_a
            b = prev_b
            current_likelihood = prev_likelihood
            #indicate that we failed in that direction
            if direction:
                lower[1] = True
            else:
                lower[0] = True

            if lower[0] and lower[1]: # if we have looked in both directions, lower precision, keep direction
                precision /= 2
                lower = [False, False]
            elif lower[0]: #if we failed left, change to right
                direction = True
            else: #if we failed right, change to left
                direction = False

        else:
            lower = [False, False]

        prev_likelihood = current_likelihood
        prev_a = a
        prev_b = b

    print "BEST EM: Loc: %s A: %s B: %s. Likelihood: %s\n"%(loc, a, b, current_likelihood)
    return loc, a, b



def likelihood(loc, a, b, deltadata):
    beta_pdfs = []
    first = True
    for value in deltadata:
        if first:
            result = log(beta.pdf(x=value, a=a, b=b, loc=loc))
        else:
            result = logaddexp(result, log(beta.pdf(x=value, a=a, b=b, loc=loc)))

    return result

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
    parser.add_argument('--iter', default=5, type=int, help='Number of iterations of the EM')
    args = parser.parse_args()
    deltapsi = pickle.load(open(args.deltapsi))

    subplot(2,1,1)
    xlim(-1, 1)
    hist(deltapsi, bins = 60, histtype='step')    


    """Divide the deltapsi values in the three groups we want to fit a beta distribution to: left, center and right. 

    """
    left_limit = -0.1
    right_limit = 0.1

    deltapsi_left = deltapsi[deltapsi < left_limit]
    deltapsi_right = deltapsi[deltapsi > right_limit]
    #TODO this 2 steps can be done in one, dont remember how
    deltapsi_center = deltapsi[deltapsi > left_limit]
    deltapsi_center = deltapsi_center[deltapsi_center < right_limit]
    #END TODO

    ######  PARAMETERS    #####
    #loc values are the center of the interval 0,1 (default 0) in a beta distribution for the numpy function
    loc_left   = -1 
    loc_center = -0.5
    loc_right  = 0
    a_left = b_left = a_center = b_center = a_right = b_right = 5 #A and B of the beta distribution


    loc_left, a_left, b_left = EM(loc_left, a_left, b_left, deltapsi_left, args.iter)
    loc_right, a_right, b_right = EM(loc_right, a_right, b_right, deltapsi_right, args.iter)
    loc_center, a_center, b_center = EM_center(loc_center, a_center, b_center, deltapsi_center, args.iter)

    center_points, x_pos_center = calc_beta_pdf(a_center, b_center, loc_center)
    right_points, x_pos_right  = calc_beta_pdf(a_right, b_right, loc_right)
    left_points, x_pos_left = calc_beta_pdf(a_left, b_left, loc_left)

    subplot(2,1,2)

    plot(x_pos_center, center_points, label="nochange")
    plot(x_pos_right, right_points, label="right")
    plot(x_pos_left, left_points, label="left")
    xlim(-1, 1)
    legend()
    show()

if __name__ == '__main__':
    main()

