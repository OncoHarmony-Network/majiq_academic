from collections import defaultdict
from scipy.stats import beta
import analysis.psi as majiq_psi
import numpy as np
from pylab import *


def _l1(p, q):

    res = abs(p-q)
    res = res / 2.0
    return  res.sum(axis=0)
    #    return (abs(p - q)).sum(axis=0)

def _kullback_lieber(p, q):
    """
    Calculates distance between matrices with paired distributions, one distance per row (event normally)

    Dkl(P|Q) = sum_i ln(P(i)/Q(i))*P(i)
    """
    """
    pseudo = 0
    p = p+pseudo #add pseudocount
    p /= p.sum(axis=0) #normalize to 1
    q = q+pseudo #add pseudocount  
    q /= q.sum(axis=0) #normalize to 1
    print "P", p
    print p.shape
    print p.sum(axis=1)
    print "Q", q
    print q.shape
    print q.sum(axis=1)
    """
    return (log(p / q)*p).sum(axis=1)


def _local_distance(a, b, l1, inv=False):
    if l1:
        res =  _l1(a, b)
        if inv: res = 1 -res
    else:
        res = _kullback_lieber(a, b)
    return res








def local_weight_eta_nu_per_lsv ( group, nexp ):

    alpha = bta = 0.5
    x = np.array([majiq_psi.BINS_CENTER, 1-majiq_psi.BINS_CENTER]).T
    jeffreys = majiq_psi.dirichlet_pdf(x,np.array([alpha,bta]))

    max_eta_idx = np.zeros(shape=(len(group), nexp), dtype=np.int)
    eta = np.zeros(shape=(len(group), nexp), dtype=np.float)
    median_psi = np.zeros(shape=(len(group), nexp), dtype=np.float)
    print " SIZES nexx", nexp, "num_lsv", len(group)
    for lidx, lsv in enumerate(group):
        for eidx, exp_lsv in enumerate(lsv):
            mpsi = []
            max_junc_eta = 0
            for jidx in xrange(exp_lsv.shape[0]):
                eta_e = np.zeros(shape=(exp_lsv.shape[1]), dtype=np.float)
                for sample in xrange(exp_lsv.shape[1]):
                    total = exp_lsv[:,sample].sum()
                    cov   = exp_lsv[jidx, sample]
                    smpl  = np.array([cov, total-cov]) + 0.5
                    mpsi.append( majiq_psi.dirichlet_pdf( x , smpl) )
                    eta_e[sample] = _local_distance( mpsi[-1], jeffreys, l1 = True )

                d_eta =  eta_e.sum()/float(exp_lsv.shape[1]) 
                if max_junc_eta < d_eta: 
                    max_junc_eta = d_eta
                    max_eta_idx[lidx, eidx] = jidx

            median_psi[lidx, eidx] = median(array(mpsi[max_eta_idx[lidx, eidx]]), axis=0)
            eta[lidx, eidx] = max_junc_eta


    nu = np.zeros(shape=(len(group), nexp), dtype=np.float)
    for lidx, lsv_exp in enumerate(group):
        for eidx, lsv in enumerate(lsv_exp):
            nu_e = np.zeros(shape=(lsv.shape[1]), dtype=np.float)
            jidx = max_eta_idx[lidx, eidx]
            for sample in xrange(lsv.shape[1]):
                    total = lsv[:,sample].sum()
                    cov = lsv[jidx, sample]
                    smpl = np.array([cov, total-cov]) + 0.5
                    psi = majiq_psi.dirichlet_pdf(x, smpl)
                    nu_e[sample] = _local_distance(psi, median_psi[lidx, eidx], l1 = True, inv=True)
            nu[lidx, eidx] = (nu_e.sum()/float(lsv.shape[1]))

    return eta, nu



def local_weight_eta_nu ( group, nexp ):

    alpha = bta = 0.5
    x = np.array([majiq_psi.BINS_CENTER, 1-majiq_psi.BINS_CENTER]).T
    jeffreys = majiq_psi.dirichlet_pdf(x,np.array([alpha,bta]))

    max_eta_idx = np.zeros(shape=(len(group), nexp), dtype=np.int)
    eta = np.zeros(shape=(len(group), nexp), dtype=np.float)
    median_psi = np.zeros(shape=(len(group), nexp), dtype=np.float)
    print " SIZES nexx", nexp, "num_lsv", len(group)

    # loop over the experiments
    for lidx, lsv in enumerate(group):
        #loop over the LSVs
        for eidx, exp_lsv in enumerate(lsv):
            junc_eta  = []
            junc_mpsi = []
            mpsi = np.zeros(shape = (exp_lsv.shape[0],majiq_psi.nbins), dtype=np.float)
            #loop over the junctions
            for jidx in xrange(exp_lsv.shape[0]):
                eta_e = 0.0
                for sample in xrange(exp_lsv.shape[1]):
                    total = exp_lsv[:,sample].sum()
                    cov   = exp_lsv[jidx, sample]
                    smpl  = np.array([cov, total-cov]) + 0.5
                    mpsi[jidx] = majiq_psi.dirichlet_pdf( x , smpl)
                    eta_e += _local_distance( mpsi[-1], jeffreys, l1 = True )

                d_eta = eta_e/float(exp_lsv.shape[1]) 
                junc_eta.append( d_eta )
                d_eta = median(mpsi, axis=0)
                junc_mpsi.append( d_eta )

            median_psi[lidx, eidx] = junc.mpsi 
            eta[lidx, eidx] = junc_eta

#    import pdb
#    pdb.set_trace()

    # loop over the experiments
    nu = np.zeros(shape=(len(group), nexp), dtype=np.float)
    for lidx, lsv_exp in enumerate(group):
        #loop over the LSVs
        for eidx, lsv in enumerate(lsv_exp):
            junc_nu = []
            #loop over the junctions
            for jidx in xrange(exp_lsv.shape[0]):
                nu_e = 0.0
                for sample in xrange(lsv.shape[1]):
                    total = lsv[:,sample].sum()
                    cov = lsv[jidx, sample]
                    smpl = np.array([cov, total-cov]) + 0.5
                    psi = majiq_psi.dirichlet_pdf(x, smpl)
                    nu_e += _local_distance(psi, median_psi[lidx, eidx][jidx], l1 = True, inv=True)
                d_nu = nu_e / float(lsv.shape[1])
                junc_nu.append( d_nu )
            nu[lidx, eidx] = junc_nu

    return eta, nu

def global_weight_ro ( group , num_exp, n=100 ):

    ''' Calculate psi for global_weight '''
    
    median_ref = [] 
    #ev = [None]*len(group)
    ev = []
    psis = np.zeros(shape=(num_exp, len(group)), dtype=np.dtype('object'))

    for eidx in xrange( num_exp):
        lsv_list = group[:,eidx]
        psi = majiq_psi.lsv_psi( lsv_list,  0.5, n, 0 )
        psis[eidx,:] = psi

    for lidx in xrange( len(group)):
        tmp = np.array([xx[0] for xx in psis[:,lidx]])
        median_ref.append(np.median(tmp, axis=0))

    ro = np.zeros( shape=(num_exp), dtype= np.float)
    ro_lsv = np.zeros(shape=num_exp, dtype=np.float)
    for eidx, exp in enumerate(psis):
        for lidx, psi in enumerate(exp):
            ro_lsv[eidx] += _local_distance(psi[0],median_ref[lidx], l1 = True)

    ro = ro_lsv / ro_lsv.sum()
    return ro

def local_weights(replicas, l1=False, median_ref=array([])):
    """
    Using either L1 or DKL, calculate the weight for every event in every replica in a group of replicas
    """
    pseudo = 0.00000000000000000000001
    distances = []
    if median_ref.size: #if we have a median reference, compare all 
        for i, replica_i in enumerate(replicas):
            distances.append(_local_distance(replica_i, median_ref, l1))   
        
        distances = array(distances)
        distances += pseudo
        distances /= distances.sum(axis=0)
        #print "DISTANCES", distances
        #print distances.shape
        weights = (1 - distances)
        #print "WEIGHTS", weights
        weights /= weights.sum(axis=0)
        #print "WEIGHTS NORM", weights

    else:
        divergences = defaultdict(list)
        for i, replica_i in enumerate(replicas):
            for j, replica_j in enumerate(replicas):
                if i == j:
                    divergences[i].append([0.]*len(replica_i)) #for performance
                else:
                    divergences[i].append(_local_distance(replica_i, replica_j, l1))  

        #with all the divergences per event and pair, acumulate weights 
        weights = zeros(shape=array(divergences[0]).shape) #generates a matrix with zeroes to acumulate al the local KL/L1 divergences
        for key in divergences.keys():
            weight_matrix = array(divergences[key])
            weight_matrix += pseudo
            weight_matrix /= weight_matrix.sum(axis=0)
            weights += (1 - weight_matrix)

        weights /= weights.sum(axis=0)
    
    return weights


def global_weights(experiments=None, locweights=None, l1=False):
    "Using a KL divergence, calculate the weight for every replica in a group of replicas"
    if locweights == None: 
        locweights = local_weights(experiments, l1)

    return locweights.sum(axis=1) / locweights.shape[1] 
    1

if __name__ == '__main__':
    #some simple tests
    #zero
    #print _kullback_lieber(array([0.2,0.3,0.4,0.1]), array([0.2,0.3,0.4,0.1]))
    #non symetric big values
    #print _kullback_lieber(array([0.2,0.3,0.4,0.1]), array([0.1,0.1,0.1,0.7]))
    #print _kullback_lieber(array([0.1,0.1,0.1,0.7]), array([0.2,0.3,0.4,0.1]))    
    #symetric small values
    #print _kullback_lieber(array([0.2,0.3,0.4,0.1]), array([0.3,0.2,0.4,0.1]))
    #print _kullback_lieber(array([0.3,0.2,0.4,0.1]), array([0.2,0.3,0.4,0.1]))

    a = array([[0.2,0.3,0.4,0.1], [0.2,0.3,0.4,0.1], [0.2,0.3,0.4,0.1], [0.2,0.3,0.4,0.1]])
    b = array([[0.8,0.1,0.1,0.0], [0.1,0.3,0.4,0.2], [0.2,0.3,0.4,0.1], [0.2,0.3,0.4,0.1]])
    c = array([[0.8,0.1,0.1,0.0], [0.8,0.1,0.1,0.0], [0.8,0.1,0.1,0.0], [0.8,0.1,0.1,0.0]])

    print "A"
    print a
    print "B"
    print b
    print "C"
    print c
    print
    print "WEIGHTS (a VS a)"
    print "----------------"
    print "- DKL"
    print local_weights([a, a])
    print global_weights([a, a])
    print "- L1"
    print local_weights([a, a], l1=True)
    print global_weights([a, a], l1=True)
    print

    print "WEIGHTS (a VS b)"
    print "- DKL"
    print local_weights([a, b])
    print global_weights([a, b])
    print "- L1"    
    print local_weights([a, b], l1=True)
    print global_weights([a, b], l1=True)
    print
    
    print "WEIGHTS (a VS a VS b)"
    print local_weights([a, a, b])
    print global_weights([a, a, b])
    print

    print "WEIGHTS (a VS a VS a VS b VS c)"
    print local_weights([a, a, a, b, c])
    print
    print local_weights([a, a, a, b, c], median_ref=a)
    #print global_weights([a, a, a, b])
    
