from collections import defaultdict

from pylab import *


def _kullback_lieber(p, q):
    """
    Calculates distance between matrices with paired distributions, one distance per row (event normally)

    Dkl(P|Q) = sum_i ln(P(i)/Q(i))*P(i)
    """
    pseudo = 1./sys.maxint
    pq = log(p+pseudo) - log(q+pseudo)
    return (pq*p).sum(axis=1)


def local_weights(replicas):
    """
    Using a KL divergence, calculate the weight for every event in every replica in a group of replicas
    """
    divergences = defaultdict(list)
    for i, replica_i in enumerate(replicas):
        for j, replica_j in enumerate(replicas):
            if i == j:
                divergences[i].append([0.]*len(replica_i)) #for performance
            else:
                divergences[i].append(_kullback_lieber(replica_i, replica_j))

    #with all the KL divergences per event and pair, acumulate weights 
    weights = zeros(shape=array(divergences[0]).shape) #generates a matrix with zeroes to acumulate al the local KL divergences
    for key in divergences.keys():
        pseudo = 0.000000000000000001
        weight_matrix = array(divergences[key])
        weight_matrix += pseudo
        weight_matrix /= weight_matrix.sum(axis=0)
        weights += (1 - weight_matrix)

    weights /= weights.sum(axis=0)
    return weights


def global_weights(experiments=None, locweights=None):
    "Using a KL divergence, calculate the weight for every replica in a group of replicas"
    if locweights == None: 
        locweights = local_weights(experiments)

    return locweights.sum(axis=1) / locweights.shape[1] 
    

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
    """
    print "WEIGHTS (a VS a)"
    print local_weights([a, a])
    print global_weights([a, a])
    print
    """
    print "WEIGHTS (a VS b)"
    print local_weights([a, b])
    print global_weights([a, b])
    print
    
    print "WEIGHTS (a VS a VS b)"
    print local_weights([a, a, b])
    print global_weights([a, a, b])
    print 
    print "WEIGHTS (a VS a VS a VS b)"
    print local_weights([a, a, a, b])
    print global_weights([a, a, a, b])
    