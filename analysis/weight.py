from collections import defaultdict

from pylab import *



def _l1(p, q):
    return (abs(p - q)).sum(axis=1)

def _kullback_lieber(p, q):
    """
    Calculates distance between matrices with paired distributions, one distance per row (event normally)

    Dkl(P|Q) = sum_i ln(P(i)/Q(i))*P(i)
    """
    pseudo = 1./sys.maxint
    pq = log(p+pseudo) - log(q+pseudo)
    return (pq*p).sum(axis=1)


def _local_distance(a, b, l1):
        if l1:
            return _l1(a, b)
        else:
            return _kullback_lieber(a, b)    

def local_weights(replicas, l1=False, median_ref=array([])):
    """
    Using a KL divergence, calculate the weight for every event in every replica in a group of replicas
    """
    pseudo = 0.00000000000000000000001
    distances = []
    if median_ref.size: #if we have a median reference, compare all 
        for i, replica_i in enumerate(replicas):
            distances.append(_local_distance(replica_i, median_ref, l1))   
        
        distances = array(distances)
        distances += pseudo
        weights = (1 - distances)
                
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
    """ 
    print "WEIGHTS (a VS a VS a VS b)"
    print local_weights([a, a, a, b], median_ref=a)
    #print global_weights([a, a, a, b])
    