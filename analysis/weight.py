from collections import defaultdict

from pylab import *

from sklearn.preprocessing import normalize #this function uses oldnumeric, and will be deprecated in Numpy 1.9


def _kullback_lieber(p, q):
    """
    Calculates distance between two distributions: 

    Dkl(P|Q) = sum_i ln(P(i)/Q(i))*P(i)

    Using absolute values
    """
    pseudo = 0.001
    p = p + pseudo
    q = q + pseudo
    return abs((log(p/q)*p).sum(axis=1))


def local_weights(*replicas):
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
        pseudo = 0.00001
        weight_matrix = array(divergences[key])
        weight_matrix += pseudo
        weight_matrix /= weight_matrix.sum(axis=0)
        weights += (1 - weight_matrix)

    weights /= weights.sum(axis=0)

    return weights


def global_weights(*experiments):
    "Using a KL divergence, calculate the weight for every replica in a group of replicas"
    return local_weights(*experiments).sum(axis=1)
   

