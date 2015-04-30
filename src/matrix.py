"""
Functions to handle the matrices operations over the delta PSI 
"""
import numpy as np


def print_matrix(matrix):
    "Print MAJIQ delta PSI matrix in screen"
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            print "%.4f" % matrix[i][j],
        print

    print


def collapse_matrix(matrix):
    "Collapse the diagonals probabilities in 1-D and return them"
    collapse = []
    # FOR TEST matrix = array([[0, 1, 2, 3, 4, 500], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [100, 1, 2, 3, 4, 5], ])

    matrix_corner = matrix.shape[0] + 1
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def _find_delta_border(V, numbins):
    "Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have"
    delta_space = list(np.linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i
            # if nothing hit, V = 1
    return numbins


def matrix_area(matrix, V=0.2, absolute=True):
    """
    Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute
    """
    collapse = collapse_matrix(matrix)
    # get the delta psi histogram borders based on the size of 'collapse'
    border = _find_delta_border(V, collapse.shape[0])
    # grab the values inside the area of interest
    area = []
    if V < 0:
        area.append(collapse[0:border + 1])
        if absolute:  #if absolute V, pick the other side of the array
            area.append(collapse[-border - 1:])
    else:
        area.append(collapse[border:])
        if absolute:  #if absolute V, pick the other side of the array
            area.append(collapse[0:len(collapse) - border])

    return sum(area)


def matrix_prob_e(matrix):
    """
    Calculate sum_v v*P(Delta PSI > V)
    """
    absolute = True
    ret = 0.
    for v in np.arange(0, 1, 0.1):
        ret += matrix_area(matrix, V=v, absolute=absolute) * v

    return ret


def matrix_e(matrix):
    "Calculates the expected value of delta PSI E()"
    collapse = collapse_matrix(matrix)  # one dimensional discrete distribution of psi
    delta_space = list(
        np.linspace(-1, 1, num=len(collapse)))  # the values that correspond to the different delta psi [-1, 1]
    e = 0
    for i, value in enumerate(collapse):
        e += value * delta_space[i]

    return e


def rank_deltas(matrices, names, V=0.2, absolute=True, E=False, ranknochange=False):
    "Rank all deltas in an event by level of change. V sets a threshold for change, E overrides V and calculates an average of V values"
    rank = []
    for i, dmatrix in enumerate(matrices):
        if E:
            v_prob = matrix_prob_e(dmatrix)
            rank.append([names[i], v_prob])
        else:
            area = matrix_area(dmatrix, V, absolute)
            if ranknochange:  # P(Delta PSI < V) = 1 - P(Delta PSI > V)
                area = 1 - area

            rank.append([names[i], area])

    rank.sort(key=lambda x: -x[1])
    return rank


def rank_deltas_lsv(matrices, names, V=0.2, absolute=True, E=False, ranknochange=False):
    "Rank all deltas in an event by level of change. V sets a threshold for change, E overrides V and calculates an average of V values"
    rank = []
    for lidx, lsv in enumerate(matrices):
        v_prob = []
        import pdb
        # pdb.set_trace()
        if len(lsv) > 2: continue
        for dmatrix in lsv:
            if E:
                v_prob.append(matrix_prob_e(dmatrix))

            else:
                area = matrix_area(dmatrix, V, absolute)
                if ranknochange:  # P(Delta PSI < V) = 1 - P(Delta PSI > V)
                    area = 1 - area
                v_prob.append(area)
        rank.append([names[lidx][1], max(v_prob)])

    rank.sort(key=lambda x: -x[1])
    return rank


def rank_empirical_delta_lsv(list_events, info):
    rank = []
    for lidx, lsv in enumerate(list_events):
        rank.append([info[lidx][1], max(lsv)])

    rank.sort(key=lambda x: -x[1])
    return rank