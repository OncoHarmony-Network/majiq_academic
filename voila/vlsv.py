import numpy as np
import scipy.interpolate as sinter



class Het:
    def __init__(self):
        """
        Het stat data. 
        :param experiment_names: List of experiment names. 
        """
        super(Het, self).__init__()
        self.groups = []
        self.junction_stats = []

    def add_group(self, expected_psi, median):
        """
        Add per LSV het group stats. 
        :param expected_psi: Expected psi 2d array, where x axis is experiments and y axis is junctions
        :param median: Median psi bins
        :param group_name: The name of the group being added
        :return: 
        """
        self.groups.append(HetGroup(expected_psi, median))

    def add_junction_stats(self, stats):
        """
        Add per junction het stat data.
        :param stats: 2d array of group comparison stats 
        :param stat_name: Name of stat
        :param junction_id: Junction id
        :return: 
        """
        self.junction_stats.append(stats)

    def cls_dict(self):
        return {'groups': HetGroup}

    def __iter__(self):
        for k, v in self.__dict__.items():
            yield k, v


class HetGroup:
    def __init__(self, expected_psi, median):
        super(HetGroup, self).__init__()
        self.expected_psi = expected_psi
        self.median = median

    def get_psi(self, experiment_index, junction_index):
        return self.expected_psi[experiment_index][junction_index]

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None).from_hdf5(h)

    def __iter__(self):
        for k, v in self.__dict__.items():
            yield k, v


def get_expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1 + 1. / len(bins), 1., 2. / len(bins)))


def get_expected_psi(bins):
    bins = np.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return np.sum(projection_prod)


def collapse_matrix(matrix):
    """
    Collapse the diagonals probabilities in 1-D and return them
    """
    collapse = []
    matrix = np.array(matrix)
    matrix_corner = matrix.shape[0]
    for i in range(-matrix_corner + 1, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def _find_delta_border(V, numbins):
    """
    Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have
    :param V:
    :param numbins:
    :return:
    """
    delta_space = list(np.linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i
    # if nothing hit, V = 1
    return numbins


def matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
    """
    Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute.
    :param collapsed_mat:
    :param V:
    :param absolute:
    :param matrix:
    :return:
    """
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    collapse = np.concatenate(([0], collapse))
    collapse = np.cumsum(collapse)
    xbins = np.linspace(-1, 1, num=collapse.size)
    if absolute:
        Vabs = abs(V)
        left, right = np.interp([-Vabs, Vabs], xbins, collapse, left=0, right=1)
        area = left + (1 - right)
    else:
        area = np.interp(V, xbins, collapse, left=0, right=1)
        if V >= 0:
            area = 1 - area
    return area


def matrix_area_spline(matrix, V=0.2, absolute=True, collapsed_mat=False, kind=2):
    """
    Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute.
    :param collapsed_mat: if True, matrix must be 1d, else matrix must be 2d square
    :param V: threshold
    :param absolute: if True, calculate two-tailed area
    :param matrix: list or array of probabilities, must sum to 1
    :param kind: parameter kind from scipy.interpolate.interp1d
    :return: probability of delta psi exceeding a certain threshold
    """
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    collapse = np.concatenate(([0], collapse))
    collapse = np.cumsum(collapse)
    xbins = np.linspace(-1, 1, num=collapse.size)
    spline = sinter.interp1d(xbins, collapse, kind=kind, fill_value=(0, 1), bounds_error=False)
    if absolute:
        Vabs = abs(V)
        left, right = spline([-Vabs, Vabs])
        area = left + (1 - right)
    else:
        area = spline(V)
        if V >= 0:
            area = 1 - area
    return float(area)


def is_lsv_changing(means, threshold):
    """
    Return true if lsv is changing based on threshold.
    :param threshold: lsv threshold value
    :return: bool
    """
    means = np.array(tuple(means))
    means_gt_zero = means[means > 0]
    means_sum = means_gt_zero.sum()
    max_value = max(means_sum, abs(means_sum))
    return max_value >= threshold
