import numpy as np


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


def get_expected_value(bins, left=0, right=1):
    step = (right - left) / len(bins)
    return np.arange(left + step / 2, right, step).dot(bins)


def get_expected_dpsi(bins):
    return get_expected_value(bins, left=-1)


def get_expected_psi(bins):
    return get_expected_value(bins)


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


def matrix_area(matrix, threshold=0.2, non_changing=False):
    """
    Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is
    absolute.

    :param non_changing:
    :param threshold:
    :param matrix:
    :return:
    """

    collapse = np.cumsum(matrix)
    collapse /= collapse[-1]
    xbins = np.linspace(-1, 1, num=collapse.size)
    abs_threshold = abs(threshold)
    left, right = np.interp([-abs_threshold, abs_threshold], xbins, collapse, left=0, right=1)
    area = left + (1 - right)
    area = np.clip(area, 0, 1)
    if non_changing:
        return 1 - area
    return area


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
