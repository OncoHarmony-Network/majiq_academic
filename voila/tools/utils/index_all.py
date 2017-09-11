# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


def index_all(array,
              element,
              nottt=False):
    """Return all indices of array that point to an element.
        Arguments:
            array   = type: list.
            element = object that you wish to get indices for from array.
            Not     = If True, return indices of elements that do *not* match
    """
    if not isinstance(array, list):
        raise Exception("Please only give lists to this function.")
    if nottt:
        matched_indices = [i for i, x in enumerate(array) if x != element]
    else:
        matched_indices = [i for i, x in enumerate(array) if x == element]
    return matched_indices
