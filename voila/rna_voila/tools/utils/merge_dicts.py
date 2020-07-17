# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    https://stackoverflow.com/a/26853961/7378802
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result
