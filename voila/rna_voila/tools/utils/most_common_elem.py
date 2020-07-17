import itertools
import operator


def most_common(the_list):
    """
    return most common element of list

    http://stackoverflow.com/a/1520716
    """
    # get an iterable of (item, iterable) pairs
    SL = sorted((x, i) for i, x in enumerate(the_list))
    # print 'SL:', SL
    groups = itertools.groupby(SL, key=operator.itemgetter(0))

    # auxiliary function to get "quality" for an item
    def _auxfun(g):
        item, iterable = g
        count = 0
        min_index = len(the_list)
        for _, where in iterable:
            count += 1
            min_index = min(min_index, where)
        # print 'item %r, count %r, minind %r' % (item, count, min_index)
        return count, -min_index

    # pick the highest-count/earliest item
    return max(groups, key=_auxfun)[0]
