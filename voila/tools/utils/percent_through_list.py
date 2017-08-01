import math

from voila.tools.utils.calebs_xrange import calebs_xrange


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


def percent_through_list(the_list,
                         update_at_this_perc=0.1):
    """
    Given a List (or length of the list), determine which
        indices in the list corrsepond to By percent along the list.

        for example, List = [A,B,C,D,E,F,G,H,I,J] (length = 10)
        if By = 0.2, the dictionary would be:
        {0: 0.0, 2: 20.0, 4: 40.0, 6: 60.0, 8: 80.0, 10: 100.0}

    Use this function to help you determine how many intervals of a given
        percentage along a list loop you are...
    """
    if isinstance(the_list, list):
        list_n = len(the_list)
    elif isinstance(the_list, int):
        list_n = the_list
        # print "Assuming you gave len(List) instead of the list itself..."
    else:
        raise ValueError("List needs to be a list (or the length of one ...)")
    if not isinstance(update_at_this_perc, float) or update_at_this_perc <= 0 or update_at_this_perc >= 1:
        raise ValueError("By needs to be float >0 and <1")
    # get items at each By percent interval
    by_percent_index = dict()
    for by in calebs_xrange(0, 1, update_at_this_perc):
        index = int(math.floor(by * list_n))
        by_percent_index[index] = by * 100.0
    return by_percent_index
