# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


def all_pairwise(list_1, list_2, order_matters=False):
    """
    Given two lists, return list of lists of all
        possible pairwise comparisons.

    If List_1 == List_2, this function is effectively performing:
        List CHOOSE/COMBINATION 2, if Order_Matters = False
        List PICK/PERMUTATION 2, if Order_Matters = True

    Returns [[a,b],[a,c],[etc,etc]] of combinations/permutations
    """
    if type(list_1) is not list:
        raise ValueError("List_1 needs to be a list...")

    if type(list_2) is not list:
        raise ValueError("List_2 needs to be a list...")

    if type(order_matters) is not bool:
        raise ValueError("Order_Matters needs to be True or False...")

    tuples = [(x, y) for x in list_1 for y in list_2 if x != y]
    for entry in tuples:
        if not order_matters:
            if (entry[1], entry[0]) in tuples:
                tuples.remove((entry[1], entry[0]))

    matrix = list()

    for t in tuples:
        matrix.append(list(t))
    return matrix
