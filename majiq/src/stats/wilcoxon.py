import numpy as np
from scipy.stats import ranksums
from majiq.src.stats import Stats


class Wilcoxon(Stats):

    __pval_cache__ = {}

    class Factory:
        def create(self): return Wilcoxon()

    @staticmethod
    def operator(neg, pos):
        wilcoxon_score, wilcoxon_pval = ranksums(neg, pos)
        return np.log(wilcoxon_pval)
