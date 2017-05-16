import numpy as np
import scipy.special as sc
from majiq.src.stats import Stats


class Tnom(Stats):

    __pval_cache__ = {}

    class Factory:
        def create(self): return Tnom()

    @staticmethod
    def operator(psi, labels):#psi, labels):
        # nsamps = csamps[clabels == 0]
        # psamps = csamps[clabels == 1]
        tnom_pval, tnom_score = Tnom.__operator(psi, labels, assume_sorted=True)
        return np.log(tnom_pval)

    @staticmethod
    def __lpathn(l, o):
        k = (l + abs(o)) // 2
        return sc.gammaln(l + 1) - sc.gammaln(k + 1) - sc.gammaln(l - k + 1)

    @staticmethod
    def __pvalue(neg, pos, score):
        """
        Arguments: neg, pos, score
        Returns p-value for TNOM score given neg negative and pos positive labels.
        """
        pval = self.__pval_cache__.get((neg, pos, score)) or self.__pval_cache__.get((pos, neg, score))
        if pval is None:
            if (pos <= score) or (neg <= score):
                pval = 1.
            else:
                l = pos + neg
                o = pos - neg
                nums = -np.inf * np.ones(2)
                delta = 2 * (np.array([pos, neg]) - score)
                cts = delta.copy()
                sn = 0
                while abs(cts[1 - sn] + o) <= l or abs(cts[sn] - o) <= l:
                    nums[sn] = np.logaddexp(nums[sn], Tnom.__lpathn(l, cts[sn] - o))
                    nums[sn] = np.logaddexp(nums[sn], Tnom.__lpathn(l, cts[1 - sn] + o))
                    sn ^= 1
                    cts[0] += delta[sn]
                    cts[1] += delta[1 - sn]
                npath = Tnom.__lpathn(l, o)
                pval = np.exp(nums[0] - npath) - np.exp(nums[1] - npath)
                pval = np.clip(pval, 0., 1.)
            Tnom.__pval_cache__[(neg, pos, score)] = pval
            Tnom.__pval_cache__[(pos, neg, score)] = pval
        return pval

    @staticmethod
    def __operator(data, truth, assume_sorted=False, return_threshold=False):
        """TNOM operator

        :param data: numpy array of values (may be sorted)
        :param truth: numpy array of boolean labels
        :param assume_sorted: if true, don't sort data
        :return: p value and score for TNOM
        """
        if not assume_sorted:
            asort = np.argsort(data)
            truth = truth[asort]
            data = data[asort]
        n = data.size
        cumpos = truth.cumsum()
        pos = cumpos[-1]
        neg = n - pos
        best_loss = min(neg, pos)
        threshold = data[-1] + 1
        for idx, value, cumnpos in zip(xrange(n), data, cumpos):
            if idx > 0 and value == data[idx - 1]:
                continue
            loss = idx + 1 + pos - cumnpos * 2
            loss = min(loss, n - loss)
            if loss < best_loss:
                if idx == 0:
                    threshold = value - 1
                elif idx < n - 1:
                    threshold = float(value + data[idx + 1]) / 2
                best_loss = loss
        pval = Tnom.__pvalue(neg, pos, best_loss)
        if return_threshold:
            retval = pval, best_loss, threshold
        else:
            retval = pval, best_loss
        return retval

