import numpy as np
import scipy.special as sc
from scipy.special import gammaln
from majiq.src.stats import Stats

class Infoscore(Stats):

    __info_cache__ = {}

    class Factory:
        def create(self):
            return Infoscore()

    @staticmethod
    def operator(psi, labels):
        # asort = samples.argsort()
        # psi = samples[asort]
        # labels = labels[asort]
        info_pval, info_score = Infoscore.__operator(psi, labels, assume_sorted=True)
        return info_pval

    @staticmethod
    def __lchoose(k, n):
        return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)

    @staticmethod
    def __entropy(vector):
        sm = vector.sum()
        if sm:
            sm = sc.entr(vector.astype(float) / sm).sum() * sm
        return sm

    @staticmethod
    def __compute_path_score(neg, pos, k, offset):
        left = (np.ones(2) * k + np.array([-offset, offset])) / 2
        right = np.array([neg, pos]) - left
        return Infoscore.__entropy(left) + Infoscore.__entropy(right)

    @staticmethod
    def __pvalue(neg, pos, score):
        pval = Infoscore.__info_cache__.get((neg, pos, score)) or Infoscore.__info_cache__.get((pos, neg, score))
        if pval is None:
            num_events = pos + neg
            offset = pos - neg
            old_counts = np.ones(num_events + 1) * -np.inf
            old_counts[neg] = 0
            new_counts = np.zeros(num_events + 1)
            last_max = 0
            last_min = 0
            bad_paths = -np.inf
            for k in range(1, num_events + 1):
                mn = min(last_max + 1, pos, offset + num_events - k)
                mx = max(last_min - 1, -neg, offset - (num_events - k))
                new_min = pos
                new_max = -neg
                for i in range(mx, mn + 1):
                    new_counts[i + neg] = -np.inf
                    ct = -np.inf
                    if (i - 1 >= last_min) and (i - 1 <= last_max):
                        ct = np.logaddexp(ct, old_counts[i - 1 + neg])
                    if (i + 1 >= last_min) and (i + 1 <= last_max):
                        ct = np.logaddexp(ct, old_counts[i + 1 + neg])
                    if ct > -np.inf:
                        sc = Infoscore.__compute_path_score(neg, pos, k, i)
                        if sc > score:
                            new_counts[i + neg] = ct
                            new_min = min(new_min, i)
                            new_max = max(new_max, i)
                        else:
                            dx = num_events - k
                            dy = offset - i
                            p = (dx + dy) / 2
                            assert abs(dy) <= dx
                            assert (p >= 0) and (p <= dx)
                            ct += Infoscore.__lchoose(p, dx)
                            bad_paths = np.logaddexp(bad_paths, ct)
                if new_min > new_max:
                    pval = 1.0
                    break
                last_min = new_min
                last_max = new_max
                temp = old_counts.copy()
                old_counts = new_counts.copy()
                new_counts = temp.copy()
                del temp
            else:
                pval = np.clip(np.exp(bad_paths - Infoscore.__lchoose(neg, neg + pos)), 0, 1)
            Infoscore.__info_cache__[(neg, pos, score)] = pval
            Infoscore.__info_cache__[(pos, neg, score)] = pval
        return pval

    @staticmethod
    def __operator(data, truth, assume_sorted=False, return_threshold=False):
        """INFO operator

        :param data: numpy array of values (may be sorted)
        :param truth: numpy array of boolean labels
        :param assume_sorted: if true, don't sort data
        :return: p value and score for INFO
        """
        if not assume_sorted:
            asort = np.argsort(data)
            truth = truth[asort]
            data = data[asort]
        n = data.size
        cumpos = truth.cumsum()
        pos = cumpos[-1]
        neg = n - pos
        left = np.zeros(2)
        right = np.array([neg, pos])
        best_loss = Infoscore.__entropy(right)
        threshold = data[-1] + 1
        for idx, value, cumnpos in zip(range(n), data, cumpos):
            if idx == 0 or value > data[idx - 1]:
                loss = Infoscore.__entropy(left) + Infoscore.__entropy(right - left)
                if loss < best_loss:
                    if idx == 0:
                        threshold = value - 1
                    elif idx < n - 1:
                        threshold = float(value + data[idx + 1]) / 2
                    best_loss = loss
            left = np.array([idx + 1 - cumnpos, cumnpos])
        pval = Infoscore.__pvalue(neg, pos, best_loss)
        if return_threshold:
            retval = pval, best_loss, threshold
        else:
            retval = pval, best_loss
        return retval