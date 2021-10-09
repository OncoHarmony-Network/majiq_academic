/**
 * TNOM.hpp
 *
 * Statistics relative to Total Number of Mistakes.
 * See Ben-Dor, Friedman, Yakhini 2002 ("Overabundance Analysis and Class
 * Discovery in Gene Expression Data")
 * for details on exact computation of statistic p-value using reflection
 * principle
 *
 * Copyright 2020 <University of Pennsylvania
 *
 * Author: Joseph K. Aicher, Jorge Vaquero-Garcia
 */

#ifndef MAJIQINCLUDE_STATS_TNOM_HPP
#define MAJIQINCLUDE_STATS_TNOM_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>
#include <limits>
#include <map>
#include <utility>

#include <majiqinclude/LogMath.hpp>

namespace MajiqInclude {
namespace TNOM {

struct TNOMRecord {
  // counts of two classes. We use neg/pos rather than 1, 2 motivated by
  // analysis of test statistic using reflection principle (random walk in
  // positive, negative direction depending on next label)
  int n_neg;
  int n_pos;
  // tnom score
  int score;
};
inline bool operator<(const TNOMRecord& x, const TNOMRecord& y) {
  return std::tie(x.n_neg, x.n_pos, x.score)
    < std::tie(y.n_neg, y.n_pos, y.score);
}
inline bool operator==(const TNOMRecord& x, const TNOMRecord& y) {
  return std::tie(x.n_neg, x.n_pos, x.score)
    == std::tie(y.n_neg, y.n_pos, y.score);
}
inline bool operator!=(const TNOMRecord& x, const TNOMRecord& y) {
  return !(x == y);
}

class TNOMCache {
 private:
  std::map<TNOMRecord, double> pvalue_;  // cached p-value for given score

  /**
   * (M = length, C = offset)
   * The number of length M sequences of +/- 1 that when added end in C.
   *
   * These can be seen as discrete paths starting at (0, 0) ending in (M, C).
   * Note that n_neg + n_pos = M, n_pos - n_neg = C considering all possible
   * paths.
   *
   * More generally, if a path starts at t and ends at c, we could count using
   * the offset: C = c - t.
   *
   * This computation is used in Lemma 3.3 of Ben-Dor 2002.
   * Gamma(w) for alternating pattern w is defined considering paths that start
   * at t(w) defined appropriately by considering if the path last hits its
   * positive/negative threshold while still being more extreme (TNOM closer to
   * 0) than the actual score.
   *
   * @param length (M): how many items in the sequence (e.g. n_neg + n_pos)
   * @param offset (C): final offset after ending sequence (e.g. n_pos - n_neg)
   */
  static double LogDistinctPaths(int length, int offset) {
    // always do positive offset (symmetric)
    if (offset < 0) {
      offset = -offset;
    }
    // not possible to have offset greater than length and evenness must be same
    if (offset > length || (offset + length) % 2 != 0) {
      return std::numeric_limits<double>::lowest();  // effectively log(0)
    }
    return detail::lchoose((offset + length) / 2, length);
  }

 public:
  /**
   * Compute proportion of sequences with n_neg, n_pos that have TNOM <= score
   *
   * @param n_neg number of "negative" class
   * @param n_pos number of "positive" class
   * @param score observed TNOM we are computing p-value for
   */
  double CalculatePValue(int n_neg, int n_pos, int score) {
    if (n_pos < n_neg) {
      // by symmetry these are identical, but only need to cache one
      std::swap(n_neg, n_pos);  // n_pos >= n_neg
    }
    // use cached result if present
    TNOMRecord key{n_neg, n_pos, score};
    auto lb = pvalue_.lower_bound(key);
    if (lb != pvalue_.end() && lb->first == key) {
      // there was a cached result
      return lb->second;
    } else {
      // no cached result, compute it now
      double result{1};

      /**
       * Ben-Dor 2002 Proposition 3.1 says that a path has TNOM <= score if and
       * only if it hits n_pos - score (A, positive threshold) or score - n_neg
       * (-B, negative threshold) at least once.
       */
      const int A = n_pos - score;
      const int B = n_neg - score;
      /**
       * We want to count the paths that start at (0, 0), end at (M, C)
       * (M = length of path, C = final offset at end of path)
       * that hit A or -B at least once.
       */
      const int M = n_pos + n_neg;
      const int C = n_pos - n_neg;  // non-negative because of swap at beginning

      // so, if A or B are 0, this is all paths.
      if (A <= 0 || B <= 0) {
        result = double{1};
      } else {
        /**
         * Lemma 3.2 tells us how to count how many paths there are as an
         * alternating sum:
         * + first term: count all paths that hit A, all paths that hit B
         * - second term: this double counts paths that hit A and B, so
         *   subtract paths that hit AB, that hit BA
         * + but this doubly subtracted paths that hit ABA or BAB, so add. And
         *   so on.
         * This terminates because we can only bounce between A and B for so
         * long.
         *
         * These terms are represented by Lambda(w) (for alternating pattern
         * between A and B) (Lemma 3.2), which are themselves path counts
         * from offsets t(w) to C, with t(w) defined in Lemma 3.3.
         */

        // accumulate positive and negative terms separately
        double logsum_pos = std::numeric_limits<double>::lowest();
        double logsum_neg = std::numeric_limits<double>::lowest();

        // starting at patterns with length 1:
        // Ni = -l(w) when w ends with B, Pi = l(w) when ends with A (Lemma 3.3)
        int Ni = 2 * B;  // negative offset on starting position, passes -B last
        int Pi = 2 * A;  // positive offset on starting position, passes A last

        bool add = true;
        while (Ni + C <= M || std::abs(Pi - C) <= M) {
          // which term are we updating?
          double& logsum_update = add ? logsum_pos : logsum_neg;
          // offset used is distance between offset and Pi, offset and -Ni
          logsum_update = detail::logadd(
              logsum_update, LogDistinctPaths(M, C - Pi));
          logsum_update = detail::logadd(
              logsum_update, LogDistinctPaths(M, C + Ni));
          // flip sign for next update
          add = !add;
          // update Ni, Pi using old values (of the other, since alternating)
          const int oldNi = Ni;
          const int oldPi = Pi;
          Ni = 2 * B + oldPi;
          Pi = 2 * A + oldNi;
        }  // done accumulating alternating sum

        // normalize positive, negative terms by total number of paths
        const double logZ = LogDistinctPaths(M, C);
        logsum_pos -= logZ;
        logsum_neg -= logZ;

        // compute result
        result = std::exp(logsum_pos);
        if (logsum_neg > std::numeric_limits<double>::lowest()) {
          result -= std::exp(logsum_neg);
        }
      }  // end nontrivial case (A, B != 0)
      pvalue_.insert(lb, std::make_pair(key, result));
      return result;
    }
  }

  TNOMCache() = default;
  TNOMCache(const TNOMCache&) = default;
  TNOMCache(TNOMCache&&) = default;
  TNOMCache& operator=(const TNOMCache&) = default;
  TNOMCache& operator=(TNOMCache&&) = default;
};

// perform TNOM test on provided random-access iterators
template <typename ItX, typename ItSort, typename ItLabels>
inline double Test(
    TNOMCache& tester, ItX x, ItSort sortx, ItLabels labels, int64_t d) {
  // counts for labels on RHS
  int rhs1 = 0;
  int rhs2 = 0;
  for (int64_t idx = 0; idx < d; ++idx) {
    const auto& j = sortx[idx];
    const auto& xj = x[j];
    if (std::isnan(xj)) {
      // missing values sorted to end, so done on this pass
      break;
    }
    // update counts on rhs
    if (labels[j]) {
      ++rhs1;
    } else {
      ++rhs2;
    }
  }  // count instances of label 1 vs label 2, set on rhs
  if (rhs1 == 0 || rhs2 == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const int64_t n = static_cast<int64_t>(rhs1 + rhs2);  // total quantified
  // counts for labels on LHS
  int lhs1 = 0;
  int lhs2 = 0;
  // right now partitioning is nothing on left, everything on right
  int min_score = std::min(rhs1, rhs2);
  // consider all other partitions
  for (int64_t idx = 0; idx < n; ++idx) {
    if (labels[sortx[idx]]) {
      // label is 1, move from right to left
      ++lhs1;
      --rhs1;
    } else {
      // move label 2 from right to left
      ++lhs2;
      --rhs2;
    }
    // update min_score considering this partition
    min_score = std::min(
        min_score, std::min(lhs1, lhs2) + std::min(rhs1, rhs2));
  }  // loop over partitions to calculate TNOM score
  // lhs has counts for both groups now, so:
  return tester.CalculatePValue(lhs1, lhs2, min_score);
}

}  // namespace TNOM
}  // namespace MajiqInclude


#endif  // MAJIQINCLUDE_STATS_TNOM_HPP
