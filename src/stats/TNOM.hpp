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

#ifndef MAJIQSTATS_TNOM_HPP
#define MAJIQSTATS_TNOM_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <algorithm>
#include <limits>
#include <map>
#include <utility>

#include "LogMath.hpp"

namespace MajiqStats {

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
    auto [result_iter, inserted]
      = pvalue_.try_emplace({n_neg, n_pos, score}, 1);
    if (inserted) {
      // result was just created, so compute it

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
        result_iter->second = 1.;
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
        double result = std::exp(logsum_pos);
        if (logsum_neg > std::numeric_limits<double>::lowest()) {
          result -= std::exp(logsum_neg);
        }
        result_iter->second = result;
      }  // end nontrivial case (A, B != 0)
    }  // end computing result when not cached
    return result_iter->second;
  }

  TNOMCache() = default;
  TNOMCache(const TNOMCache&) = default;
  TNOMCache(TNOMCache&&) = default;
  TNOMCache& operator=(const TNOMCache&) = default;
  TNOMCache& operator=(TNOMCache&&) = default;
};

/**
 * Compute p-values for TNOM test on input data
 *
 * @param x: 2D array, test for each row over observations in columns
 * @param sortx: 2D array, indexes that sort x per row
 * @param labels: 2D array, class labels for each observation
 *
 * @return p-value for test statistic for each row (1D array)
 */
template <typename RealT>
pybind11::array_t<double> TNOM(
    pybind11::array_t<RealT> x,
    pybind11::array_t<ssize_t> sortx,
    pybind11::array_t<bool> labels) {
  // check for correct number of dimensions
  if (x.ndim() != 2) {
    throw std::runtime_error("x is not 2-dimensional");
  } else if (sortx.ndim() != 2) {
    throw std::runtime_error("sortx is not 2-dimensional");
  } else if (labels.ndim() != 2) {
    throw std::runtime_error("labels is not 2-dimensional");
  }
  // check that dimensions match
  if (x.shape(0) != sortx.shape(0) || x.shape(1) != sortx.shape(1)) {
    throw std::runtime_error("x.shape does not match sortx.shape");
  } else if (x.shape(0) != labels.shape(0) || x.shape(1) != labels.shape(1)) {
    throw std::runtime_error("x.shape does not match labels.shape");
  }

  // create output array, 1D with length equal to rows of x
  pybind11::array_t<double> result(x.shape(0));
  // unchecked access to the array values
  auto _x = x.template unchecked<2>();
  auto _sortx = sortx.template unchecked<2>();
  auto _labels = labels.template unchecked<2>();
  auto _result = result.template mutable_unchecked<1>();

  // create object that caches previously computed statistics
  TNOMCache tester;

  // calculate statistics per row of x
  for (pybind11::ssize_t i = 0; i < _x.shape(0); ++i) {
    // counts for labels on RHS
    int rhs1 = 0;
    int rhs2 = 0;
    for (pybind11::ssize_t idx = 0; idx < _x.shape(1); ++idx) {
      const pybind11::ssize_t j = _sortx(i, idx);
      if (std::isnan(_x(i, j))) {
        // missing values are sorted to end, so we are done on this pass
        break;
      }
      // update counts on rhs
      if (_labels(i, j)) {
        ++rhs1;
      } else {
        ++rhs2;
      }
    }  // count instances of label 1 vs label 2, set on rhs
    // if either group has none quantified, no p-value
    if (rhs1 == 0 || rhs2 == 0) {
      _result(i) = std::numeric_limits<double>::quiet_NaN();
      continue;  // compute for next test
    }
    const pybind11::ssize_t n = static_cast<pybind11::ssize_t>(rhs1 + rhs2);
    // counts for labels on LHS
    int lhs1 = 0;
    int lhs2 = 0;
    // right now partitioning is nothing on left, everything on right
    int min_score = std::min(rhs1, rhs2);
    // consider all other partitions
    for (pybind11::ssize_t idx = 0; idx < n; ++idx) {
      if (_labels(i, _sortx(i, idx))) {
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
    _result(i) = tester.CalculatePValue(lhs1, lhs2, min_score);
  }

  // return result
  return result;
}

}  // namespace MajiqStats


#endif  // MAJIQSTATS_TNOM_HPP