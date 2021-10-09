/**
 * InfoScore.hpp
 *
 * Statistics relative to InfoScore (see ScoreGenes)
 *
 * Copyright 2020 <University of Pennsylvania
 */

#ifndef MAJIQINCLUDE_STATS_INFOSCORE_HPP
#define MAJIQINCLUDE_STATS_INFOSCORE_HPP

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include <majiqinclude/LogMath.hpp>

namespace MajiqInclude {
namespace InfoScore {

struct InfoScoreRecord {
  // counts of two classes. We use neg/pos motivated by analysis of sequences
  // as discrete walks with labels {-1, +1}
  int n_neg;
  int n_pos;
  // score associated
  double score;
};
inline bool operator<(const InfoScoreRecord& x, const InfoScoreRecord& y) {
  return std::tie(x.n_neg, x.n_pos, x.score)
    < std::tie(y.n_neg, y.n_pos, y.score);
}
inline bool operator==(const InfoScoreRecord& x, const InfoScoreRecord& y) {
  return std::tie(x.n_neg, x.n_pos, x.score)
    == std::tie(y.n_neg, y.n_pos, y.score);
}
inline bool operator!=(const InfoScoreRecord& x, const InfoScoreRecord& y) {
  return !(x == y);
}

// empirical entropy for specified counts of observations multiplied by the
// total number of observations
template <std::size_t num_classes = 2>
inline double InfoScoreEntropy(std::array<int, num_classes> class_ct) {
  const int total_ct = std::accumulate(class_ct.begin(), class_ct.end(), int{});
  return std::accumulate(
      class_ct.begin(), class_ct.end(), double{0},
      [total_ct](double sum, int ct) -> double {
      if (ct == 0) {
        return sum;
      } else {
        return sum - ct * std::log(static_cast<double>(ct) / total_ct);
      } });
}


/**
 * Compute entropy for partition with k labels on left with current offset
 * Entropy(lhs) + Entropy(rhs)
 * Note that Entropy is actually weighted by the number of observations too
 */
inline double InfoScorePathScore(int n_neg, int n_pos, int k, int offset) {
  std::array<int, 2> lhs;
  std::array<int, 2> rhs;
  lhs[0] = (k - offset) / 2;  // negative examples
  lhs[1] = k - lhs[0];  // positive examples
  rhs[0] = n_neg - lhs[0];
  rhs[1] = n_pos - lhs[1];

  return InfoScoreEntropy(lhs) + InfoScoreEntropy(rhs);
}


class InfoScoreCache {
 private:
  std::map<InfoScoreRecord, double> pvalue_;

 public:
  /**
   * Calculate p-value for sequence with n_neg, n_pos negative/positive labels
   * and given info score
   *
   * @param n_neg: number of negative labels
   * @param n_pos: number of positive labels
   * @param score: infoscore, the minimum path score over all possible partitions
   *
   * @return p-value (Pr{score(X) <= observed score : sequences(n_neg, n_pos)})
   */
  double CalculatePValue(int n_neg, int n_pos, double score) {
    if (n_pos < n_neg) {
      // by symmetry these are identical, but only need to cache one
      std::swap(n_neg, n_pos);
    }
    const InfoScoreRecord key{n_neg, n_pos, score};
    // search for key in pvalue_ records
    auto lb = pvalue_.lower_bound(key);
    if (lb != pvalue_.end() && lb->first == key) {
      // use cached result
      return lb->second;
    } else {
      // not cached, compute it

      // length of path, final offset
      const int N = n_pos + n_neg;
      const int O = n_pos - n_neg;

      // Dynamic programming O(N^2) compute
      // naively O(N^2) grid, but use two O(N) vectors instead
      // Note: these counts are in *logspace*
      std::vector<double> OldCounts(
          1 + N, std::numeric_limits<double>::lowest());
      std::vector<double> NewCounts(
          1 + N, std::numeric_limits<double>::lowest());
      // Indexes 0, ..., L map to offsets -n_neg, ..., n_pos
      // Base case (length 0 path): there is exactly one path (point) that is
      // valid, which is at (0, 0)
      OldCounts[0 + n_neg] = 0;

      // what are the offsets of our last minimum/maximum value within
      // acceptable region? (paths that have smaller info score)
      // at k = 0, it's only possible to have offset 0.
      int last_max = 0;
      int last_min = 0;

      // an alternative count of bad paths that achieve the score
      double bad_paths = std::numeric_limits<double>::lowest();

      // loop over offsets
      for (int k = 1; k <= N; ++k) {
        // highest possible offset on this new step
        int M = std::min(1 + last_max, std::min(n_pos, O + N - k));
        // lowest possible offset on this new step
        int m = std::max(last_min - 1, std::max(-n_neg, O - (N - k)));

        // initialize new_min, new_max to invalid range (we hope that there are
        // paths afterward that are valid
        int new_min = n_pos;
        int new_max = -n_neg;

        double S = 0;
        // loop over possible offsets
        for (int i = m; i <= M; ++i) {
          // initialize new counts to log(0)
          NewCounts[i + n_neg] = std::numeric_limits<double>::lowest();
          if (i + k % 2 == 1) {
            // not possible to have odd/even or even/odd position/offset
            continue;
          }
          // how many valid paths can end at this position/offset (k, i)?
          double C = std::numeric_limits<double>::lowest();
          // consider i - 1, i + 1 from previous counts
          if (i - 1 >= last_min && i - 1 <= last_max) {
            C = detail::logadd(C, OldCounts[(i - 1) + n_neg]);
          }
          if (i + 1 >= last_min && i + 1 <= last_max) {
            C = detail::logadd(C, OldCounts[(i + 1) + n_neg]);
          }
          // how many valid paths (up until this point) could get to (k, i)?
          if (C > std::numeric_limits<double>::lowest()) {
            // what's the score at this point?
            S = InfoScorePathScore(n_neg, n_pos, k, i);
            // so is this within the allowable region
            if (S > score) {
              // valid path within allowable region
              NewCounts[i + n_neg] = C;  // so we update new counts
              // update new_min, new_max
              if (i < new_min) {
                new_min = i;
              }
              if (i > new_max) {
                new_max = i;
              }
            } else {  // S <= score
              // we hit the boundary of allowable region. Every path from here
              // to the end is bad
              // How many paths start at 0, stop here, and end at (N, O)?
              // multiply num paths to here and num paths from here to end
              // (add in log space)
              const int dx = N - k;
              const int dy = O - i;
              const int p = (dx + dy) / 2;
              C += detail::lchoose(p, dx);
              // update count of bad paths
              bad_paths = detail::logadd(bad_paths, C);
            }  // were the paths valid or not?
          }  // any paths (up until this point) that could get to (k, i)?
        }  // loop over possible offsets i at fixed position k?

        if (new_min > new_max) {
          // there are no valid paths (all paths have <= score)
          double result{1};
          pvalue_.insert(lb, std::make_pair(key, result));
          return result;
        }

        // update old values using new values in preparation for next offset
        last_min = new_min;
        last_max = new_max;
        std::swap(NewCounts, OldCounts);
      }  // iterate over positions until reach end of sequence
      auto result = std::exp(bad_paths - detail::lchoose(n_neg, N));
      pvalue_.insert(lb, std::make_pair(key, result));
      return result;
    }
  }
};

// perform InfoScore test on provided random-access iterators
template <typename ItX, typename ItSort, typename ItLabels>
inline double Test(
    InfoScoreCache& tester, ItX x, ItSort sortx, ItLabels labels, int64_t d) {
  // first pass: accumulate count of negative/positive examples
  int n_neg = 0;
  int n_pos = 0;  // label == true
  for (int64_t idx = 0; idx < d; ++idx) {
    const auto& j = sortx[idx];
    const auto& xj = x[j];
    if (std::isnan(xj)) {
      // missing values are sorted to end, so we are done on this pass
      break;
    }
    // update counts on rhs
    if (labels[j]) {
      ++n_pos;
    } else {
      ++n_neg;
    }
  }  // count instances of label 1 vs label 2
  if (n_pos == 0 || n_neg == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const int64_t n = static_cast<int64_t>(n_neg + n_pos);
  double min_score = std::numeric_limits<double>::max();
  int offset = 0;  // initial offset
  for (int64_t k = 0; k < n; ++k) {
    offset += labels[sortx[k]] ? 1 : -1;
    min_score = std::min(
        min_score,
        InfoScorePathScore(n_neg, n_pos, 1 + k, offset));
  }  // loop over offsets to compute minimum score
  // lhs has counts for both groups now, so:
  return tester.CalculatePValue(n_neg, n_pos, min_score);
}

}  // namespace InfoScore
}  // namespace MajiqInclude

#endif  // MAJIQINCLUDE_STATS_INFOSCORE_HPP
