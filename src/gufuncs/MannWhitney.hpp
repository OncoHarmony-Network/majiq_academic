/**
 * MannWhitney.hpp
 *
 * Implementation of Mann-Whitney U (two-sample Wilcoxon) two-sided test for
 * use in MAJIQ.
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Mann-Whitney U is a two-sample test under the null hypothesis that the
 * distributions generating two samples with sizes n1 and n2 are the same.
 * The procedure calculates the ranks of each observation, appropriately
 * handling ties. The sum of the ranks will always be N*(N+1)/2, where N=n1+n2,
 * and the test statistic is derived from the sum of the ranks in either of the
 * samples (one determines the other given that the sum of the ranks is
 * invariant), say R1, R2. The minimum value for Ri (i=1 or 2) is ni*(ni+1)/2,
 * so the test-statistic Ui provides a shift to Ri to make the minimum value 0.
 * Then, U1+U2 is equal to n1*n2, which is the maximum value for either
 * statistic. Note then that there are exactly 2*n1*n2 + 1 possible values that
 * they can take, which makes it amenable to dynamic programming/cached
 * evaluation.
 *
 * Mann and Whitney, in their 1947 paper "On a Test of Whether one of Two
 * Random Variables is Stochastically Larger than the Other," note that, if we
 * can ignore ties, the null hypothesis/distribution can be represented by a
 * uniform distribution over binary sequences of length N with n1 0's and n2
 * 1's. The test-statistic U1 is then the number of times that a 1 preceds a 0.
 * This enables a recurrence relation to be created for the number of sequences
 * with a given value of U. This can be used to compute exact probabilities for
 * realizations of U under the null hypothesis, allowing us to compute p-values.
 * This is not the most efficient way of calculating p-values, but it is
 * straightforward and easy to implement. This is only feasible for small
 * samples, but for larger samples there exists a normal asymptotic
 * approximation that is good enough.
 *
 * The relevant summary statistics from the data are calculated in the
 * class/struct MajiqGufuncs::MannWhitneySummary, which counts the sum
 * of ranks and number of samples in each group with valid quantifications (>=
 * 0), and tie correction term for normal approximation.
 * The pvalue is chosen from exact calculation vs asymptotic calculation (using
 * normal approximation with continuity correction) depending on whether the
 * total sample size is greater than MANNWHITNEY_MAXSIZE.
 * Asymptotic calculation is done as described on Wikipedia (continuity
 * correction as described in SciPy docs/source)
 * Exact calculation uses exact Wilcoxon distribution as described in Mann and
 * Whitney's 1947 paper. Uses dynamic programming to cache calculation of
 * event cardinalities in null distribution for single or cumulative values of
 * test statistic, and for binomial coefficients, with keys/hashes as
 * implemented in MajiqStats::details.
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_MANNWHITNEY_HPP
#define MAJIQGUFUNCS_MANNWHITNEY_HPP

#include <numpy/ndarraytypes.h>

#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include "helpers.hpp"

#include <boost/math/distributions/normal.hpp>


namespace MajiqGufuncs {
namespace MannWhitney {

constexpr int_fast64_t MANNWHITNEY_MAXSIZE = 64;

/**
 * Summary statistics describing two samples for test statistic
 */
struct MannWhitneySummary {
  int_fast64_t n1;  // number of samples with first label
  int_fast64_t n2;  // number of samples with second label
  double U1;  // test statistic
  double tie_correction;  // correction term for ties in asymptotic approximation

  MannWhitneySummary(
      int_fast64_t _n1, int_fast64_t _n2, double _U1, double _tie_correction)
    : n1{_n1}, n2{_n2}, U1{_U1}, tie_correction{_tie_correction} { }
  MannWhitneySummary() : n1{0}, n2{0}, U1{0}, tie_correction{0} { }
  MannWhitneySummary(const MannWhitneySummary&) = default;
  MannWhitneySummary(MannWhitneySummary&&) = default;
  MannWhitneySummary& operator=(const MannWhitneySummary&) = default;
  MannWhitneySummary& operator=(MannWhitneySummary&&) = default;
};

struct MannWhitneyStatRecord {
  int_fast64_t n1;
  int_fast64_t n2;
  int_fast64_t U1;

  MannWhitneyStatRecord(int_fast64_t n1, int_fast64_t n2, int_fast64_t U1)
    : n1(n1), n2(n2), U1(U1) { }
};
inline bool operator<(
    const MannWhitneyStatRecord& x, const MannWhitneyStatRecord& y) {
  return std::tie(x.n1, x.n2, x.U1) < std::tie(y.n1, y.n2, y.U1);
}
inline bool operator==(
    const MannWhitneyStatRecord &x, const MannWhitneyStatRecord &y) {
  return std::tie(x.n1, x.n2, x.U1) == std::tie(y.n1, y.n2, y.U1);
}

struct MannWhitneyChooseRecord {
  int_fast64_t n;
  int_fast64_t k;

  MannWhitneyChooseRecord(int_fast64_t n, int_fast64_t k) : n(n), k(k) { }
};
inline bool operator<(
    const MannWhitneyChooseRecord& x, const MannWhitneyChooseRecord& y) {
  return std::tie(x.n, x.k) < std::tie(y.n, y.k);
}
inline bool operator==(
    const MannWhitneyChooseRecord &x, const MannWhitneyChooseRecord &y) {
  return std::tie(x.n, x.k) == std::tie(y.n, y.k);
}


class MannWhitneyCache {
 private:
  // introduce caching variables and locks on them
  std::map<MannWhitneyStatRecord, int_fast64_t> count_cache_;
  std::map<MannWhitneyChooseRecord, int_fast64_t> choose_cache_;
  std::map<MannWhitneyStatRecord, int_fast64_t> cumcount_cache_;
  // standard normal distribution set by constructor (mean 0, var 1)
  const boost::math::normal_distribution<double> z_dist_;

 public:
  MannWhitneyCache() : z_dist_(0., 1.) { }

  /** number of outcomes in null distribution sample space with test-statistic
   *
   * @param n1: sample size of one of the samples
   * @param n2: sample size of the other samples
   * @param U: integer-valued test statistic
   *
   * @note the null distribution is uniform over the set of binary
   * sequences of length (n1+n2) with n1 0's, and the test statistic
   * U is the number of times a 1 precedes a 0 in a sequence. By
   * symmetry, U is exchangable with `n1*n2 - U`, and n1 is
   * exchangable with n2. We use these symmetries to reduce the
   * number of values we need to save
   *
   * @note should not overflow if n1 + n2 <= 64, but can overflow
   * after that. We should not be computing exact counts for these
   * larger sample sizes
   */
  int_fast64_t CountWithStatistic(
          int_fast64_t n1, int_fast64_t n2, int_fast64_t U
  ) {
    // apply symmetry on n1, n2 --> pick n1 <= n2
    if (n1 > n2) {
      std::swap(n1, n2);
    }
    // apply symmetry on U --> pick smaller value of U.
    U = std::min(U, n1 * n2 - U);

    // handle base cases
    if ((n2 <= 0) || (n1 < 0) || (U < 0) || (n1 == 0 && U > 0)) {
      // n2 <= 0: 0 total samples, sequence length 0
      // n1 < 0: invalid sample size
      // U < 0: invalid value of test statistic
      // n1 == 0 && U > 0: n1 == 0 implies U = 0.
      return 0;
    } else if (n1 == 0 /*&& U == 0*/) {
      // previous conditions (U < 0), (n1 == 0 && U > 0) make
      // n1 == 0 imply that U == 0
      return 1;
    }

    // use cached result if present
    MannWhitneyStatRecord record{n1, n2, U};
    // get iterator to record (insert if doesn't exist)
    auto [result_iter, inserted] = count_cache_.try_emplace(record, 0);
    if (inserted) {
      // result was just created, so compute it
      result_iter->second = CountWithStatistic(n1, n2 - 1, U)
        + CountWithStatistic(n1 - 1, n2, U - n2);
    }
    return result_iter->second;
  }

  /** Cached computation of n choose k using Pascal's Triangle
   *
   * @param n, k arguments to n choose k := n! / (k! (n-k)!)
   *
   * @note should not overflow if n <= 64, but can overflow after
   * that
   * @note values n < 0, k < 0, k > n are chosen to return 0
   * @note uses symmetry relation n choose k == n choose n-k
   */
  int_fast64_t Choose(int_fast64_t n, int_fast64_t k) {
    // invalid value of n
    if (n < 0) {
      return 0;
    }
    // symmetry on k, so keep k <= n / 2
    if (k >= n / 2) {
      k = n - k;
    }
    // invalid value of k
    if (k < 0) {
      return 0;
    } else if (k == 0) {
      // base case
      return 1;
    }

    // use cached result if present
    MannWhitneyChooseRecord record{n, k};
    // get iterator to record (insert if doesn't exist)
    auto [result_iter, inserted] = choose_cache_.try_emplace(record, 0);
    if (inserted) {
      // result was just created, so compute it
      result_iter->second = Choose(n - 1, k) + Choose(n - 1, k - 1);
    }
    return result_iter->second;
  }

  /** Number of outcomes in null distribution with statistic less
   * than or equal to U
   *
   * Number of outcomes in null distribution with test statistic less
   * than or equal to U. Note that this is inclusive of U.
   *
   * @param n1: sample size of one of the samples
   * @param n2: sample size of the other samples
   * @param U: integer-valued test statistic
   *
   * @note cum_count_left(U) + cum_count_right(U+1) = choose(n1+n2, n1)
   * will be used for large U
   */
  int_fast64_t CumulativeCountFromLeft(
      int_fast64_t n1, int_fast64_t n2, int_fast64_t U) {
    // invalid values of n1, n2 return 0
    if (n1 < 0 || n2 < 0) {
      return 0;
    }
    // count from right if that would be less values to add
    if (U > n1 * n2 / 2) {
       // cum_count_left(U) + cum_count_right(U+1) = choose(n1+n2, n1)
       // return Choose(n1 + n2, n1) - CumulativeCountFromRight(n1, n2, U + 1);
       return Choose(n1 + n2, n1)
         - CumulativeCountFromLeft(n1, n2, n1 * n2 - (U + 1));
    }
    // invalid value of U return 0
    if (U < 0) {
      return 0;
    }

    // use cached result if present
    MannWhitneyStatRecord record{n1, n2, U};
    // get iterator to record (insert if doesn't exist)
    auto [result_iter, inserted] = cumcount_cache_.try_emplace(record, 0);
    if (inserted) {
      // result was just created, so compute it
      result_iter->second = CumulativeCountFromLeft(n1, n2, U - 1)
        + CountWithStatistic(n1, n2, U);
    }
    return result_iter->second;
  }

  /** Number of outcomes in null distribution with statistic greater
   * than or equal to U
   *
   * Number of outcomes in null distribution with test statistic greater
   * than or equal to U. Note that this is inclusive of U.
   *
   * @param n1: sample size of one of the samples
   * @param n2: sample size of the other samples
   * @param U: integer-valued test statistic
   *
   * @note cumcount_right(U) = cumcount_left(n1*n2 - U)
   */
  int_fast64_t CumulativeCountFromRight(
      int_fast64_t n1, int_fast64_t n2, int_fast64_t U) {
    return CumulativeCountFromLeft(n1, n2, n1 * n2 - U);
  }

  /** Calculate two-sided exact p-value
   *
   * @note return exact value
   */
  double ExactPValue(int_fast64_t n1, int_fast64_t n2, double U1) {
    // maximum possible value of U
    const int_fast64_t max_U = n1 * n2;
    // mean value of U under null distribution
    const double mean_U_null = static_cast<double>(max_U) / 2.;
    // always compute pvalue from left -- flip statistic around
    if (U1 > mean_U_null) {
      return ExactPValue(n1, n2, max_U - U1);
    }
    // cumulative probability of U1 under null distribution
    const int_fast64_t p_numerator = CumulativeCountFromLeft(
        n1, n2, static_cast<int_fast64_t>(std::floor(U1)));
    const int_fast64_t p_denominator = Choose(n1 + n2, n1);
    const double prob_U1 = static_cast<double>(p_numerator) / p_denominator;
    // We want two-sided p-value, so multiply by 2
    return std::min(2 * prob_U1, 1.);
  }

  /** Calculate two-sided asymptotic p-value using normal approximation
   *
   * @note See <https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction>
   * for details, but in particular:
   * mean_U = n1*n2 / 2
   * std_U = sqrt(n1*n2/12 * ((n1+n2+1) - tie_correction))
   *
   * @note Use continuity correction -- 1/2 towards mean
   */
  double AsymptoticPValue(
      int_fast64_t n1, int_fast64_t n2, double U1, double tie_correction) {
    // shared intermediates
    const int_fast64_t n1n2 = n1 * n2;
    const int_fast64_t n = n1 + n2;
    // mean of null distribution
    const double mean_U = static_cast<double>(n1n2) / 2.;
    // standard deviation of null distribution (min value is 0)
    const double var_U = std::max(n1n2 * ((n + 1) - tie_correction) / 12., 0.);
    const double std_U = std::sqrt(var_U);
    // z = (U1 - mean_U) / (std_U + eps)
    // continuity correction of 0.5 in numerator towards 0
    // use negative numerator so we can just use CDF
    const double z_numerator = -std::max(std::fabs(U1 - mean_U) - 0.5, 0.);
    // calculate p-value using standard normal distribution
    const double pval = 2. * boost::math::cdf(z_dist_, z_numerator / std_U);
    return std::max(0., std::min(pval, 1.));
  }

  double CalculatePValue(const MannWhitneySummary& summary) {
    // determine whether we are doing asymptotic p value or not
    if (summary.n1 == 0 || summary.n2 == 0) {
      return std::numeric_limits<double>::quiet_NaN();
    } else if (summary.n1 + summary.n2 <= MANNWHITNEY_MAXSIZE) {
      return ExactPValue(summary.n1, summary.n2, summary.U1);
    } else {
      return AsymptoticPValue(
          summary.n1, summary.n2, summary.U1, summary.tie_correction);
    }
  }
};

// perform MannWhitney test on provided random-access iterators
template <typename ItX, typename ItSort, typename ItLabels>
inline double Test(
    MannWhitneyCache& tester, ItX x, ItSort sortx, ItLabels labels, npy_intp d) {
  static_assert(
      std::is_same_v<std::random_access_iterator_tag,
      typename std::iterator_traits<ItX>::iterator_category>,
      "MannWhitney::Test requires x to be random-access iterator");
  static_assert(
      std::is_same_v<std::random_access_iterator_tag,
      typename std::iterator_traits<ItSort>::iterator_category>,
      "MannWhitney::Test requires sortx to be random-access iterator");
  static_assert(
      std::is_same_v<std::random_access_iterator_tag,
      typename std::iterator_traits<ItLabels>::iterator_category>,
      "MannWhitney::Test requires labels to be random-access iterator");

  // get MannWhitneySummary for the row
  MannWhitneySummary summary;
  {
    // NOTE: n1 ~ label = True, n2 ~ label = False
    // running accumulator of the sum of ranks for group 1
    double R1 = 0.;
    // running accumulator for numerator in tie correction term
    int_fast64_t tie_numerator = 0;

    // loop through data/labels, noting ties, ignoring missing data
    npy_intp idx = 0;
    int last_rank = 0;
    while (idx < d) {
      const auto& val = x[sortx[idx]];
      if (npy_isnan(val)) {
        // missing values should be sorted to end, so we can break
        break;
      }
      // iterate over data to count number of samples with current value
      int_fast64_t val_n1 = 0;
      int_fast64_t val_n2 = 0;
      for (; idx < d; ++idx) {
        const auto& j = sortx[idx];
        if (val != x[j]) {
          // no longer tied
          break;
        }
        if (labels[j]) {
          ++val_n1;
        } else {
          ++val_n2;
        }
      }
      // accumulate n1, n2, R1, last_rank, tie_numerator
      summary.n1 += val_n1;
      summary.n2 += val_n2;
      const int_fast64_t val_n = val_n1 + val_n2;  // total with value
      if (val_n == 1) {
        // no ties, straightforward to update last rank, R1
        ++last_rank;
        if (val_n1 == 1) {
          R1 += last_rank;
        }
      } else {
        // we have ties
        // last_rank + 0.5 * (val_n + 1) is the value we use for tied rank
        R1 += val_n1 * (last_rank + 0.5 * (val_n + 1));
        last_rank += val_n;
        // accumulate for ties
        tie_numerator = val_n * (val_n * val_n - 1);
      }
    }  // while (idx < _x.shape(1))

    if (summary.n1 > 0 && summary.n2 > 0) {
      // get U1(R1, summary.n1)
      summary.U1 = R1 - (summary.n1 * (summary.n1 + 1)) / 2;
      // calculate tie correction term
      summary.tie_correction = static_cast<double>(tie_numerator)
        / ((summary.n1 + summary.n2) * (summary.n1 + summary.n2 - 1));
    }
  }
  return tester.CalculatePValue(summary);
}

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp stride_x = *steps++;
  const npy_intp stride_sortx = *steps++;
  const npy_intp stride_labels = *steps++;
  const npy_intp stride_out = *steps++;
  // core dimensions
  const npy_intp dim_core = dimensions[0];
  // inner strides
  const npy_intp inner_stride_x = steps[0];
  const npy_intp inner_stride_sortx = steps[1];
  const npy_intp inner_stride_labels = steps[2];
  // pointers to data
  char* x_ptr = args[0];
  char* sortx_ptr = args[1];
  char* labels_ptr = args[2];
  char* out_ptr = args[3];

  // create cache of previous tests for MannWhitney
  MannWhitneyCache tester;

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i,
      x_ptr += stride_x, sortx_ptr += stride_sortx,
      labels_ptr += stride_labels, out_ptr += stride_out) {
    detail::get_value<RealT>(out_ptr) = Test(
        tester,
        detail::CoreIt<RealT>::begin(x_ptr, inner_stride_x),
        detail::CoreIt<npy_intp>::begin(sortx_ptr, inner_stride_sortx),
        detail::CoreIt<bool>::begin(labels_ptr, inner_stride_labels),
        dim_core);
  }
  return;
}

static char name[] = "mannwhitneyu";
static char signature[] = "(n),(n),(n)->()";
static char doc[] = R"pbdoc(
Compute p-values for Mann-Whitney U test on input data

Compute p-values for Mann-Whitney U test on input data, using a two-sided
alternative hypothesis and omitting nan values.

Parameters
----------
x: array[float]
    test over observations in last axis
sortx: array[int]
    Previously computed values of np.argsort(x, axis=-1) that will not be
    checked
labels: array[bool]
    test over labels in last axis

Returns
-------
array[float]
    broadcast p-values for observations/labels. Invalid tests are nan
)pbdoc";
constexpr int ntypes = 2;
constexpr int nin = 3;
constexpr int nout = 1;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>)
};
static char types[
  ntypes * (nin + nout)
] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_INTP, NPY_BOOL, NPY_DOUBLE,
  // for use with npy_double func
  NPY_DOUBLE, NPY_INTP, NPY_BOOL, NPY_DOUBLE,
};

}  // namespace MannWhitney
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_MANNWHITNEY_HPP
