// MannWhitney.hpp
//
// Implementation of Mann-Whitney U (two-sample Wilcoxon) two-sided test for
// use in MAJIQ.
//
// Mann-Whitney U is a two-sample test under the null hypothesis that the
// distributions generating two samples with sizes n1 and n2 are the same.
// The procedure calculates the ranks of each observation, appropriately
// handling ties. The sum of the ranks will always be N*(N+1)/2, where N=n1+n2,
// and the test statistic is derived from the sum of the ranks in either of the
// samples (one determines the other given that the sum of the ranks is
// invariant), say R1, R2. The minimum value for Ri (i=1 or 2) is ni*(ni+1)/2,
// so the test-statistic Ui provides a shift to Ri to make the minimum value 0.
// Then, U1+U2 is equal to n1*n2, which is the maximum value for either
// statistic. Note then that there are exactly 2*n1*n2 + 1 possible values that
// they can take, which makes it amenable to dynamic programming/cached
// evaluation.
//
// Mann and Whitney, in their 1947 paper "On a Test of Whether one of Two
// Random Variables is Stochastically Larger than the Other," note that, if we
// can ignore ties, the null hypothesis/distribution can be represented by a
// uniform distribution over binary sequences of length N with n1 0's and n2
// 1's. The test-statistic U1 is then the number of times that a 1 preceds a 0.
// This enables a recurrence relation to be created for the number of sequences
// with a given value of U. This can be used to compute exact probabilities for
// realizations of U under the null hypothesis, allowing us to compute p-values.
// This is not the most efficient way of calculating p-values, but it is
// straightforward and easy to implement. This is only feasible for small
// samples, but for larger samples there exists a normal asymptotic
// approximation that is good enough.
//
// This header implements MajiqStats::MannWhitney, which provides the
// Calc_pval() interface from MajiqStats::TestStat.
// The relevant summary statistics from the data are calculated in the
// class/struct MajiqStats::details::MannWhitneySummary, which counts the sum
// of ranks and number of samples in each group with valid quantifications (>=
// 0), and tie correction term for normal approximation.
// The pvalue is chosen from exact calculation vs asymptotic calculation (using
// normal approximation with continuity correction) depending on whether the
// total sample size is greater than MANNWHITNEY_MAXSIZE.
// Asymptotic calculation is done as described on Wikipedia (continuity
// correction as described in SciPy docs/source)
// Exact calculation uses exact Wilcoxon distribution as described in Mann and
// Whitney's 1947 paper. Uses dynamic programming to cache calculation of
// event cardinalities in null distribution for single or cumulative values of
// test statistic, and for binomial coefficients, with keys/hashes as
// implemented in MajiqStats::details.
//
// Author: Joseph K Aicher

#ifndef MANNWHITNEY_H
#define MANNWHITNEY_H
// includes
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <omp.h>
#include "testStats.hpp"

// maximum total sample size (n1 + n2) used for exact p-value calculation
// the maximal null sample space comes when n1 = n2 = MAXSIZE / 2. This must
// have less than 2 ** 63 sequences, which allows us to do calculations with
// exact integer arithmetic using signed 64-bit integers.
// This can be set to a lower value to use asymptotic pvalues earlier for
// faster speed (in exchange for loss of accuracy), but it cannot be safely
// increased above 66 without loss of precision
#define MANNWHITNEY_MAXSIZE 64

namespace MajiqStats {
    namespace details {
        /**
         * Summary statistics describing two samples for test statistic
         */
        struct MannWhitneySummary {
            int_fast64_t n1;  // number of samples with first label
            int_fast64_t n2;  // number of samples with second label
            double U1;  // test statistic
            double tie_correction;  // correction term for ties in asymptotic approximation

            /** Process data and labels to produce summary of data for testing
             *
             * @param data vector of observations in sorted order
             * @param labels vector of binary labels corresponding to data
             *
             * @note label 0 --> n1, label 1 --> n2
             */
            MannWhitneySummary(
                    std::vector<float>& data,
                    std::vector<int>& labels
            ) {
                // running count of n1, n2 as we loop through data/labels
                int_fast64_t acc_n1 = 0;
                int_fast64_t acc_n2 = 0;
                // running accumulator of the sum of ranks for group 1
                double R1 = 0.;
                // running accumulator for numerator in tie correction term
                int_fast64_t tie_numerator = 0;

                // loop through data/labels, noting ties, ignoring missing data
                const unsigned int data_size = data.size();  // NOTE: not n1 + n2 because missing data
                unsigned int idx = 0;  // current index of array
                unsigned int last_rank = 0;  // 1 less than next rank (ranks start at 1)
                while (idx < data_size) {
                    // current value -- check if missing, count ties
                    double x = data[idx];
                    // handle missing observation
                    if (x < 0) {
                        ++idx;
                        continue;
                    }
                    // iterate over data to count number of samples with current value
                    int_fast64_t x_n1 = 0;  // number observations in group 1 with value x
                    int_fast64_t x_n2 = 0;  // number observations in group 2 with value x
                    for (; idx < data_size && data[idx] == x; ++idx) {
                        switch(labels[idx]) {
                            case 0:
                                ++x_n1;
                                break;
                            case 1:
                                ++x_n2;
                                break;
                            default:
                                std::cerr << "Invalid label passed to MannWhitney Summary\n";
                                break;
                        }
                    }
                    // acumulate n1, n2, R1, last_rank, tie_numerator
                    acc_n1 += x_n1;
                    acc_n2 += x_n2;
                    const int_fast64_t x_n = x_n1 + x_n2;  // total with value x
                    if (x_n == 1) {
                        // no ties, straightforward to update last_rank, R1
                        ++last_rank;
                        if (x_n1 == 1) {
                            R1 += last_rank;
                        }
                    } else {
                        // we have ties
                        // last_rank + 0.5 * (x_n + 1) is the value we use for tied rank
                        R1 += x_n1 * (last_rank + 0.5 * (x_n + 1));
                        last_rank += x_n;
                        // accumulate for ties
                        tie_numerator = x_n * (x_n * x_n - 1);
                    }
                }
                // so set final values
                n1 = acc_n1;
                n2 = acc_n2;
                U1 = R1 - (n1 * (n1 + 1)) / 2;
                tie_correction = static_cast<double>(tie_numerator) / ((n1 + n2) * (n1 + n2 - 1));
            }
        };

        struct MannWhitneyStatRecord {
            int_fast64_t n1;
            int_fast64_t n2;
            int_fast64_t U1;

            MannWhitneyStatRecord(int_fast64_t n1, int_fast64_t n2, int_fast64_t U1)
                    : n1(n1), n2(n2), U1(U1) {
                // nothing to do
            }

            bool operator<(const MannWhitneyStatRecord& other) const {
                return n1 < other.n1
                    || (n1 == other.n1 && n2 < other.n2)
                    || (n1 == other.n1 && n2 == other.n2 && U1 < other.U1);
            }

            bool operator==(const MannWhitneyStatRecord &other) const {
                return n1 == other.n1 && n2 == other.n2 && U1 == other.U1;
            }
        };

        struct MannWhitneyStatHash {
            std::size_t operator()(const MannWhitneyStatRecord &record) const {
                size_t res = 0;  // initialize value that will hold hash
                hash_combine(res, record.n1);
                hash_combine(res, record.n2);
                hash_combine(res, record.U1);
                return res;
            }
        };

        struct MannWhitneyChooseRecord {
            int_fast64_t n;
            int_fast64_t k;

            MannWhitneyChooseRecord(int_fast64_t n, int_fast64_t k) : n(n), k(k) {
                // nothing to do
            }

            bool operator<(const MannWhitneyChooseRecord& other) const {
                return n < other.n || (n == other.n && k < other.k);
            }

            bool operator==(const MannWhitneyChooseRecord &other) const {
                return n == other.n && k == other.k;
            }
        };

        struct MannWhitneyChooseHash {
            std::size_t operator()(const MannWhitneyChooseRecord &record) const {
                size_t res = 0;  // initialize value that will hold hash
                hash_combine(res, record.n);
                hash_combine(res, record.k);
                return res;
            }
        };
    }


    class MannWhitney : public MajiqStats::TestStat {
        private:
            // introduce caching variables and locks on them
            std::unordered_map<
                details::MannWhitneyStatRecord, int_fast64_t, details::MannWhitneyStatHash
                > count_cache_;
            std::unordered_map<
                details::MannWhitneyChooseRecord, int_fast64_t, details::MannWhitneyChooseHash
                > choose_cache_;
            std::unordered_map<
                details::MannWhitneyStatRecord, int_fast64_t, details::MannWhitneyStatHash
                > cumcount_cache_;
            omp_lock_t lock_count_cache_;
            omp_lock_t lock_choose_cache_;
            omp_lock_t lock_cumcount_cache_;

        public:
            MannWhitney() {
                // set up locks
                omp_init_lock(&lock_count_cache_);
                omp_init_lock(&lock_choose_cache_);
                omp_init_lock(&lock_cumcount_cache_);
            }

            ~MannWhitney() {
                // remove locks
                omp_destroy_lock(&lock_count_cache_);
                omp_destroy_lock(&lock_choose_cache_);
                omp_destroy_lock(&lock_cumcount_cache_);
            }

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

                // use cached result if present, handle locks for multithreading
                details::MannWhitneyStatRecord record(n1, n2, U);
                int_fast64_t result;
                omp_set_lock(&lock_count_cache_);
                if (count_cache_.count(record) > 0) {  // have result already
                    result = count_cache_[record];
                    omp_unset_lock(&lock_count_cache_);
                } else {  // compute from scratch
                    // unset lock while computing result
                    omp_unset_lock(&lock_count_cache_);

                    // compute result
                    result = CountWithStatistic(n1, n2 - 1, U)
                        + CountWithStatistic(n1 - 1, n2, U - n2);

                    // insert result into cache, handling locks
                    omp_set_lock(&lock_count_cache_);
                    count_cache_[record] = result;
                    omp_unset_lock(&lock_count_cache_);
                }

                // return result
                return result;
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

                // use cached result if present, handle locks for multithreading
                details::MannWhitneyChooseRecord record(n, k);
                int_fast64_t result;
                omp_set_lock(&lock_choose_cache_);
                if (choose_cache_.count(record) > 0) {  // have result already
                    result = choose_cache_[record];
                    omp_unset_lock(&lock_choose_cache_);
                } else {  // compute from scratch
                    // unset lock while computing result
                    omp_unset_lock(&lock_choose_cache_);

                    // compute result using preceding terms
                    result = Choose(n - 1, k) + Choose(n - 1, k - 1);

                    // insert result into cache, handling locks
                    omp_set_lock(&lock_choose_cache_);
                    choose_cache_[record] = result;
                    omp_unset_lock(&lock_choose_cache_);
                }

                // return result
                return result;
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
                    int_fast64_t n1, int_fast64_t n2, int_fast64_t U
            ) {
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

                // use cached result if present, handle locks for multithreading
                details::MannWhitneyStatRecord record(n1, n2, U);
                int_fast64_t result;
                omp_set_lock(&lock_cumcount_cache_);
                if (cumcount_cache_.count(record) > 0) {  // have result already
                    result = cumcount_cache_[record];
                    omp_unset_lock(&lock_cumcount_cache_);
                } else {  // compute from scratch
                    // unset lock while computing result
                    omp_unset_lock(&lock_cumcount_cache_);

                    // compute result, use preceding cumulative count and current count
                    result = CumulativeCountFromLeft(n1, n2, U - 1)
                        + CountWithStatistic(n1, n2, U);

                    // insert result into cache, handling locks
                    omp_set_lock(&lock_cumcount_cache_);
                    cumcount_cache_[record] = result;
                    omp_unset_lock(&lock_cumcount_cache_);
                }

                // return result
                return result;
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
                    int_fast64_t n1, int_fast64_t n2, int_fast64_t U
            ) {
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
                const int_fast64_t p_numerator =
                    CumulativeCountFromLeft(n1, n2, static_cast<int_fast64_t>(std::floor(U1)));
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
            double AsymptoticPValue(int_fast64_t n1, int_fast64_t n2, double U1, double tie_correction) {
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
                // calculate p-value using normal distribution
                // TODO use nicer implementation of 2 * standard normal CDF
                const double pval = 1.
                    + std::erf(
                            z_numerator
                            / (M_SQRT2 * std_U + std::numeric_limits<double>::epsilon())
                    );
                return std::max(0., std::min(pval, 1.));
            }

            /** Calculate two-sided Mann-Whitney p-value for specified labels
             *
             * @param data vector of observations, in sorted order
             * @param labels vector of binary values (0 or 1) indicating class
             * label, corresponding to values of data
             * @param score pointer to float that will hold the test statistic
             * (for class 0). This is currently reserved for use by TNOM
             * statistic and unused by MannWhitney test
             *
             * If the sample size is greater than MANNWHITNEY_MAXSIZE, use
             * normal approximation
             *
             */
            double Calc_pval(
                    std::vector<float>& data,
                    std::vector<int>& labels,
                    float* score
            ) {
                // get summary statistic from data
                details::MannWhitneySummary summary(data, labels);
                // // set score to value of U1
                // *score = static_cast<float>(summary.U1);
                // determine whether we are doing asymptotic p value or not
                if (summary.n1 + summary.n2 <= MANNWHITNEY_MAXSIZE) {
                    return ExactPValue(summary.n1, summary.n2, summary.U1);
                } else {
                    return AsymptoticPValue(summary.n1, summary.n2, summary.U1,
                            summary.tie_correction);
                }
            }
    };
}

#endif
