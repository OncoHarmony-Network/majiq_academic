// tTest.hpp
//
// Implementation of Welch's t-test (two-sample t-test assuming unequal
// variances) for use in MAJIQ.
//
// Authors: Joseph K Aicher, Jorge Vaquero-Garcia

#ifndef TTEST_H
#define TTEST_H
#include <cmath>
#include <vector>
#include <iostream>
#include "testStats.hpp"
#include "scythestat/distributions.h"

namespace MajiqStats{
    namespace details {
        struct tTestSummary {
            private:
                /** Estimate pooled degrees of freedom for Welch t-test
                 *
                 * Estimate pooled degrees of freedom for Welch t-test using
                 * Welch-Satterthwaite equation
                 *
                 * @param var1 sample variance of first sample
                 * @param var2 sample variance of second sample
                 * @param n1 sample size of first sample
                 * @param n2 sample size of second sample
                 *
                 * @note Modeled after _unequal_var_ttest_denom in SciPy
                 *
                 * @returns degrees of freedom, denominator for t-statistic
                 */
                inline std::pair<double, double> DOFWelchSatterthwaite(
                        double var1, double var2,
                        int n1, int n2
                ) {
                    // compute ratio of variance to sample size
                    const double ratio1 = var1 / n1;
                    const double ratio2 = var2 / n2;
                    // compute components of dof
                    const double numerator_term = ratio1 + ratio2;
                    const double denom1 = ratio1 * ratio1 / (n1 - 1);
                    const double denom2 = ratio2 * ratio2 / (n2 - 1);
                    // estimated degrees of freedom are
                    const double dof = numerator_term * numerator_term / (denom1 + denom2);
                    // denominator for t-test is
                    const double t_denom = std::sqrt(ratio1 + ratio2);
                    // return pair of dof, t_denom
                    return std::make_pair(dof, t_denom);
                }

            public:
                double dof;  // estimate of degrees of freedom of t-statistic
                double t;  // t-statistic

                /** Process data and labels to produce summary of data for testing
                 *
                 * @param data vector of observations in sorted order
                 * @param labels vector of binary labels corresponding to data
                 *
                 * @note calculate mean/variance per group in single-pass using
                 * Welford's algorithm, which is numerically stable
                 * @note t-statistic/dof are calculated assuming unequal
                 * variances as per Welch's t-test
                 */
                tTestSummary(
                        std::vector<float>& data,
                        std::vector<int>& labels
                ) {
                    // determine n, sample mean, sample variance (with Bessel's correction)
                    // first pass: get n1, n2, sum1, sum2 to caclulate first moment
                    int n1 = 0;
                    int n2 = 0;
                    double mean1 = 0.;
                    double mean2 = 0.;
                    double rss1 = 0.;
                    double rss2 = 0.;
                    for (unsigned int idx = 0; idx < data.size(); ++idx) {
                        // handle missing observation
                        if (data[idx] < 0) {
                            continue;
                        }
                        // increment ni/sumi associated with observation
                        switch(labels[idx]) {
                            case 0: {
                                ++n1;
                                const double dx = data[idx] - mean1;
                                mean1 += dx / n1;
                                rss1 += dx * (data[idx] - mean1);
                                break;
                            }
                            case 1: {
                                ++n2;
                                const double dx = data[idx] - mean2;
                                mean2 += dx / n2;
                                rss2 += dx * (data[idx] - mean2);
                                break;
                            }
                            default: {
                                std::cerr << "Invalid label passed to tTestSummary\n";
                                break;
                            }
                        }
                    }
                    // if not enough data, no need to go further
                    if (n1 < 2 || n2 < 2) {
                        // we do not have enough samples
                        dof = 0.;
                        t = 0.;
                        return;
                    }
                    // // get sample variance using Bessel's correction
                    const double var1 = rss1 / (n1 - 1);
                    const double var2 = rss2 / (n2 - 1);
                    // pool variance terms together
                    const std::pair<double, double> pooled_variance_stats =
                        DOFWelchSatterthwaite(var1, var2, n1, n2);
                    // final summary statistics
                    dof = pooled_variance_stats.first;
                    t = (mean1 - mean2) / pooled_variance_stats.second;
                }
        };

    }

    class tTest: public MajiqStats::TestStat{
        public:
            /** two-sided p-value for t statistic with given degrees of freedom
             *
             * @note does not check for invalid degrees of freedom
             */
            double TwoSidedPValue(double dof, double t) {
                return 2 * scythe::pt(-std::fabs(t), dof);
            }

            /** Calculate two-sided p-value for t-test given data/labels
             */
            double Calc_pval(
                    std::vector<float>& data, std::vector<int>& labels, float* score
            ) {
                // extract summary statistics from data
                details::tTestSummary summary(data, labels);
                // // set score to value of t statistic (unfortunately useless
                // // without degrees of freedom)
                // *score = static_cast<float>(summary.t);
                // determine p-value using statistics
                if (summary.dof > 0) {
                    return TwoSidedPValue(summary.dof, summary.t);
                } else {
                    return 1.;  // if not enough data, set p-value to 1
                }
            }
    };
}
#endif
