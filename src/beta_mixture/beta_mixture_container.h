#ifndef MIXTURE_CONTAINER_H
#define MIXTURE_CONTAINER_H

#include <limits>
#include <cmath>
#include <vector>
#include "numpy/ndarraytypes.h"
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>

#define CDF_PRECISION 17 // bits of precision, has to be higher than QUANTILE_PRECISION
#define QUANTILE_PRECISION 14  // bits of precision
#define QUANTILE_MAXITER 40  // maximum number of iterations to do


namespace BetaMixture {
    namespace detail {
        template<typename _npy_type>
        struct Moments {
            const _npy_type mean;
            const _npy_type variance;
            Moments(_npy_type mean, _npy_type variance)
                : mean(mean), variance(variance) {
            }
        };

        template<typename _npy_type>
        Moments<_npy_type> _beta_moments(_npy_type a, _npy_type b) {
            const auto mean = a / (a + b);
            const auto variance = mean * (1 - mean) / (1 + a + b);
            return Moments<_npy_type>(mean, variance);
        }

        template<typename _npy_type>
        inline _npy_type& _access(char *x) {
            return *reinterpret_cast<_npy_type*>(x);
        }
        inline char* _offset(char *x, npy_intp stride) {
            return x + stride;
        }
        inline char* _offset(char *x, npy_intp stride, npy_intp offset) {
            return _offset(x, stride * offset);
        }

        template<typename _npy_type>
        class BetaDistributionMixture {
            typedef boost::math::policies::policy<boost::math::policies::digits2<CDF_PRECISION>> policy_cdf;
            typedef boost::math::beta_distribution<_npy_type, policy_cdf> component_beta;

            private:
                npy_intp _size;
                char *a_buff;
                const npy_intp a_stride;
                char *b_buff;
                const npy_intp b_stride;

            public:
                BetaDistributionMixture(npy_intp size,
                        char *a_buff, npy_intp a_stride,
                        char *b_buff, npy_intp b_stride)
                        : _size(size),
                          a_buff(a_buff),
                          a_stride(a_stride),
                          b_buff(b_buff),
                          b_stride(b_stride) {
                }
                npy_intp size() const { return _size; }

            // access of individual components
            public:
                inline _npy_type alpha(npy_intp idx) const {
                    return _access<_npy_type>(_offset(a_buff, a_stride, idx));
                }
                inline _npy_type beta(npy_intp idx) const {
                    return _access<_npy_type>(_offset(b_buff, b_stride, idx));
                }
                inline component_beta operator()(npy_intp idx) const {
                    return component_beta(alpha(idx), beta(idx));
                }
                Moments<_npy_type> idx_moments(npy_intp idx) const {
                    return _beta_moments<_npy_type>(alpha(idx), beta(idx));
                }
                /**
                 * log of beta function for distribution parameters
                 */
                inline _npy_type idx_logbeta(npy_intp idx) const {
                    using boost::math::lgamma;
                    const auto a = alpha(idx);
                    const auto b = beta(idx);
                    return lgamma<_npy_type>(a) + lgamma<_npy_type>(b)
                        - lgamma<_npy_type>(a + b);
                }

            private:
                /**
                 * quick check if all alpha/beta values are valid
                 */
                bool is_invalid() const {
                    for (npy_intp idx = 0; idx < size(); ++idx) {
                        if (alpha(idx) < 0 || beta(idx) < 0) {
                            return true;
                        }
                    }
                    return false;
                }

            // moments of distribution
            public:
                /**
                 * first two central moments of mixture distribution
                 */
                Moments<_npy_type> moments() const {
                    // accumulator on means of component distribution
                    _npy_type mean_mean = 0;
                    _npy_type rss_mean = 0;
                    // accumulator on variances of component distribution
                    _npy_type mean_variance = 0;
                    // accumulate values over different beta distributions
                    for (npy_intp idx = 0; idx < size(); ++idx) {
                        const auto a = alpha(idx);
                        const auto b = beta(idx);
                        if (a < 0 || b < 0) {
                            _npy_type invalid = std::numeric_limits<_npy_type>::quiet_NaN();
                            return Moments<_npy_type>(invalid, invalid);
                        }
                        const auto moments = _beta_moments<_npy_type>(a, b);
                        // update accumulators
                        const auto n = idx + 1;
                        const auto delta_mean = moments.mean - mean_mean;
                        mean_mean += delta_mean / n;
                        rss_mean += delta_mean * (moments.mean - mean_mean);
                        mean_variance += (moments.variance - mean_variance) / n;
                    }
                    return Moments<_npy_type>(mean_mean, mean_variance + rss_mean / size());
                }

            // cdf of distribution
            private:
                /**
                 * Probability of rv from mixture distribution <= x
                 */
                inline _npy_type cdf_nocheck_(_npy_type x) const {
                    // accumulate unnormalized cdf
                    _npy_type unnormalized = 0;
                    for (npy_intp idx = 0; idx < size(); ++idx) {
                        unnormalized += boost::math::cdf(operator()(idx), x);
                    }
                    return unnormalized / size();
                }
            public:
                /**
                 * Probability of rv from mixture distribution <= x
                 */
                _npy_type cdf(_npy_type x) const {
                    // handle domain
                    if (x < 0) {
                        return 0;
                    } else if (x > 1) {
                        return 1;
                    } else if (is_invalid()) {
                        return std::numeric_limits<_npy_type>::quiet_NaN();
                    } else {
                        return cdf_nocheck_(x);
                    }
                }

            private:
                /**
                 * Probability density of mixture distribution at x
                 */
                _npy_type pdf_nocheck_(_npy_type x) const {
                    // accumulate unnormalized pdf
                    _npy_type unnormalized = 0;
                    for (npy_intp idx = 0; idx < size(); ++idx) {
                        unnormalized += boost::math::pdf(operator()(idx), x);
                    }
                    return unnormalized / size();
                }
            public:
                /**
                 * Probability density of mixture distribution at x
                 */
                _npy_type pdf(_npy_type x) const {
                    // handle domain
                    if (x < 0 || x > 1) {
                        return 0;
                    } else if (is_invalid()) {
                        return std::numeric_limits<_npy_type>::quiet_NaN();
                    } else {
                        return pdf_nocheck_(x);
                    }
                }

                /**
                 * Quantile function of mixture distribution at q (x st q = cdf(x))
                 */
                _npy_type quantile(_npy_type q) const {
                    if (q < 0 || q > 1 || is_invalid()) {
                        return std::numeric_limits<_npy_type>::quiet_NaN();
                    }
                    // Step 1: find small interval that includes cdf(x)=q
                    _npy_type lb, ub;
                    _npy_type x;  // first guess
                    // use Chebyshev inequality to bound x
                    {
                        const _npy_type max_tail = q < 0.5 ? q : 1 - q;
                        const auto m = moments();  // mean/variance
                        const _npy_type max_deviation = std::sqrt(m.variance / max_tail);
                        // we know that true answer is between
                        lb = m.mean - max_deviation;
                        ub = m.mean + max_deviation;
                        if (lb < 0) {
                            lb = 0;
                        }
                        if (ub > 1) {
                            ub = 1;
                        }
                        // pick starting point on interval
                        x = lb + q * (ub - lb);
                    }
                    // cheap derivative of cdf using cached evaluation of beta denominators
                    std::vector<_npy_type> logbetas(size());
                    for (npy_intp idx = 0; idx < size(); ++idx) {
                        logbetas[idx] = idx_logbeta(idx);
                    }
                    // derivative of cdf
                    const auto f = [&](_npy_type x) -> _npy_type {
                        const _npy_type logx = std::log(x);
                        const _npy_type logmx = std::log(1 - x);
                        std::vector<_npy_type> logpdfs(size());
                        _npy_type logpdfs_max = 0;  // get max for logsumexp
                        for (npy_intp idx = 0; idx < size(); ++idx) {
                            logpdfs[idx] = (alpha(idx) - 1) * logx
                                + (beta(idx) - 1) * logmx
                                - logbetas[idx];
                            if (idx == 0 || logpdfs[idx] > logpdfs_max) {
                                logpdfs_max = logpdfs[idx];
                            }
                        }
                        // add values
                        _npy_type sum = 0;
                        for (npy_intp idx = 0; idx < size(); ++idx) {
                            sum += std::exp(logpdfs[idx] - logpdfs_max);
                        }
                        sum *= std::exp(logpdfs_max);  // put the max term back in
                        sum /= size();  // normalize it
                        return sum;
                    };
                    using boost::math::tools::newton_raphson_iterate;
                    uintmax_t maxiter = QUANTILE_MAXITER;
                    return newton_raphson_iterate(
                            [&](_npy_type x) { return std::make_pair(cdf(x) - q, f(x)); },
                            x, lb, ub, QUANTILE_PRECISION, maxiter);
                }
        };
    }
}

#endif
