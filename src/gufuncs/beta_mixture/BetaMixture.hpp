/**
 * BetaMixture.hpp
 *
 * Compute beta mixture around distribution
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <boost/math/distributions/beta.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace MajiqGufuncs {
namespace BetaMixture {

constexpr int CDF_DIGITS2 = 17;
constexpr int QUANTILE_DIGITS2 = 14;
using fast_policy = boost::math::policies::policy<
  boost::math::policies::promote_double<false>,
  boost::math::policies::promote_float<true>,
  boost::math::policies::digits2<CDF_DIGITS2>
  >;
template <typename RealT>
using DistT = boost::math::beta_distribution<RealT, fast_policy>;

// still assumes n_mixture > 0. If we have zero-length axis, handle at outer
// loop
template <typename RealT>
inline bool IsInvalidComponent(RealT a, RealT b) {
  return npy_isnan(a) || npy_isnan(b) || a <= 0 || b <= 0;
}
template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline bool IsInvalid(ItA a, ItB b, const npy_intp n_mixture) {
  for (npy_intp i = 0; i < n_mixture; ++i, ++a, ++b) {
    if (IsInvalidComponent(*a, *b)) {
      return true;
    }
  }
  return false;
}

template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _Mean(ItA a, ItB b, const npy_intp n_mixture) {
  RealT sum{0};
  for (npy_intp i = 0; i < n_mixture; ++i, ++a, ++b) {
    if (IsInvalidComponent(*a, *b)) {
      return std::numeric_limits<RealT>::quiet_NaN();
    }
    sum += *a / (*a + *b);
  }
  return sum / n_mixture;
}

template <typename RealT>
struct CentralMoments {
  RealT mean;
  RealT variance;
};

/**
 * Get mean, variance of beta mixture, assuming n_mixture > 0
 */
template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline CentralMoments<RealT> _Moments(
    ItA a1, ItB b1, const npy_intp n_mixture) {
  // iterators for second pass
  ItA a2{a1};
  ItB b2{b1};
  // mean = mean_mean. variance = mean_variance + variance_mean
  RealT sum_mean{0};
  RealT sum_variance{0};
  for (npy_intp i = 0; i < n_mixture; ++i, ++a1, ++b1) {
    if (IsInvalidComponent(*a1, *b1)) {
      return CentralMoments<RealT>{
        std::numeric_limits<RealT>::quiet_NaN(),
        std::numeric_limits<RealT>::quiet_NaN()
      };
    }
    const RealT denom = *a1 + *b1;
    const RealT mean = *a1 / denom;
    sum_mean += mean;
    sum_variance += mean * (1 - mean) / (1 + denom);
  }
  const RealT mean_mean = sum_mean / n_mixture;
  const RealT mean_variance = sum_variance / n_mixture;
  // get variance_mean in second pass
  RealT rss_mean{0};
  for (npy_intp i = 0; i < n_mixture; ++i, ++a2, ++b2) {
    const RealT denom = *a2 + *b2;
    const RealT mean = *a2 / denom;
    const RealT residual = mean - mean_mean;
    rss_mean += residual * residual;
  }
  const RealT variance_mean = rss_mean / n_mixture;
  return CentralMoments<RealT>{mean_mean, mean_variance + variance_mean};
}

/**
 * Get beta distribution approximation to mixture matching mean and variance
 *
 * Assumes that mean, variance are possible to match moments with beta:
 * 0 < mean < 1, variance < mean * (1 - mean)
 */
template <typename RealT>
inline std::pair<RealT, RealT> _BetaApproximation(
    const CentralMoments<RealT>& moments) {
  if (npy_isnan(moments.mean)) {
    return std::make_pair(std::numeric_limits<RealT>::quiet_NaN(),
        std::numeric_limits<RealT>::quiet_NaN());
  }
  // https://www.johndcook.com/blog/2021/04/07/beta-given-mean-variance/
  const RealT a0 = moments.mean * (
      moments.mean * (1 - moments.mean) / moments.variance - 1);
  const RealT b0 = a0 * (1 - moments.mean) / moments.mean;
  return std::make_pair(a0, b0);
}
template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline std::pair<RealT, RealT> _BetaApproximation(
    ItA a, ItB b, const npy_intp n_mixture) {
  const CentralMoments<RealT> moments = _Moments(a, b, n_mixture);
  return _BetaApproximation(moments);
}

template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _CDF_unchecked(RealT x, ItA a, ItB b, const npy_intp n_mixture) {
  using boost::math::cdf;
  RealT sum{0};
  for (npy_intp i = 0; i < n_mixture; ++i, ++a, ++b) {
    DistT<RealT> dist{*a, *b};
    sum += cdf(dist, x);
  }
  return sum / n_mixture;
}

template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _CDF(RealT x, ItA a, ItB b, const npy_intp n_mixture) {
  if (npy_isnan(x) || IsInvalid(a, b, n_mixture)) {
    return std::numeric_limits<RealT>::quiet_NaN();
  } else if (x <= 0) {
    return RealT{0};
  } else if (x >= 1) {
    return RealT{1};
  } else {
    return _CDF_unchecked(x, a, b, n_mixture);
  }
}

/**
 * Discrete approximation of PDF on n_out uniformly spaced bins
 */
template <int digits2 = 17, typename ItA, typename ItB, typename ItOut,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline void _PMF(ItA a, ItB b, ItOut out,
    const npy_intp n_mixture, const npy_intp n_out) {
  RealT mean = _Mean(a, b, n_mixture);
  if (npy_isnan(mean)) {
    // if mean NaN, then PMF should be NaN
    if (n_out > 1) {
      std::fill(out, out + n_out, mean);
    } else {
      out[0] = mean;
    }
    return;
  } else if (n_out == 1) {
    out[0] = RealT{1};
  }
  // calculate CDF at *ends* of each bin. Start around mean and go out until
  // hitting desired tolerance
  // then, take differences (right to left) to get PMF
  constexpr RealT TOL_MIN{RealT{1} / (1 << digits2)};
  constexpr RealT TOL_MAX{1 - TOL_MIN};
  // initialize out to zeros
  std::fill(out, out + n_out, RealT{0});
  // easier to think about indexing endpoints[1 + n_out] with i
  // although we only work with out relative to endpoints[1:]. So indexing
  // of out will be on i - 1
  npy_intp i_mean = static_cast<npy_intp>(n_out * mean);
  // I don't check that 0 <= i_mean <= n_out because mean on [0, 1]
  {
    // carry index between loops (used to fill remaining points with 1 on rhs)
    npy_intp i;
    const RealT n_out_f = static_cast<RealT>(n_out);  // divide to get back x
    // go backwards from mean
    for (i = i_mean; i > 0; --i) {
      out[i - 1] = _CDF_unchecked(i / n_out_f, a, b, n_mixture);
      if (out[i - 1] <= TOL_MIN) {
        break;  // rest of i's close to zero, and they are already zero
      }
    }
    // go forwards from mean
    for (i = 1 + i_mean; i < n_out; ++i) {
      out[i - 1] = _CDF_unchecked(i / n_out_f, a, b, n_mixture);
      if (out[i - 1] >= TOL_MAX) {
        break;  // rest of i's close to one
      }
    }
    // fill remaining points with 1
    for (i = std::min(i, n_out - 1); i < n_out; ++i) {
      out[i] = RealT{1};
    }
  }
  // take differences to get discretized PMF
  for (npy_intp i = n_out - 1; i > 0; --i) {
    out[i] -= out[i - 1];
  }
  return;
}

template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _PDF_unchecked(RealT x, ItA a, ItB b, const npy_intp n_mixture) {
  using boost::math::pdf;
  RealT sum{0};
  for (npy_intp i = 0; i < n_mixture; ++i, ++a, ++b) {
    DistT<RealT> dist{*a, *b};
    sum += pdf(dist, x);
  }
  return sum / n_mixture;
}

template <typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _PDF(RealT x, ItA a, ItB b, const npy_intp n_mixture) {
  if (npy_isnan(x) || IsInvalid(a, b, n_mixture)) {
    return std::numeric_limits<RealT>::quiet_NaN();
  } else if (x <= 0 || x >= 1) {
    return RealT{0};
  } else {
    return _PDF_unchecked(x, a, b, n_mixture);
  }
}

template <uintmax_t MAX_ITER = 20, typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _Quantile(
    const RealT q, const ItA a, const ItB b, const npy_intp n_mixture) {
  // special cases to return early
  if (npy_isnan(q)) {
    return std::numeric_limits<RealT>::quiet_NaN();
  }
  const CentralMoments<RealT> moments = _Moments(a, b, n_mixture);
  if (npy_isnan(moments.mean)) {
    // if mean is nan, quantiles should be nan
    return std::numeric_limits<RealT>::quiet_NaN();
  } else if (q <= 0) {
    return RealT{0};
  } else if (q >= 1) {
    return RealT{1};
  } else if (n_mixture == 1) {
    return boost::math::quantile(DistT<RealT>{a[0], b[0]}, q);
  }

  // otherwise, get bound on x, good first guess
  RealT x_lb{1}, x_ub{0}, x{};
  {
    // use Chebyshev's inequality to bound x
    const RealT max_tail{q < 0.5 ? q : 1 - q};
    const RealT max_deviation = std::sqrt(moments.variance / max_tail);
    x_lb = std::max(RealT{0}, moments.mean - max_deviation);
    x_ub = std::min(RealT{1}, moments.mean + max_deviation);
    x = x_lb + q * 2 * (x_ub - x_lb);
    if (x_ub - x_lb <= RealT{1} / (1 << QUANTILE_DIGITS2)) {
      // our bound is tight enough that we can just return it
      return x;
    }
  }
  // otherwise, set up pdf/cdf for root finding
  const auto cdf = [&a, &b, n_mixture](RealT x) {
    if (x <= 0) {
      return RealT{0};
    } else if (x >= 1) {
      return RealT{1};
    } else {
      return _CDF_unchecked(x, a, b, n_mixture);
    }
  };
  // to compute pdf, we will reuse the denominator terms many times
  std::vector<RealT> logbetas(n_mixture);
  for (npy_intp i = 0; i < n_mixture; ++i) {
    logbetas[i]
      = std::lgamma(a[i]) + std::lgamma(b[i])
      - std::lgamma(a[i] + b[i]);
  }
  const auto logpdf_component = [&a, &b, &logbetas](
      npy_intp i, RealT logx, RealT log1mx) {
    return (a[i] - 1) * logx + (b[i] - 1) * log1mx - logbetas[i];
  };
  std::vector<RealT> _logpdf_buffer(n_mixture);
  const auto pdf = [&_logpdf_buffer, &logpdf_component, n_mixture](RealT x) {
    if (x <= 0 || x >= 1) {
      return RealT{0};
    }
    const RealT logx = std::log(x);
    const RealT log1mx = std::log1p(-x);
    // compute logpdf of each component, get max among them
    RealT logpdf_max = std::numeric_limits<RealT>::lowest();
    for (npy_intp i = 0; i < n_mixture; ++i) {
      _logpdf_buffer[i] = logpdf_component(i, logx, log1mx);
      logpdf_max = std::max(logpdf_max, _logpdf_buffer[i]);
    }
    // logsumexp trick
    RealT sum{0};
    for (npy_intp i = 0; i < n_mixture; ++i) {
      sum += std::exp(_logpdf_buffer[i] - logpdf_max);
    }
    sum *= std::exp(logpdf_max);
    return sum / n_mixture;
  };
  // newton argument
  const auto f = [&cdf, &pdf, q](RealT x) {
    return std::make_pair(
        cdf(x) - q,
        pdf(x) + std::numeric_limits<RealT>::epsilon());
  };
  using boost::math::tools::newton_raphson_iterate;
  uintmax_t maxiter = MAX_ITER;
  x = newton_raphson_iterate(f, x, x_lb, x_ub, QUANTILE_DIGITS2, maxiter);
  return x;
}

template <typename Generator, typename RealT>
inline RealT _Sample_unchecked(Generator& g, RealT a, RealT b) {
  // assumes a, b finite and positive
  using boost::random::beta_distribution;
  beta_distribution<RealT> dist{a, b};
  return dist(g);
}
template <typename Generator, typename ItA, typename ItB,
         typename RealT = typename std::iterator_traits<ItA>::value_type>
inline RealT _SampleMixture_unchecked(
    Generator& g, ItA a, ItB b, const npy_intp n_mixture) {
  // assumes a, b finite and positive
  using boost::random::uniform_int_distribution;
  uniform_int_distribution<npy_intp> dist_source(0, n_mixture - 1);
  // which beta distribution to sample from?
  const auto i = dist_source(g);
  return _Sample_unchecked(g, a[i], b[i]);
}


}  // namespace BetaMixture
}  // namespace MajiqGufuncs


#endif  // MAJIQGUFUNCS_BETAMIXTURE_HPP
