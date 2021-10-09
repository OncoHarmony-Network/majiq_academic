/**
 * DPsiDiscrete.hpp
 *
 * Compute discretized posterior distribution from two input distributions
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQINCLUDE_DPSIDISCRETE_HPP
#define MAJIQINCLUDE_DPSIDISCRETE_HPP

#include "BetaMixture.hpp"
#include "LogMath.hpp"

namespace MajiqInclude {
namespace DPsiDiscrete {

/**
 * Collapse log probabilities for psi1, psi2 to overlapping probabilities for
 * difference in PSI (i.e. logsumexp over diagonals).
 *
 * If there are n bins for PSI, we consider them a uniform distribution on that
 * bin. Then, taking diagonals of the outer product of discretized
 * distribution, we get overlapping triangles on [-1, -1 + 2 / n], [-1 + 1 / n,
 * -1 + 3 / n], ...
 *  Thus, we expect exactly 2 * n_psibins - 1 outputs starting at logp_dpsi.
 *
 *  This function gets the weight on each triangle. It does not attempt to
 *  interpolate out the overlaps to independent bins.
 */
template <typename ItPsi1, typename ItPsi2, typename ItDPsi>
inline void DPsiConvolve(
    ItPsi1 logp_psi1, ItPsi2 logp_psi2, const int64_t n_psibins,
    // fills 2 * n_psibins - 1 values of logp_dpsi
    ItDPsi logp_dpsi) {
  for (int64_t offset = 1 - n_psibins; offset < n_psibins;
      ++offset, ++logp_dpsi) {
    using MajiqInclude::detail::logsumexp_diag;
    *logp_dpsi = logsumexp_diag(logp_psi1, logp_psi2, n_psibins, offset);
  }
}

/**
 * Collapse psi distributions to dpsi distribution treating as independent
 *
 * Assumes that output dpsi bins are exactly 2 * input psi bins, which allows
 * us to interpolate by simple averaging (odd number of bins requires
 * per-position averaging of adjacent weights after DPsiConvolve)
 *
 * Does not normalize result because needs to be multiplied by prior anyway.
 */
template <typename ItPsi1, typename ItPsi2, typename ItDPsi>
inline void DPsiIndependent(
    ItPsi1 logp_psi1, ItPsi2 logp_psi2, const int64_t n_psibins,
    // fills 2 * n_psibins values of logp_dpsi
    ItDPsi logp_dpsi) {
  // DPsiConvolve into bins after first one
  DPsiConvolve(logp_psi1, logp_psi2, n_psibins, 1 + logp_dpsi);
  {  // average adjacent bins
    logp_dpsi[0] = logp_dpsi[1];
    for (int64_t i = 1; i < 2 * n_psibins - 1; ++i) {
      using MajiqInclude::detail::logadd;
      logp_dpsi[i] = logadd(logp_dpsi[i], logp_dpsi[1 + i]);
    }
  }
  return;
}

/**
 * Apply prior to distribution and normalize
 *
 * logp_dpsi at input should be independent dpsi, at output, normalized after
 * prior applied
 */
template <typename ItDPsi, typename ItPrior>
inline void DPsiPosterior(ItDPsi logp_dpsi, ItPrior logp_prior,
    const int64_t n_dpsibins) {
  // apply prior
  for (int64_t i = 0; i < n_dpsibins; ++i) {
    logp_dpsi[i] += logp_prior[i];
  }
  // normalize
  using MajiqInclude::detail::lognormalize;
  lognormalize(logp_dpsi, n_dpsibins);
  return;
}

}  // namespace DPsiDiscrete
}  // namespace MajiqInclude
#endif  // MAJIQINCLUDE_DPSIDISCRETE_HPP
