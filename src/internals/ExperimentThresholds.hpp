/**
 * ExperimentThresholds.hpp
 *
 * Class for determining thresholds for acceptance of MAJIQ introns of various
 * lengths by matching junction acceptance threshold probability with identical
 * per-position readrate.
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Model: each position on a transcript feature (junction, intron) has the same
 * underlying readrate. We observe reads for each position as a Poisson random
 * variable with this underlying readrate.
 *
 * For junctions, there are P (aka eff_len) positions mapped (by identity) to P
 * bins. However, for introns of length L, there are L additional positions
 * that must be distributed to the P positions.
 *
 * We accept junctions on the condition that they have at least minreads reads
 * observed across all bins and minpos bins with at least 1 read. For some
 * acceptance probability q, we can backsolve for minimum readrate required to
 * pass thresholds with that probability.
 *
 * Then, we can backsolve the intron thresholds for each length L that pass
 * with just under the same probability q when given equivalent per-position
 * readrate.
 */

#ifndef MAJIQ_EXPERIMENTTHRESHOLDS_HPP
#define MAJIQ_EXPERIMENTTHRESHOLDS_HPP

#include <algorithm>
#include <cmath>
#include <utility>
#include <limits>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/tools/roots.hpp>

#include "MajiqTypes.hpp"


namespace majiq {
namespace detail {

/**
 * Readrate with at least minreads with sufficient probability
 */
inline real_t ReadrateMinreads(
    junction_ct_t minreads,
    real_t match_probability) {
  // edge cases
  if (match_probability >= 1) {
    return std::numeric_limits<real_t>::infinity();
  } else if (match_probability <= 0 || minreads <= 0) {
    return 0;
  }

  // use root-finding algorithm to solve otherwise
  using boost::math::poisson_distribution;
  using boost::math::pdf;
  using boost::math::cdf;
  using boost::math::complement;
  using boost::math::tools::newton_raphson_iterate;

  // objective, derivative of objective to use root-finding methods
  auto objective_pair = [minreads, match_probability](real_t rate) {
    poisson_distribution<real_t> dist(rate);  // current distribution
    real_t objective =
      // sf(minreads - 1) - match_probability = 0
      cdf(complement(dist, minreads - 1))
      - match_probability;
    // turns out that derivative of cdf wrt lambda is -pmf, so...
    real_t dobjective = pdf(dist, minreads - 1);
    return std::make_tuple(objective, dobjective);
  };
  // Newton iteration to solve total readrate with match_probability
  const real_t guess = minreads;
  const real_t bound_min = 0;
  const real_t bound_max = 10 * guess;
  uintmax_t niter = 100;  // number of iterations to limit to
  const real_t total_readrate = newton_raphson_iterate(
      objective_pair, guess, bound_min, bound_max,
      0.4 * std::numeric_limits<real_t>::digits, niter);
  // divide total readrate evenly acrosss bins
  return total_readrate;
}

/**
 * Per-position readrate meeting minreads with sufficient probability
 */
inline real_t ReadrateJunctionMinreads(
    junction_ct_t minreads,
    junction_pos_t numbins,
    real_t acceptance_probability) {
  return ReadrateMinreads(minreads, acceptance_probability) / numbins;
}

/**
 * Probability of a bin having at least mincov reads that gives desired
 * acceptance_probability when requiring specified minbins no less than
 * mincov
 */
inline real_t ProbBinMincovFromBinomial(
    junction_pos_t minbins,
    junction_pos_t numbins,
    real_t acceptance_probability) {
  using boost::math::binomial_distribution;
  return binomial_distribution<real_t>::find_upper_bound_on_p(
      numbins, minbins - 1, 1 - acceptance_probability);
}

/**
 * Per-position readrate meeting minpos with sufficient probability
 */
inline real_t ReadrateJunctionMinpos(
    junction_pos_t minpos,
    junction_pos_t numbins,
    real_t acceptance_probability) {
  // the nonzero probability that we need per position is:
  real_t nonzero_prob = ProbBinMincovFromBinomial(minpos, numbins,
      acceptance_probability);
  // we know that this is equal to 1 - exp(-readrate), so invert:
  return -std::log(1 - nonzero_prob);
}

/**
 * Per-position readrate meeting intron thresholds for given length
 */
inline real_t ReadrateIntronThresholds(
    junction_pos_t minbins,
    junction_ct_t mincov,
    junction_pos_t numbins,
    junction_pos_t intron_length,
    real_t acceptance_probability) {
  // per-bin probability for minimum coverage needed
  real_t bin_prob = ProbBinMincovFromBinomial(minbins, numbins,
      acceptance_probability);
  // per-bin readrate
  real_t bin_readrate = ReadrateMinreads(mincov, bin_prob);
  // rescale to positions
  return (bin_readrate * numbins) / (numbins + intron_length);
}

/**
 * Longest intron length with nonzero coverage probability less than input
 */
inline junction_pos_t MaxLengthNonzeroCoverage(
    real_t readrate,
    junction_pos_t numbins,
    real_t prob_mincov) {
  // get rate for which Poisson probability of nonzero is prob_mincov
  // (solve: 1 - exp(-lambda) = prob_mincov)
  const real_t exact_binrate = -std::log(1 - prob_mincov);
  // solve for greatest intron length where the per-bin readrate is no
  // greater than exact_binrate.
  // The per-bin readrate is (P+L)/P*readrate, so we have:
  return static_cast<junction_pos_t>(
      std::floor(numbins * exact_binrate / readrate - numbins));
}

/**
 * total readrate aggregated over all intron positions
 */
inline real_t TotalReadrate(real_t readrate, junction_pos_t numbins,
    junction_pos_t intron_length) {
  return readrate * (numbins + intron_length);
}

/**
 * readrate for intron bin given positional readrate assuming fractional
 * positions per bin: (numbins + intron_length) / numbins
 */
inline real_t BinReadrate(
    real_t readrate, junction_pos_t numbins, junction_pos_t intron_length) {
  return TotalReadrate(readrate, numbins, intron_length) / numbins;
}

/**
 * The probability that an intron bin has at least mincov reads given Poisson
 * rate per position, number of bins, and intron length
 */
inline real_t ProbBinMincovFromPoisson(
    real_t readrate,
    junction_pos_t numbins,
    junction_pos_t intron_length,
    junction_ct_t mincov) {
  // edge case
  if (mincov <= 0) {
    return 1;  // will always have at least 0 reads
  }

  // calculate Poisson survival function with appropriate bin readrate
  using boost::math::poisson_distribution;
  using boost::math::cdf;
  using boost::math::complement;
  // distribution associated with the readrate for an intron bin
  real_t bin_readrate = BinReadrate(readrate, numbins, intron_length);
  poisson_distribution<real_t> dist(bin_readrate);
  return cdf(complement(dist, mincov - 1));
}

/**
 * Smallest value of minbins that has acceptance probability no greater than
 * acceptance_probability given readrate, total number of bins, intron
 * length, and mincov = 1
 */
inline junction_pos_t MinbinsFromAcceptance(
    real_t readrate,
    junction_pos_t numbins,
    real_t acceptance_probability,
    junction_pos_t intron_length) {
  // get probability that a bin is nonzero
  constexpr junction_ct_t mincov_for_nonzero = 1;  // mincov for nonzero
  real_t prob_nonzero = ProbBinMincovFromPoisson(readrate, numbins,
      intron_length, mincov_for_nonzero);

  // consider distribution of bins with nonzero coverage
  // round complement quantiles down (we want under acceptance_probability)
  using boost::math::binomial_distribution;
  using boost::math::policies::policy;
  using boost::math::policies::discrete_quantile;
  using boost::math::policies::integer_round_up;
  using boost::math::quantile;
  using boost::math::complement;
  binomial_distribution<
    real_t, policy<discrete_quantile<integer_round_up>>
    > dist_nonzero(numbins, prob_nonzero);
  // min { k : P(NonzeroBins >= k) <= q } (1 + since isf for > k, not >= k)
  return 1 + quantile(complement(dist_nonzero, acceptance_probability));
}

/**
 * Smallest value of minreads that has probability no greater than
 * acceptance_probability given total readrate
 */
inline junction_ct_t MinreadsFromAcceptance(
    real_t total_readrate,
    real_t acceptance_probability) {
  // consider the total count distribution
  // round complement quantiles down (want just under acceptance_probability)
  using boost::math::poisson_distribution;
  using boost::math::policies::policy;
  using boost::math::policies::discrete_quantile;
  using boost::math::policies::integer_round_up;
  using boost::math::quantile;
  using boost::math::complement;
  poisson_distribution<
    real_t, policy<discrete_quantile<integer_round_up>>
    > dist(total_readrate);
  // min { k : P(X >= k) <= prob } (1 + since isf for > k, not >= k)
  return 1 + quantile(complement(dist, acceptance_probability));
}

/**
 * Smallest value of mincov that has probability no greater than
 * prob_bin_mincov given readrate, total number of bins, intron length
 */
inline junction_ct_t MincovFromAcceptance(
    real_t readrate,
    junction_pos_t numbins,
    junction_pos_t intron_length,
    real_t prob_bin_mincov) {
  const real_t bin_readrate = BinReadrate(readrate, numbins, intron_length);
  return MinreadsFromAcceptance(bin_readrate, prob_bin_mincov);
}

/**
 * Intron acceptance probability for provided thresholds
 */
inline real_t IntronAcceptanceProbability(
    real_t readrate,
    junction_pos_t numbins,
    junction_pos_t intron_length,
    junction_pos_t minbins,
    junction_ct_t mincov) {
  // probability of a bin having at least minimum coverage
  real_t prob_mincov = ProbBinMincovFromPoisson(readrate, numbins,
      intron_length, mincov);
  // distribution of number bins with at least min coverage
  using boost::math::binomial_distribution;
  using boost::math::cdf;
  using boost::math::complement;
  binomial_distribution<real_t> dist_numpassed(numbins, prob_mincov);
  return cdf(complement(dist_numpassed, minbins - 1));
}

/**
 * Junction minreads acceptance probability
 */
inline real_t JunctionMinreadsAcceptanceProbability(
    real_t readrate,
    junction_pos_t numbins,
    junction_ct_t minreads) {
  // total readrate
  constexpr junction_pos_t equivalent_intron_length = 0;
  const real_t total_readrate = TotalReadrate(readrate, numbins,
      equivalent_intron_length);
  // Poisson distribution with total readrate...
  using boost::math::poisson_distribution;
  using boost::math::cdf;
  using boost::math::complement;
  poisson_distribution<real_t> dist(total_readrate);
  // sf(minreads - 1)
  return cdf(complement(dist, minreads - 1));
}

/**
 * Junction minpos acceptance probability
 */
inline real_t JunctionMinposAcceptanceProbability(
    real_t readrate,
    junction_pos_t numbins,
    junction_pos_t minpos) {
  // per-bin/position probability of nonzero reads
  constexpr junction_pos_t equivalent_intron_length = 0;
  constexpr junction_ct_t mincov = 1;
  return IntronAcceptanceProbability(readrate, numbins,
      equivalent_intron_length, minpos, mincov);
}

}  // namespace detail

struct IntronThresholds {
  junction_ct_t minreads_;
  junction_pos_t minbins_;  // minimum bins with at least
  junction_ct_t mincov_;  // mincov reads for intron to pass
};

class IntronThresholdsGenerator {
 private:
  // total number of bins present
  const junction_pos_t total_bins_;
  // maximum value for minbins_
  const junction_pos_t max_minbins_;
  // acceptance probability being matched
  const real_t acceptance_probability_;
  // per-position readrate being matched
  const real_t readrate_;
  // binomial probability for max_minbins to match acceptance
  const real_t prob_bin_mincov_;
  // maximum length where mincov = 1 (use different algorithm afterwards)
  const junction_pos_t max_length_nonzerocov_;
  // remember minreads because we also use as intron threshold
  const junction_ct_t minreads_;

 public:
  /**
   * @param minreads, minpos junction thresholds
   * @param total_bins number of bins for each feature (also positions for
   * junctions)
   * @param max_pctbins maximum percent of total bins for minbins
   * @param junction_acceptance_probability determine positional readrate so
   * that junctions have this acceptance probability given junction thresholds
   * @param intron_acceptance_probability intron thresholds set so that
   * determined positional readrate has just under this acceptance probability
   */
  IntronThresholdsGenerator(
      junction_ct_t minreads, junction_pos_t minpos,
      junction_pos_t total_bins, real_t max_pctbins,
      real_t junction_acceptance_probability,
      real_t intron_acceptance_probability)
      : total_bins_{total_bins},
        max_minbins_{
          std::max(1,
              static_cast<junction_pos_t>(
                total_bins_ * std::clamp(max_pctbins, real_t{0}, real_t{1})))},
        acceptance_probability_{intron_acceptance_probability},
        readrate_{std::max(
            detail::ReadrateJunctionMinreads(
              minreads, total_bins_, junction_acceptance_probability),
            detail::ReadrateJunctionMinpos(
              minpos, total_bins_, junction_acceptance_probability))},
        prob_bin_mincov_{detail::ProbBinMincovFromBinomial(
            max_minbins_, total_bins_, acceptance_probability_)},
        max_length_nonzerocov_{detail::MaxLengthNonzeroCoverage(
            readrate_, total_bins_, prob_bin_mincov_)},
        minreads_{minreads} { }

  /**
   * intron thresholds that match junction threshold acceptance_probability for
   * the same positional readrate
   */
  IntronThresholds operator()(junction_pos_t intron_length) const {
    if (intron_length <= max_length_nonzerocov_) {
      return {
        minreads_,
        detail::MinbinsFromAcceptance(
            readrate_, total_bins_, acceptance_probability_, intron_length),
        1};
    } else {
      return {
        minreads_,
        max_minbins_,
        detail::MincovFromAcceptance(
            readrate_, total_bins_, intron_length, prob_bin_mincov_)};
    }
  }
};

class ExperimentThresholds {
 public:
  const junction_ct_t minreads_;
  const junction_ct_t mindenovo_;
  const junction_pos_t minpos_;

  ExperimentThresholds(
      junction_ct_t minreads, junction_ct_t mindenovo, junction_pos_t minpos)
      : minreads_{std::max(1, minreads)},
        mindenovo_{std::max(minreads_, mindenovo)},
        minpos_{std::max(1, minpos)} { }


  IntronThresholdsGenerator intron_thresholds_generator(
      junction_pos_t total_bins, real_t max_pctbins,
      real_t junction_acceptance_probability,
      real_t intron_acceptance_probability) {
    return IntronThresholdsGenerator(
        mindenovo_, minpos_, total_bins, max_pctbins,
        junction_acceptance_probability, intron_acceptance_probability);
  }
};

}  // namespace majiq
#endif  // MAJIQ_EXPERIMENTTHRESHOLDS_HPP
