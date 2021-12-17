/**
 * constants.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PYBIND_CONSTANTS_HPP
#define MAJIQ_PYBIND_CONSTANTS_HPP

#include "../MajiqTypes.hpp"
#include "../ExperimentThresholds.hpp"

namespace majiq {
namespace bindings {

constexpr bool DEFAULT_BUILD_PROCESS_IR = true;
constexpr junction_ct_t DEFAULT_BUILD_MINREADS = 3;
constexpr junction_ct_t DEFAULT_BUILD_MINDENOVO = 5;
constexpr junction_pos_t DEFAULT_BUILD_MINPOS = 2;
constexpr real_t DEFAULT_BUILD_MAX_PCTBINS = 0.6;
constexpr real_t DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY = 0.5;
constexpr real_t DEFAULT_BUILD_MATCH_INTRON_PROBABILITY = 0.95;
constexpr real_t DEFAULT_BUILD_MINEXPERIMENTS = 0.5;
constexpr bool DEFAULT_BUILD_DENOVO_JUNCTIONS = true;
constexpr bool DEFAULT_BUILD_DENOVO_IR = true;
constexpr bool DEFAULT_BUILD_KEEP_ANNOTATED_IR = false;
constexpr size_t DEFAULT_BUILD_NUM_BOOTSTRAPS = 30;
constexpr real_t DEFAULT_BUILD_STACK_PVALUE = 1e-7;
constexpr ExperimentStrandness DEFAULT_BAM_STRANDNESS
  = ExperimentStrandness::NONE;
constexpr int DEFAULT_BAM_NTHREADS = 1;
constexpr real_t DEFAULT_BUILD_SIMPL_MINPSI = 0.01;
constexpr real_t DEFAULT_BUILD_SIMPL_MINREADS_ANNOTATED_JUNCTION = 0;
constexpr real_t DEFAULT_BUILD_SIMPL_MINREADS_DENOVO_JUNCTION = 0;
constexpr real_t DEFAULT_BUILD_SIMPL_MINREADS_INTRON = 0;
//
// default value of ExperimentThresholds
static const auto DEFAULT_THRESHOLDS = ExperimentThresholds(
    DEFAULT_BUILD_MINREADS, DEFAULT_BUILD_MINDENOVO, DEFAULT_BUILD_MINPOS,
    DEFAULT_BUILD_MAX_PCTBINS, DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
    DEFAULT_BUILD_MATCH_INTRON_PROBABILITY);

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_CONSTANTS_HPP
