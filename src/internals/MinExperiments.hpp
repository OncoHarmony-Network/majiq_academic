/**
 * MinExperiments.hpp
 *
 * How MAJIQ handles min-experiments
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_MINEXPERIMENTS_HPP
#define MAJIQ_MINEXPERIMENTS_HPP

#include <stdexcept>
#include <algorithm>


namespace majiq {
namespace detail {

inline size_t min_experiments_from_float(
    size_t num_experiments, real_t min_experiments_f) {
  // get min_experiments as size_t from input
  if (min_experiments_f < 0) {
    throw std::invalid_argument("min_experiments must be non-negative");
  } else if (min_experiments_f < 1) {
    // less than 1 ~ percentage of number of experiments
    min_experiments_f *= num_experiments;
  }
  // go to next number of experiments, max being group.num_experiments_
  return std::min(
      static_cast<size_t>(std::ceil(min_experiments_f)),
      num_experiments);
}

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_MINEXPERIMENTS_HPP
