/**
 * MaijqTypes.hpp
 *
 * base types for MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_TYPES_HPP
#define MAJIQ_TYPES_HPP

#include <cstdint>

namespace majiq {

using real_t = float;  // for numerics
using position_t = int64_t;  // for coordinates, difference in coordinates

enum class GeneStrandness : char {
  FORWARD = '+',
  REVERSE = '-',
  AMBIGUOUS = '.',
};

enum class ExperimentStrandness : char {
  FORWARD = 'F',  // read 1 is in forward direction, Salmon ISF
  REVERSE = 'R',  // read 2 is in forward direction, Salmon ISR
  NONE = 'N',  // could be either way, Salmon IU
};

}  // namespace majiq

#endif  // MAJIQ_TYPES_HPP
