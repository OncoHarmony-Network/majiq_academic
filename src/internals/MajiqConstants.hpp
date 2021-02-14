/**
 * MajiqConstants.hpp
 *
 * constants for MAJIQ internals
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONSTANTS_HPP
#define MAJIQ_CONSTANTS_HPP

#include "MajiqTypes.hpp"

namespace majiq {

// how far can an individual denovo junction extend an exon?
constexpr position_t MAX_DENOVO_DIFFERENCE = 400;

}  // namespace majiq


#endif  // MAJIQ_CONSTANTS_HPP
