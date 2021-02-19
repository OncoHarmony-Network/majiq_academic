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
// how many bases does alignment have to overlap exons?
constexpr uint32_t USE_MIN_OVERHANG = 8;

}  // namespace majiq


#endif  // MAJIQ_CONSTANTS_HPP
