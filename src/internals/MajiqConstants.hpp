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
static_assert(
    MAX_DENOVO_DIFFERENCE >= 0, "MAX_DENOVO_DIFFERENCE must be non-negative");
// how many bases does alignment have to overlap exons?
constexpr junction_pos_t USE_MIN_OVERHANG = 8;
static_assert(USE_MIN_OVERHANG > 0, "USE_MIN_OVERHANG must be positive");
// for historical reaasons, we store event descriptions in fixed-width
// character array of this size
// TODO (jaicher) maybe use variablle-length serialization?
constexpr int64_t EVENT_DESCRIPTION_WIDTH = 250;
static_assert(EVENT_DESCRIPTION_WIDTH > 6,
    "EVENT_DESCRIPTION_WIDTH must be long enough to hold s|na|i");

}  // namespace majiq


#endif  // MAJIQ_CONSTANTS_HPP
