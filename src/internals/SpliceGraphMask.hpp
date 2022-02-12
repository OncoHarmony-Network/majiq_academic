/**
 * SpliceGraphMask.hpp
 *
 * General boolean mask over splicegraph junctions and introns
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_SPLICEGRAPHMASK_HPP
#define MAJIQ_SPLICEGRAPHMASK_HPP

#include "SpliceGraphValues.hpp"

namespace majiq {

// for now, concrete specialization, in the future may want to make child class
static_assert(sizeof(bool) == sizeof(char),
    "Using char instead of bool (vector<bool> is not a container)"
    " but did not expect size mismatch between bool and char");
using SpliceGraphMask = detail::SpliceGraphValues<char>;

}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPHMASK_HPP
