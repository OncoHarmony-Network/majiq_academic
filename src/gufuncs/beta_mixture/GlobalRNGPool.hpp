/**
 * GlobalRNGPool.hpp
 *
 * Static global pool of random number generators. Defines
 * MajiqGufuncs::BetaMixture::global_rng_pool
 *
 * Copyright 2021 <University of Pennsylvania>
 */

#ifndef MAJIQGUFUNCS_BETAMIXTURE_GLOBALRNGPOOL_HPP
#define MAJIQGUFUNCS_BETAMIXTURE_GLOBALRNGPOOL_HPP

#include <majiqinclude/ResourcePool.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace MajiqGufuncs {
namespace BetaMixture {

using RNGT = boost::random::mt19937;
using RNGPoolT = MajiqInclude::ResourcePool<RNGT>;
using RNGPtrT = typename RNGPoolT::ExternalPtrT;
static RNGPoolT global_rng_pool{};

}  // namespace BetaMixture
}  // namespace majiqGufuncs

#endif  // MAJIQGUFUNCS_BETAMIXTURE_GLOBALRNGPOOL_HPP
