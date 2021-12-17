/**
 * globalRNG.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PYBIND_GLOBALRNG_HPP
#define MAJIQ_PYBIND_GLOBALRNG_HPP

#include "../MajiqTypes.hpp"
#include <majiqinclude/ResourcePool.hpp>

namespace majiq {
namespace bindings {

static MajiqInclude::ResourcePool<rng_t> global_rng_pool{};

inline void init_RNGfunctions(pybind11::module_& m) {
  m.def("rng_seed",
      [](int64_t x) { global_rng_pool.seed(x); },
      "Set seed for pool of RNGs in new_majiq.internals",
      pybind11::arg("seed"));
  m.def("rng_resize",
      [](int64_t n) { global_rng_pool.resize(n); },
      "Resize pool of RNGs for at least n simultaneous threads",
      pybind11::arg("n"));
  return;
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_GLOBALRNG_HPP
