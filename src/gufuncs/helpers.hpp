/**
 * helpers.hpp
 *
 * MajiqGufuncs::detail helper functions
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_HELPERS_HPP
#define MAJIQGUFUNCS_HELPERS_HPP

#include <numpy/ndarraytypes.h>

namespace MajiqGufuncs {
namespace detail {
/**
 * Get x[idx] given type and stride of array
 */
template <typename T>
inline T& get_value(char* x_ptr, npy_intp idx, npy_intp stride) {
  return *reinterpret_cast<T*>(x_ptr + stride * idx);
}
}  // namespace detail
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_HELPERS_HPP
