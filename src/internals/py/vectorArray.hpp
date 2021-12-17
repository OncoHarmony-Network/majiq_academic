/**
 * vectorArray.hpp
 *
 * Bindings for std::vector as pybind11::array_t
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_VECTORARRAY_HPP
#define MAJIQ_PYBIND_VECTORARRAY_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>

namespace majiq {
namespace bindings {

/*
 * Create read-only array view into vector with offset (i.e. for struct member)
 */
template <class OutputT, class InputT>
inline pybind11::array_t<OutputT> ArrayFromVectorAndOffset(
    const std::vector<InputT>& src,
    size_t offset,
    pybind11::object& handle) {
  // pointer to first element of src after adding specified memory offset
  const OutputT* first = reinterpret_cast<const OutputT*>(
      reinterpret_cast<const char*>(src.data()) + offset);
  // construct array
  pybind11::array_t<OutputT> result = pybind11::array_t(
      // shape
      {src.size()},
      // strides
      {sizeof(InputT)},
      // pointer to first element
      first,
      // object to handle array
      handle);
  // set array to readonly -- discourage Python from editing it
  pybind11::detail::array_proxy(result.ptr())->flags
    &= ~pybind11::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  // return resulting array
  return result;
}

/*
 * Create read-only array view into offsets vector for start vs end
 */
template <class OutputT>
inline pybind11::array_t<OutputT> ArrayFromOffsetsVector(
    const std::vector<OutputT>& src,
    bool is_start,
    pybind11::object& handle) {
  // pointer to first element of src after adding offset for start or end
  const OutputT* first = src.data() + (is_start ? 0 : 1);
  // construct array
  pybind11::array_t<OutputT> result = pybind11::array_t(
      // shape
      {src.size() - 1},
      // strides
      {sizeof(OutputT)},
      // pointer to first element
      first,
      // object to handle array
      handle);
  // set array to readonly -- discourage Python from editing it
  pybind11::detail::array_proxy(result.ptr())->flags
    &= ~pybind11::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  // return resulting array
  return result;
}


}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_VECTORARRAY_HPP
