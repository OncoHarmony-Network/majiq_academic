/**
 * OffsetOr.hpp
 *
 * Implementation to summarize logical_or between offsets without changing
 * dimension
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_OFFSETOR_HPP
#define MAJIQGUFUNCS_OFFSETOR_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>

#include <gufuncs/helpers.hpp>

namespace MajiqGufuncs {
namespace OffsetOr {

template <typename ItX, typename ItOffsets, typename ItOut>
inline void Inner(
    ItX x, ItOffsets offsets, ItOut out,
    const npy_intp d_xout, const npy_intp d_offsets) {
  // indexes into x, out
  npy_intp i_x = 0;
  // get offset to start with
  npy_intp next_offset = d_xout;
  if (d_offsets > 0) {
    next_offset = std::min(next_offset, offsets[0]);
  }
  // before first offset, no groups --> identify
  for (; i_x < next_offset; ++i_x) {
    out[i_x] = x[i_x];
  }  // loop over x, out before first offset
  // between offsets, we will get sum of values
  for (npy_intp i_offsets = 1; i_offsets < d_offsets; ++i_offsets) {
    if (i_x == d_xout) {
      return;  // no more values to check/update
    }
    next_offset = std::max(next_offset, std::min(d_xout, offsets[i_offsets]));
    // loop to fill out happens separately than first pass on x
    npy_intp i_out = i_x;
    // accumulate x between offsets
    npy_bool acc_x{NPY_FALSE};
    for (; i_x < next_offset; ++i_x) {
      if (x[i_x]) {
        acc_x = NPY_TRUE;
        i_x = next_offset;
        break;
      }
    }  // done accumulating x between offsets
    // update out
    for (; i_out < next_offset; ++i_out) {
      out[i_out] = acc_x;
    }
  }  // done looping between offsets
  // after last offset, no groups --> identify
  for (; i_x < d_xout; ++i_x) {
    out[i_x] = x[i_x];
  }  // loop over x, out after last offset
  return;
}

static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp stride_x = *steps++;
  const npy_intp stride_offsets = *steps++;
  const npy_intp stride_out = *steps++;
  // core dimensions
  const npy_intp dim_xout = dimensions[0];
  const npy_intp dim_offsets = dimensions[1];
  // inner strides
  const npy_intp inner_stride_x = steps[0];
  const npy_intp inner_stride_offsets = steps[1];
  const npy_intp inner_stride_out = steps[2];
  // pointers to data
  char* x_ptr = args[0];
  char* offsets_ptr = args[1];
  char* out_ptr = args[2];

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i,
      x_ptr += stride_x, offsets_ptr += stride_offsets, out_ptr += stride_out) {
    Inner(
        detail::CoreIt<npy_bool>::begin(x_ptr, inner_stride_x),
        detail::CoreIt<npy_intp>::begin(offsets_ptr, inner_stride_offsets),
        detail::CoreIt<npy_bool>::begin(out_ptr, inner_stride_out),
        dim_xout, dim_offsets);
  }
}

static char name[] = "offset_logical_or";
static char signature[] = "(n),(k)->(n)";
static char doc[] = R"pbdoc(
logical or over groups defined by offsets

Parameters
----------
x1: array_like
    Booleans to be summarized
x2: array_like
    Offsets to be used
)pbdoc";
constexpr int ntypes = 1;
constexpr int nin = 2;
constexpr int nout = 1;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer),
};
static char types[ntypes * (nin + nout)] = { NPY_BOOL, NPY_INTP, NPY_BOOL };

}  // namespace OffsetOr
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_OFFSETOR_HPP
