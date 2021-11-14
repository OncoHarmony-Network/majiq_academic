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

#include <gufuncs/CoreIt.hpp>

namespace MajiqGufuncs {
namespace OffsetOr {

// assumes d_xout > 1, d_offsets >= 2
template <typename ItX, typename ItOffsets, typename ItOut>
inline void Inner(
    ItX x, ItOffsets offsets, ItOut out,
    const npy_intp d_xout, const npy_intp d_offsets) {
  npy_intp i = 0;  // index over x/out
  // before first offset, input to output
  for (; i < std::min(*offsets, d_xout); ++i, ++x, ++out) {
    *out = *x;
  }
  if (i == d_xout) { return; }
  // iterate to next offset for each offset
  ++offsets;
  for (npy_intp j = 1; j < d_offsets; ++j, ++offsets) {
    const auto& next_offset = std::min(*offsets, d_xout);
    if (i >= next_offset) { continue; }  // offsets that come before irrelevant
    // accumulate OR over x
    npy_bool acc_x{NPY_FALSE};
    for (npy_intp i_x = i; i_x < next_offset; ++i_x, ++x) {
      if (*x) {
        acc_x = NPY_TRUE;
        // get x ready for next offset, no need for next_offset - i_x iterations
        x += next_offset - i_x;
        break;
      }
    }
    // set out for these offsets to the accumulated value from x
    for (; i < next_offset; ++i, ++out) {
      *out = acc_x;
    }
    if (i >= d_xout) { return; }
  }
  // after last offset, input to output
  for (; i < d_xout; ++i, ++x, ++out) {
    *out = *x;
  }
  return;
}

static char name[] = "offset_logical_or";
constexpr int nin = 2;
constexpr int nout = 1;
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

static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp* outer_stride = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_xout = dimensions[0];
  const npy_intp dim_offsets = dimensions[1];
  // inner strides
  const npy_intp inner_stride_x = steps[0];
  const npy_intp inner_stride_offsets = steps[1];
  const npy_intp inner_stride_out = steps[2];

  // pointers to data
  auto x = detail::CoreIt<npy_bool>::begin(args[0], outer_stride[0]);
  auto offsets = detail::CoreIt<npy_intp>::begin(args[1], outer_stride[1]);
  auto out = detail::CoreIt<npy_bool>::begin(args[2], outer_stride[2]);

  // trivial cases
  if (dim_xout < 1 || dim_broadcast < 1) {
    return;
  } else if (dim_xout == 1) {
    // there is only one value, so we can just directly copy
    if (dim_broadcast > 1) {
      std::copy(x, x + dim_broadcast, out);
    } else {
      *out = *x;
    }
    return;
  } else if (dim_offsets < 2) {
    // no actual offsets, so we can just directly copy
    // note that dim_xout > 1
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++out) {
      auto x_x = x.with_stride(inner_stride_x);
      auto out_out = out.with_stride(inner_stride_out);
      std::copy(x_x, x_x + dim_xout, out_out);
    }
    return;
  }

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++offsets, ++out) {
    Inner(
        x.with_stride(inner_stride_x),
        offsets.with_stride(inner_stride_offsets),
        out.with_stride(inner_stride_out),
        dim_xout, dim_offsets);
  }
}

constexpr int ntypes = 1;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer),
};
static char types[ntypes * (nin + nout)] = { NPY_BOOL, NPY_INTP, NPY_BOOL };

}  // namespace OffsetOr
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_OFFSETOR_HPP
