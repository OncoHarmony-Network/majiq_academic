/**
 * ClipAndNormalize.hpp
 *
 * Implementation to clip and normalize input vector between offsets
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_CLIPANDNORMALIZE_HPP
#define MAJIQGUFUNCS_CLIPANDNORMALIZE_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>
#include <limits>

#include <gufuncs/helpers.hpp>


namespace MajiqGufuncs {
namespace ClipAndNormalize {

struct NaNZeroCov {
  template <typename RealT>
  RealT operator()(RealT) {
    return std::numeric_limits<RealT>::quiet_NaN();
  }
};
struct ZeroZeroCov {
  template <typename RealT>
  RealT operator()(RealT) {
    return RealT{0};
  }
};

template <typename StrictT, typename ItX, typename ItOffsets, typename ItOut,
         typename RealT = typename std::iterator_traits<ItX>::value_type>
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
  // before first offset, no groups
  for (; i_x < next_offset; ++i_x) {
    const auto& xi = x[i_x];
    if (npy_isnan(xi)) {
      out[i_x] = xi;
    } else if (xi > 0) {
      out[i_x] = RealT{1};
    } else {
      out[i_x] = StrictT{}(RealT{});
    }
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
    RealT acc_x{0};
    for (; i_x < next_offset; ++i_x) {
      const auto& xi = x[i_x];
      if (npy_isnan(xi)) {
        acc_x = xi;
        i_x = next_offset;
        break;
      } else if (xi > 0) {
        out[i_x] = xi;
        acc_x += xi;
      } else {
        out[i_x] = RealT{0};
      }
    }  // done accumulating x between offsets
    // update out
    if (npy_isnan(acc_x)) {
      for (; i_out < next_offset; ++i_out) {
        out[i_out] = acc_x;
      }
    } else if (acc_x > 0) {
      for (; i_out < next_offset; ++i_out) {
        out[i_out] /= acc_x;
      }
    } else {
      for (; i_out < next_offset; ++i_out) {
        out[i_out] = StrictT{}(RealT{});
      }
    }
  }  // done looping between offsets
  // after last offset, no groups
  for (; i_x < d_xout; ++i_x) {
    const auto& xi = x[i_x];
    if (npy_isnan(xi)) {
      out[i_x] = xi;
    } else if (xi > 0) {
      out[i_x] = RealT{1};
    } else {
      out[i_x] = StrictT{}(RealT{});
    }
  }  // loop over x, out after last offset
  return;
}

static char name[] = "clip_and_normalize";
static char strictname[] = "clip_and_normalize_strict";
constexpr int nin = 2;
constexpr int nout = 1;
static char signature[] = "(n),(k)->(n)";
static char doc[] = R"pbdoc(
Clipped (non-negative), permissive normalized values for offsets (0 total -> 0)

Signature: (n),(k)->(n)

Computes sum of positive values between adjacent cummax(offsets). Return ratios
of clipped values to sum. NaN values are propagated within groups, and a sum
of zero leads to a zero result.

Notes
-----
This can be seen as a group sum where groups are defined as slices between
offsets.

Parameters
----------
x1: array_like
    Values to be clipped/normalized
x2: array_like
    Offsets to be used
)pbdoc";
static char strictdoc[] = R"pbdoc(
Clipped (non-negative), strict normalized values for offsets (0 total -> NaN)

Signature: (n),(k)->(n)

Computes sum of positive values between adjacent cummax(offsets). Return ratios
of clipped values to sum. NaN values are propagated within groups, and a sum
of zero leads to a NaN result.

Notes
-----
This can be seen as a group sum where groups are defined as slices between
offsets.

Parameters
----------
x1: array_like
    Values to be clipped/normalized
x2: array_like
    Offsets to be used
)pbdoc";

// implement clip_and_normalize(x: np.ndarray, offsets: np.ndarray)
template <typename RealT, typename StrictT>
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
  auto x = detail::CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto offsets = detail::CoreIt<npy_intp>::begin(args[1], outer_stride[1]);
  auto out = detail::CoreIt<RealT>::begin(args[2], outer_stride[2]);

  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++offsets, ++out) {
    Inner<StrictT>(
        x.with_stride(inner_stride_x),
        offsets.with_stride(inner_stride_offsets),
        out.with_stride(inner_stride_out),
        dim_xout, dim_offsets);
  }
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, ZeroZeroCov>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, ZeroZeroCov>)
};
PyUFuncGenericFunction strictfuncs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, NaNZeroCov>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, NaNZeroCov>)
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_INTP, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_INTP, NPY_DOUBLE,
};


}  // namespace ClipAndNormalize
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_CLIPANDNORMALIZE_HPP
