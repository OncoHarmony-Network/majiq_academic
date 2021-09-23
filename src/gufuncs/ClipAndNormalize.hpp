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

#include "helpers.hpp"


namespace MajiqGufuncs {
namespace ClipAndNormalize {
/**
 * Inner loop. strict indicates if 0 / 0 = NaN (True) or 0 (False)
 */
template <typename RealT, bool strict>
inline void Inner(
    char* x, char* offsets, char* out,
    const npy_intp d_xout, const npy_intp d_offsets,
    const npy_intp s_x, const npy_intp s_offsets, const npy_intp s_out) {
  using detail::get_value;
  // indexes into x, out
  npy_intp i_x = 0;
  // get offset to start with
  npy_intp next_offset = d_xout;
  if (d_offsets > 0) {
    next_offset = std::min(
        next_offset, get_value<npy_intp>(offsets, 0, s_offsets));
  }
  // before first offset, we just clip values to be non-negative
  for (; i_x < next_offset; ++i_x) {
    const RealT& xi = get_value<RealT>(x, i_x, s_x);
    RealT& outi = get_value<RealT>(out, i_x, s_out);
    if (npy_isnan(xi)) {
      outi = xi;
    } else if (xi > 0) {
      outi = RealT{1};
    } else {
      if constexpr(strict) {
        outi = std::numeric_limits<RealT>::quiet_NaN();
      } else {
        outi = RealT{0};
      }
    }
  }  // loop over x, out before first offset
  // between offsets, we will get sum of values
  for (npy_intp i_offsets = 1; i_offsets < d_offsets; ++i_offsets) {
    if (i_x == d_xout) {
      // there are no more values to be checked
      return;
    }
    next_offset = std::max(
        next_offset, std::min(
          d_xout, get_value<npy_intp>(offsets, i_offsets, s_offsets)));
    // we will loop over out separately, track where it started
    npy_intp i_out = i_x;
    // accumulate x between offsets
    RealT acc_x{0};
    for (; i_x < next_offset; ++i_x) {
      const RealT& xi = get_value<RealT>(x, i_x, s_x);
      RealT& outi = get_value<RealT>(out, i_x, s_out);
      if (npy_isnan(xi)) {
        acc_x = xi;
        i_x = next_offset;  // make it as if i_x reached end of loop
        break;
      } else if (xi > 0) {
        outi = xi;
        acc_x += outi;
      } else {
        outi = RealT{0};
      }
    }  // done accumulating x between offsets
    // update out
    if (npy_isnan(acc_x)) {
      for (; i_out < next_offset; ++i_out) {
        get_value<RealT>(out, i_out, s_out) = acc_x;
      }
    } else if (acc_x > 0) {
      for (; i_out < next_offset; ++i_out) {
        get_value<RealT>(out, i_out, s_out) /= acc_x;
      }
    } else {
      if constexpr(strict) {
        for (; i_out < next_offset; ++i_out) {
          get_value<RealT>(out, i_out, s_out)
            = std::numeric_limits<RealT>::quiet_NaN();
        }
      } else {
        for (; i_out < next_offset; ++i_out) {
          get_value<RealT>(out, i_out, s_out) = RealT{0};
        }
      }
    }
  }  // loop over x, out between offsets
  for (; i_x < d_xout; ++i_x) {
    const RealT& xi = get_value<RealT>(x, i_x, s_x);
    RealT& outi = get_value<RealT>(out, i_x, s_out);
    if (npy_isnan(xi)) {
      outi = xi;
    } else if (xi > 0) {
      outi = RealT{1};
    } else {
      if constexpr(strict) {
        outi = std::numeric_limits<RealT>::quiet_NaN();
      } else {
        outi = RealT{0};
      }
    }
  }  // loop over x, out after last offset
  return;
}

// implement clip_and_normalize(x: np.ndarray, offsets: np.ndarray)
template <typename RealT, bool strict>
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
    Inner<RealT, strict>(
        x_ptr, offsets_ptr, out_ptr,
        dim_xout, dim_offsets,
        inner_stride_x, inner_stride_offsets, inner_stride_out);
  }
}

static char name[] = "clip_and_normalize";
static char strictname[] = "clip_and_normalize_strict";
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
constexpr int ntypes = 2;
constexpr int nin = 2;
constexpr int nout = 1;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, false>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, false>)
};
PyUFuncGenericFunction strictfuncs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, true>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, true>)
};
static char types[
  ntypes * (nin + nout)
] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_INTP, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_INTP, NPY_DOUBLE,
};


}  // namespace ClipAndNormalize
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_CLIPANDNORMALIZE_HPP
