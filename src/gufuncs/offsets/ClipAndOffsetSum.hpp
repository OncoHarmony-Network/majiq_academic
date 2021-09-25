/**
 * ClipAndOffsetSum.hpp
 *
 * Implementation to clip and normalize input vector between offsets
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_CLIPANDOFFSETSUM_HPP
#define MAJIQGUFUNCS_CLIPANDOFFSETSUM_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>
#include <limits>

#include <gufuncs/helpers.hpp>


namespace MajiqGufuncs {
namespace ClipAndOffsetSum {

struct NoClip {
  // identity function
  template <typename RealT>
  RealT operator()(RealT x) { return x; }
};
struct Clip {
  template <typename RealT>
  RealT operator()(RealT x) { return std::max(RealT{0}, x); }
};

template <typename ClipT, typename ItX, typename ItOffsets, typename ItOut,
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
  // before first offset -- no groups
  for (; i_x < next_offset; ++i_x) {
    const auto& xi = x[i_x];
    if (npy_isnan(xi)) {
      out[i_x] = xi;
    } else {
      out[i_x] = ClipT{}(xi);
    }
  }  // reached first offset
  // loop within each group defined by offsets
  for (npy_intp i_offsets = 1; i_offsets < d_offsets; ++i_offsets) {
    if (i_x == d_xout) {
      return;  // no more values to check/add
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
      } else {
        acc_x += ClipT{}(xi);
      }
    }  // accumulated x between offsets
    // update out
    for (; i_out < next_offset; ++i_out) {
      out[i_out] = acc_x;
    }
  }  // done looping over offsets
  // after last offset -- no groups
  for (; i_x < d_xout; ++i_x) {
    const auto& xi = x[i_x];
    if (npy_isnan(xi)) {
      out[i_x] = xi;
    } else {
      out[i_x] = ClipT{}(xi);
    }
  }  // done looping over x, out
  return;
}

// implement clip_and_normalize(x: np.ndarray, offsets: np.ndarray)
template <typename RealT, typename ClipT>
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
    Inner<ClipT>(
        detail::CoreIt<RealT>::begin(x_ptr, inner_stride_x),
        detail::CoreIt<npy_intp>::begin(offsets_ptr, inner_stride_offsets),
        detail::CoreIt<RealT>::begin(out_ptr, inner_stride_out),
        dim_xout, dim_offsets);
  }
}

static char name[] = "clip_and_offsetsum";
static char noclipname[] = "offsetsum";
static char signature[] = "(n),(k)->(n)";
static char doc[] = R"pbdoc(
Sum of positive values between offsets

Signature: (n),(k)->(n)

Computes sum of positive values between adjacent cummax(offsets). NaN values
are propagated. Values before and after offsets are treated as their own
groups.

Notes
-----
This can be seen as a group sum where groups are defined as slices between
offsets.

Parameters
----------
x1: array_like
    Values to be added
x2: array_like
    Offsets to be used
)pbdoc";
static char noclipdoc[] = R"pbdoc(
Sum of values between offsets

Signature: (n),(k)->(n)

Computes sum of values between adjacent cummax(offsets). NaN values are
propagated. Values before and after offsets are treated as their own groups.

Notes
-----
This can be seen as a group sum where groups are defined as slices between
offsets.

Parameters
----------
x1: array_like
    Values to be added
x2: array_like
    Offsets to be used
)pbdoc";
constexpr int ntypes = 2;
constexpr int nin = 2;
constexpr int nout = 1;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, Clip>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, Clip>)
};
PyUFuncGenericFunction noclipfuncs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, NoClip>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, NoClip>)
};
static char types[
  ntypes * (nin + nout)
] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_INTP, NPY_FLOAT,
  // for use with npy_double func
  NPY_DOUBLE, NPY_INTP, NPY_DOUBLE,
};


}  // namespace ClipAndOffsetSum
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_CLIPANDOFFSETSUM_HPP
