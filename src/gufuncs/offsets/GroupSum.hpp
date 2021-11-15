/**
 * GroupSum.hpp
 *
 * Implementation to take sum over group indexes
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_GROUPSUM_HPP
#define MAJIQGUFUNCS_GROUPSUM_HPP

#include <numpy/ndarraytypes.h>

#include <algorithm>

#include <gufuncs/CoreIt.hpp>

namespace MajiqGufuncs {
namespace GroupSum {

template <typename ItX, typename ItG, typename ItOut,
         typename GroupT = typename std::iterator_traits<ItG>::value_type,
         typename OutT = typename std::iterator_traits<ItOut>::value_type>
inline void Inner(
    ItX x, ItG gx, const npy_intp n_x, ItOut out, const npy_intp n_out) {
  // set out to zero
  if (n_out < 1) {
    return;
  } else if (n_out == 1) {
    *out = OutT{0};
  } else {
    std::fill(out, out + n_out, OutT{0});
  }
  // fill in values
  for (npy_intp i = 0; i < n_x; ++i, ++x, ++gx) {
    if (0 <= *gx && *gx < static_cast<GroupT>(n_out)) {
      out[*gx] += *x;
    }
  }
  return;
}

static char name[] = "groupsum";
constexpr int nin = 3;
constexpr int nout = 1;
static char signature[] = "(n),(n),(g)->(g)";
static char doc[] = R"pbdoc(
per-group sum over x for unique groups

Parameters
----------
x: array
  Values to sum over
gx: array[int]
  Groups which values in x belong to
dummy: array[float]
  dummy array. shape of core axis determines the groups {0, ..., g - 1} which
  will have their sum computed

Returns
-------
array
  Sums of values of x with shared values of gx corresponding to indices
)pbdoc";

template <typename T, typename GroupT, typename OutT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp* outer_stride = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_x = dimensions[0];
  const npy_intp dim_out = dimensions[1];
  // inner strides
  const npy_intp str_x = steps[0];
  const npy_intp str_gx = steps[1];
  // ignore str_dummy = steps[2]
  const npy_intp str_out = steps[3];

  // pointers to data
  using MajiqGufuncs::detail::CoreIt;
  auto x = CoreIt<T>::begin(args[0], outer_stride[0]);
  auto gx = CoreIt<GroupT>::begin(args[1], outer_stride[1]);
  // skip dummy ~ args[2], outer_stride[2]
  auto out = CoreIt<OutT>::begin(args[3], outer_stride[3]);

  // trivial cases
  if (dim_broadcast < 1 || dim_out < 1) {
    return;
  } else if (dim_x < 1) {
    for (npy_intp i = 0; i < dim_broadcast; ++i, ++out) {
      if (dim_out == 1) {
        *out = OutT{0};
      } else {
        auto out_out = out.with_stride(str_out);
        std::fill(out_out, out_out + dim_out, OutT{0});
      }
    }
    return;
  }

  // outer loop over broadcast dimension
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++gx, ++out) {
    Inner(x.with_stride(str_x), gx.with_stride(str_gx), dim_x,
        out.with_stride(str_out), dim_out);
  }
  return;
}

constexpr int ntypes = 8;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_bool, npy_int64, npy_int64>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_bool, npy_uint64, npy_int64>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_int64, npy_int64, npy_int64>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_int64, npy_uint64, npy_int64>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, npy_int64, npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float, npy_uint64, npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, npy_int64, npy_double>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double, npy_uint64, npy_double>),
};
static char types[ntypes * (nin + nout)] = {
  NPY_BOOL, NPY_INT64, NPY_DOUBLE, NPY_INT64,
  NPY_BOOL, NPY_UINT64, NPY_DOUBLE, NPY_INT64,
  NPY_INT64, NPY_INT64, NPY_DOUBLE, NPY_INT64,
  NPY_INT64, NPY_UINT64, NPY_DOUBLE, NPY_INT64,
  NPY_FLOAT, NPY_INT64, NPY_DOUBLE, NPY_FLOAT,
  NPY_FLOAT, NPY_UINT64, NPY_DOUBLE, NPY_FLOAT,
  NPY_DOUBLE, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE,
  NPY_DOUBLE, NPY_UINT64, NPY_DOUBLE, NPY_DOUBLE,
};

}  // namespace GroupSum
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_GROUPSUM_HPP
