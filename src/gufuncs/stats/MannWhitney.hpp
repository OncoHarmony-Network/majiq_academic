/**
 * MannWhitney.hpp
 *
 * Implementation of Mann-Whitney U (two-sample Wilcoxon) two-sided test for
 * use in MAJIQ.
 *
 * Copyright 2020 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_MANNWHITNEY_HPP
#define MAJIQGUFUNCS_MANNWHITNEY_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/stats/MannWhitney.hpp>


namespace MajiqGufuncs {
namespace MannWhitney {

static char name[] = "mannwhitneyu";
constexpr int nin = 3;
constexpr int nout = 1;
static char signature[] = "(n),(n),(n)->()";
static char doc[] = R"pbdoc(
Compute p-values for Mann-Whitney U test on input data

Compute p-values for Mann-Whitney U test on input data, using a two-sided
alternative hypothesis and omitting nan values.

Parameters
----------
x: array[float]
    test over observations in last axis
sortx: array[int]
    Previously computed values of np.argsort(x, axis=-1) that will not be
    checked
labels: array[bool]
    test over labels in last axis

Returns
-------
array[float]
    broadcast p-values for observations/labels. Invalid tests are nan
)pbdoc";

template <typename RealT>
static void Outer(
    char** args, npy_intp* dimensions, npy_intp* steps, void* data) {
  // outer loop dimensions and index
  const npy_intp dim_broadcast = *dimensions++;
  // strides on each variable for outer loop
  const npy_intp* outer_stride = steps;
  steps += nin + nout;

  // core dimensions
  const npy_intp dim_core = dimensions[0];
  // inner strides
  const npy_intp inner_stride_x = steps[0];
  const npy_intp inner_stride_sortx = steps[1];
  const npy_intp inner_stride_labels = steps[2];

  // pointers to data
  auto x = detail::CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto sortx = detail::CoreIt<npy_intp>::begin(args[1], outer_stride[1]);
  auto labels = detail::CoreIt<bool>::begin(args[2], outer_stride[2]);
  auto out = detail::CoreIt<RealT>::begin(args[3], outer_stride[3]);

  if (dim_core < 1) {
    // no samples to process, so must be nan
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }

  // create cache of previous tests for MannWhitney
  using MajiqInclude::MannWhitney::MannWhitneyCache;
  MannWhitneyCache tester;
  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++sortx, ++labels, ++out) {
  using MajiqInclude::MannWhitney::Test;
    *out = Test(
        tester,
        x.with_stride(inner_stride_x),
        sortx.with_stride(inner_stride_sortx),
        labels.with_stride(inner_stride_labels),
        dim_core);
  }
  return;
}

constexpr int ntypes = 2;
PyUFuncGenericFunction funcs[ntypes] = {
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_float>),
  reinterpret_cast<PyUFuncGenericFunction>(&Outer<npy_double>)
};
static char types[ntypes * (nin + nout)] = {
  // for use with npy_float func
  NPY_FLOAT, NPY_INTP, NPY_BOOL, NPY_DOUBLE,
  // for use with npy_double func
  NPY_DOUBLE, NPY_INTP, NPY_BOOL, NPY_DOUBLE,
};

}  // namespace MannWhitney
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_MANNWHITNEY_HPP
