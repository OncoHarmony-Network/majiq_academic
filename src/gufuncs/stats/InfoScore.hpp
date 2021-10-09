/**
 * InfoScore.hpp
 *
 * Statistics relative to InfoScore (see ScoreGenes)
 *
 * Copyright 2020 <University of Pennsylvania
 */

#ifndef MAJIQGUFUNCS_INFOSCORE_HPP
#define MAJIQGUFUNCS_INFOSCORE_HPP

#include <numpy/ndarraytypes.h>

#include <limits>

#include <gufuncs/CoreIt.hpp>
#include <majiqinclude/stats/InfoScore.hpp>

namespace MajiqGufuncs {
namespace InfoScore {

static char name[] = "infoscore";
constexpr int nin = 3;
constexpr int nout = 1;
static char signature[] = "(n),(n),(n)->()";
static char doc[] = R"pbdoc(
Compute p-values for InfoScore test on input data

Compute p-values for InfoScore test on input data omitting nan values.

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
  using MajiqGufuncs::detail::CoreIt;
  auto x = CoreIt<RealT>::begin(args[0], outer_stride[0]);
  auto sortx = CoreIt<npy_intp>::begin(args[1], outer_stride[1]);
  auto labels = CoreIt<bool>::begin(args[2], outer_stride[2]);
  auto out = CoreIt<RealT>::begin(args[3], outer_stride[3]);

  if (dim_core < 1) {
    // no samples to process, so must be nan
    out.fill(dim_broadcast, std::numeric_limits<RealT>::quiet_NaN());
    return;
  }

  // create cache of previous tests for MannWhitney
  using MajiqInclude::InfoScore::InfoScoreCache;
  InfoScoreCache tester;
  // outer loop on broadcasted variables
  for (npy_intp i = 0; i < dim_broadcast; ++i, ++x, ++sortx, ++labels, ++out) {
    using MajiqInclude::InfoScore::Test;
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

}  // namespace InfoScore
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_INFOSCORE_HPP
