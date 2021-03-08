/**
 * pybind_utils.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PYBIND_UTILS_HPP
#define MAJIQ_PYBIND_UTILS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>
#include <sstream>
#include <string>
#include <memory>
#include <array>
#include <stdexcept>

#include "internals/MajiqTypes.hpp"


// groups for netcdf
constexpr char CONTIGS_NC_GROUP[] = "contigs";
constexpr char GENES_NC_GROUP[] = "genes";
constexpr char EXONS_NC_GROUP[] = "exons";
constexpr char JUNCTIONS_NC_GROUP[] = "junctions";
constexpr char INTRONS_NC_GROUP[] = "introns";
constexpr char SJ_JUNCTIONS_NC_GROUP[] = "sj_junctions";
constexpr char SJ_JUNCTIONS_RAW_NC_GROUP[] = "sj_junctions_raw";
constexpr char CONTIG_INTRONS_NC_GROUP[] = "contig_introns";

constexpr bool DEFAULT_BUILD_PROCESS_IR = true;
constexpr majiq::junction_ct_t DEFAULT_BUILD_MINREADS = 3;
constexpr majiq::junction_ct_t DEFAULT_BUILD_MINDENOVO = 5;
constexpr majiq::junction_pos_t DEFAULT_BUILD_MINPOS = 2;
constexpr majiq::real_t DEFAULT_BUILD_MAX_PCTBINS = 0.6;
constexpr majiq::real_t DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY = 0.5;
constexpr majiq::real_t DEFAULT_BUILD_MATCH_INTRON_PROBABILITY = 0.95;
constexpr majiq::real_t DEFAULT_BUILD_MINEXPERIMENTS = 0.5;
constexpr bool DEFAULT_BUILD_DENOVO_JUNCTIONS = true;
constexpr bool DEFAULT_BUILD_DENOVO_IR = true;
constexpr bool DEFAULT_BUILD_KEEP_ANNOTATED_IR = false;
constexpr size_t DEFAULT_BUILD_NUM_BOOTSTRAPS = 30;
constexpr majiq::real_t DEFAULT_BUILD_STACK_PVALUE = 1e-7;
constexpr majiq::ExperimentStrandness DEFAULT_BAM_STRANDNESS
  = majiq::ExperimentStrandness::NONE;
constexpr int DEFAULT_BAM_NTHREADS = 1;
constexpr majiq::real_t DEFAULT_BUILD_SIMPL_MINPSI = 0.01;
constexpr majiq::real_t DEFAULT_BUILD_SIMPL_MINREADS_ANNOTATED_JUNCTION = 0;
constexpr majiq::real_t DEFAULT_BUILD_SIMPL_MINREADS_DENOVO_JUNCTION = 0;
constexpr majiq::real_t DEFAULT_BUILD_SIMPL_MINREADS_INTRON = 0;



namespace majiq_pybind {
namespace py = pybind11;

/*
 * Create read-only array view into vector with offset (i.e. for struct member)
 */
template <class OutputT, class InputT>
inline py::array_t<OutputT> ArrayFromVectorAndOffset(
    const std::vector<InputT>& src,
    size_t offset,
    py::object& handle) {
  // pointer to first element of src after adding specified memory offset
  const OutputT* first = reinterpret_cast<const OutputT*>(
      reinterpret_cast<const char*>(src.data()) + offset);
  // construct array
  py::array_t<OutputT> result = py::array_t(
      // shape
      {src.size()},
      // strides
      {sizeof(InputT)},
      // pointer to first element
      first,
      // object to handle array
      handle);
  // set array to readonly -- discourage Python from editing it
  py::detail::array_proxy(result.ptr())->flags
    &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  // return resulting array
  return result;
}

/*
 * Create read-only array view into vector with offset (i.e. for struct member)
 */
template <class OutputT>
inline py::array_t<OutputT> ArrayFromOffsetsVector(
    const std::vector<OutputT>& src,
    bool is_start,
    py::object& handle) {
  // pointer to first element of src after adding offset for start or end
  const OutputT* first = src.data() + (is_start ? 0 : 1);
  // construct array
  py::array_t<OutputT> result = py::array_t(
      // shape
      {src.size() - 1},
      // strides
      {sizeof(OutputT)},
      // pointer to first element
      first,
      // object to handle array
      handle);
  // set array to readonly -- discourage Python from editing it
  py::detail::array_proxy(result.ptr())->flags
    &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  // return resulting array
  return result;
}

/**
 * Extract attrs from features to build xarray Dataset
 */
inline py::object XarrayDatasetFromObject(
    py::object features,
    py::str idx_name,
    std::initializer_list<py::str> attrs,
    std::initializer_list<py::str> data_attrs = {}) {
  py::function np_arange = py::module_::import("numpy").attr("arange");
  py::module_ xr = py::module_::import("xarray");
  py::function xr_Dataset = xr.attr("Dataset");
  py::function xr_DataArray = xr.attr("DataArray");
  // define index with appropriate length
  size_t inferred_length = 0;
  for (auto x : attrs) {
    size_t x_len = features.attr(x).attr("__len__")().cast<size_t>();
    if (inferred_length > 0 && x_len != inferred_length) {
      std::ostringstream oss;
      oss << "Dataset has arrays with mismatched lengths, including "
        << inferred_length << " and " << x_len;
      throw std::runtime_error(oss.str());
    }
    inferred_length = x_len;
  }
  for (auto x : data_attrs) {
    size_t x_len = features.attr(x).attr("__len__")().cast<size_t>();
    if (inferred_length > 0 && x_len != inferred_length) {
      std::ostringstream oss;
      oss << "Dataset has arrays with mismatched lengths, including "
        << inferred_length << " and " << x_len;
      throw std::runtime_error(oss.str());
    }
    inferred_length = x_len;
  }
  inferred_length
    = inferred_length > 0
    ? inferred_length : features.attr("__len__")().cast<size_t>();
  // coordinates for dataset in dictionary
  py::dict coordinates;
  coordinates[idx_name] = xr_DataArray(
      np_arange(inferred_length), py::arg("dims") = idx_name);
  // extract attributes
  for (auto x : attrs) {
    coordinates[x] = xr_DataArray(features.attr(x), py::arg("dims") = idx_name);
  }
  py::dict data_vars;
  for (auto x : data_attrs) {
    data_vars[x] = xr_DataArray(features.attr(x), py::arg("dims") = idx_name);
  }
  // define the index dimension
  return xr_Dataset(
      py::arg("data_vars") = data_vars,
      py::arg("coords") = coordinates);
}

inline py::object OpenXarrayDataset(py::str netcdf_path, py::str group) {
  return py::module_::import("xarray").attr("open_dataset")(
      netcdf_path, py::arg("group") = group);
}

}  // namespace majiq_pybind


#endif  // MAJIQ_PYBIND_UTILS_HPP
