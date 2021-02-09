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

#include "internals/Contigs.hpp"


// groups for netcdf
constexpr char CONTIGS_NC_GROUP[] = "contigs";
constexpr char GENES_NC_GROUP[] = "genes";
constexpr char EXONS_NC_GROUP[] = "exons";
constexpr char JUNCTIONS_NC_GROUP[] = "junctions";
constexpr char INTRONS_NC_GROUP[] = "introns";
constexpr char SJ_JUNCTIONS_NC_GROUP[] = "sj_junctions";
constexpr char SJ_JUNCTIONS_RAW_NC_GROUP[] = "sj_junctions_raw";


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
  using namespace py::literals;
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
      np_arange(inferred_length), "dims"_a = idx_name);
  // extract attributes
  for (auto x : attrs) {
    coordinates[x] = xr_DataArray(features.attr(x), "dims"_a = idx_name);
  }
  py::dict data_vars;
  for (auto x : data_attrs) {
    data_vars[x] = xr_DataArray(features.attr(x), "dims"_a = idx_name);
  }
  // define the index dimension
  return xr_Dataset("data_vars"_a = data_vars, "coords"_a = coordinates);
}

inline py::object OpenXarrayDataset(py::str netcdf_path, py::str group) {
  using namespace py::literals;
  return py::module_::import("xarray").attr("open_dataset")(
      netcdf_path, "group"_a = group);
}

}  // namespace majiq_pybind


#endif  // MAJIQ_PYBIND_UTILS_HPP
