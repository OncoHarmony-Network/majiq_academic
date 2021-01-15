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
#include <memory>
#include <array>

#include "internals/Contigs.hpp"
#include "internals/Genes.hpp"


namespace majiq_pybind {
namespace py = pybind11;

constexpr char CONTIGS_NC_GROUP[] = "contigs";
constexpr char GENES_NC_GROUP[] = "genes";

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
    std::initializer_list<py::str> attrs) {
  using namespace py::literals;
  py::function np_arange = py::module_::import("numpy").attr("arange");
  py::module_ xr = py::module_::import("xarray");
  py::function xr_Dataset = xr.attr("Dataset");
  py::function xr_DataArray = xr.attr("DataArray");
  // coordinates for dataset in dictionary
  py::dict coordinates;
  // define the index dimension
  coordinates[idx_name] = xr_DataArray(
      np_arange(features.attr("__len__")()), "dims"_a = idx_name);
  // extract remaining attributes
  for (auto x : attrs) {
    coordinates[x] = xr_DataArray(features.attr(x), "dims"_a = idx_name);
  }
  return xr_Dataset("coords"_a = coordinates);
}

inline py::object OpenXarrayDataset(py::str netcdf_path, py::str group) {
  using namespace py::literals;
  return py::module_::import("xarray").attr("open_dataset")(
      netcdf_path, "group"_a = group);
}

inline std::shared_ptr<majiq::Contigs> ContigsFromNetcdf(py::str netcdf_path) {
  auto xr_contigs = OpenXarrayDataset(netcdf_path, py::str(CONTIGS_NC_GROUP));
  auto result = std::make_shared<majiq::Contigs>();
  py::list seqids = xr_contigs
    .attr("__getitem__")("seqid").attr("values").attr("tolist")();
  for (auto seqid : seqids) {
    result->add(seqid.cast<majiq::seqid_t>());
  }
  return result;
}

inline std::shared_ptr<majiq::Genes> GenesFromNetcdf(
    std::shared_ptr<majiq::Contigs> contigs, py::str netcdf_path) {
  using majiq::position_t;
  using majiq::geneid_t;
  using majiq::genename_t;
  auto xr_genes = OpenXarrayDataset(netcdf_path, py::str(GENES_NC_GROUP));
  auto get_array = [&xr_genes](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_genes.attr("__getitem__")(key));
  };
  py::array_t<size_t> _contig_idx = get_array("contig_idx");
  auto contig_idx = _contig_idx.unchecked<1>();
  py::array_t<position_t> _start = get_array("start");
  auto start = _start.unchecked<1>();
  py::array_t<position_t> _end = get_array("end");
  auto end = _end.unchecked<1>();
  py::array_t<std::array<char, 1>> _strand = get_array("strand");
  auto strand = _strand.unchecked<1>();
  py::list geneid
    = xr_genes.attr("__getitem__")("gene_id").attr("values").attr("tolist")();
  py::list genename
    = xr_genes.attr("__getitem__")("gene_name").attr("values").attr("tolist")();
  std::vector<majiq::Gene> gene_vec{};
  gene_vec.reserve(geneid.size());
  for (size_t i = 0; i < geneid.size(); ++i) {
    gene_vec.push_back(majiq::Gene{
        majiq::KnownContig{contig_idx(i), contigs},
        majiq::ClosedInterval{start(i), end(i)},
        static_cast<majiq::GeneStrandness>(strand(i)[0]),
        geneid[i].cast<majiq::geneid_t>(),
        genename[i].cast<majiq::genename_t>()});
  }
  return std::make_shared<majiq::Genes>(gene_vec);
}

}  // namespace majiq_pybind


#endif  // MAJIQ_PYBIND_UTILS_HPP
