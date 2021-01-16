/**
 * pySJJunctions.cpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <array>
#include <vector>
#include <memory>
#include <cstddef>

#include "pybind_utils.hpp"
#include "internals/SJJunctions.hpp"
#include "internals/MajiqTypes.hpp"


namespace py = pybind11;


std::shared_ptr<majiq::SJJunctions> JunctionsFromNetcdf(py::str netcdf_path) {
  using majiq::position_t;
  using majiq::junction_ct_t;
  using majiq::junction_pos_t;
  auto xr_junctions = majiq_pybind::OpenXarrayDataset(
      netcdf_path, py::str(SJ_JUNCTIONS_NC_GROUP));
  auto get_array = [&xr_junctions](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_junctions.attr("__getitem__")(key));
  };
  py::array_t<size_t> _contig_idx = get_array("contig_idx");
  auto contig_idx = _contig_idx.unchecked<1>();
  py::array_t<position_t> _start = get_array("start");
  auto start = _start.unchecked<1>();
  py::array_t<position_t> _end = get_array("end");
  auto end = _end.unchecked<1>();
  py::array_t<std::array<char, 1>> _strand = get_array("strand");
  auto strand = _strand.unchecked<1>();
  py::array_t<junction_pos_t> _numpos = get_array("numpos");
  auto numpos = _numpos.unchecked<1>();
  py::array_t<junction_ct_t> _numreads = get_array("numreads");
  auto numreads = _numreads.unchecked<1>();

  // load contigs
  auto contigs = majiq_pybind::ContigsFromNetcdf(netcdf_path);

  // fill in junctions
  std::vector<majiq::SJJunction> sj_vec(numpos.shape(0));
  for (size_t i = 0; i < sj_vec.size(); ++i) {
    sj_vec[i] = majiq::SJJunction{
      majiq::KnownContig{contig_idx(i), contigs},
      majiq::OpenInterval{start(i), end(i)},
      static_cast<majiq::GeneStrandness>(strand(i)[0]),
      majiq::ExperimentCounts{numreads(i), numpos(i)}
    };
  }
  return std::make_shared<majiq::SJJunctions>(std::move(sj_vec));
}

majiq::SJJunctionsPositions JunctionsRawFromNetcdf(py::str netcdf_path) {
  auto xr_raw = majiq_pybind::OpenXarrayDataset(
      netcdf_path, py::str(SJ_JUNCTIONS_RAW_NC_GROUP));
  auto get_array = [&xr_raw](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_raw.attr("__getitem__")(key));
  };
  using majiq::junction_ct_t;
  using majiq::junction_pos_t;
  py::array_t<junction_ct_t> _position_reads = get_array("position_reads");
  auto position_reads = _position_reads.unchecked<1>();
  py::array_t<junction_pos_t> _position = get_array("position");
  auto position = _position.unchecked<1>();
  py::array_t<size_t> _offsets = get_array("_offsets");
  auto offsets = _offsets.unchecked<1>();
  py::dict xr_raw_attrs = xr_raw.attr("attrs");

  auto num_positions = xr_raw_attrs["num_positions"].cast<junction_pos_t>();
  std::vector<size_t> offsets_vec(offsets.shape(0));
  for (size_t i = 0; i < offsets_vec.size(); ++i) {
    offsets_vec[i] = offsets(i);
  }
  std::vector<majiq::PositionReads> pr_vec(position_reads.shape(0));
  for (size_t i = 0; i < pr_vec.size(); ++i) {
    pr_vec[i] = majiq::PositionReads{position(i), position_reads(i)};
  }
  auto junctions = JunctionsFromNetcdf(netcdf_path);

  return majiq::SJJunctionsPositions{std::move(junctions), std::move(pr_vec),
                                     std::move(offsets_vec), num_positions};
}

void init_SJJunctionsModule(py::module_& m) {
  using majiq::position_t;
  using majiq::SJJunctions;
  using majiq::SJJunctionsPositions;
  using majiq::ExperimentStrandness;
  using majiq_pybind::ArrayFromOffsetsVector;
  using majiq_pybind::ArrayFromVectorAndOffset;

  auto pySJJunctions = py::class_<SJJunctions, std::shared_ptr<SJJunctions>>(
      m, "SJJunctions", "Summarized junction counts for an experiment");
  auto pySJJunctionsPositions = py::class_<SJJunctionsPositions>(
      m, "SJJunctionsPositions",
      "Summarized and per-position counts for an experiment");

  auto pyExperimentStrandness = py::enum_<ExperimentStrandness>(
      m, "ExperimentStrandness")
    .value("FORWARD", ExperimentStrandness::FORWARD)
    .value("REVERSE", ExperimentStrandness::REVERSE)
    .value("NONE", ExperimentStrandness::NONE);

  pySJJunctions
    // NOTE majiq::Contigs needs to be bound prior to this function being called
    .def_property_readonly("_contigs", &SJJunctions::contigs,
        "Underlying contigs corresponding to contig_idx")
    .def_property_readonly("contigs",
        [](py::object& self) { return self.attr("_contigs").attr("df")(); },
        "View underlying contigs as xarray Dataset")
    .def_property_readonly("contig_idx",
        [](py::object& sj_obj) -> py::array_t<size_t> {
        SJJunctions& sj = sj_obj.cast<SJJunctions&>();
        const size_t offset = offsetof(majiq::SJJunction, contig.contig_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::SJJunction>(
            sj.data(), offset, sj_obj);
        },
        "array[int] of indexes indicating contig junction belongs to")
    .def_property_readonly("start",
        [](py::object& sj_obj) -> py::array_t<position_t> {
        SJJunctions& sj = sj_obj.cast<SJJunctions&>();
        const size_t offset = offsetof(majiq::SJJunction, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, majiq::SJJunction>(
            sj.data(), offset, sj_obj);
        },
        "array[int] for junction start")
    .def_property_readonly("end",
        [](py::object& sj_obj) -> py::array_t<position_t> {
        SJJunctions& sj = sj_obj.cast<SJJunctions&>();
        const size_t offset = offsetof(majiq::SJJunction, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, majiq::SJJunction>(
            sj.data(), offset, sj_obj);
        },
        "array[int] for junction end")
    .def_property_readonly("strand",
        [](py::object& sj_obj) -> py::array_t<std::array<char, 1>> {
        SJJunctions& sj = sj_obj.cast<SJJunctions&>();
        const size_t offset = offsetof(majiq::SJJunction, strand);
        return ArrayFromVectorAndOffset<std::array<char, 1>, majiq::SJJunction>(
            sj.data(), offset, sj_obj);
        },
        "array[char] for junction strand")
    .def_property_readonly("numreads",
        [](py::object& sj_obj) -> py::array_t<majiq::junction_ct_t> {
        SJJunctions& sj = sj_obj.cast<SJJunctions&>();
        const size_t offset = offsetof(majiq::SJJunction, data.numreads);
        return ArrayFromVectorAndOffset<majiq::junction_ct_t, majiq::SJJunction>(
            sj.data(), offset, sj_obj);
        },
        "array[int] for total number of reads")
    .def_property_readonly("numpos",
        [](py::object& sj_obj) -> py::array_t<majiq::junction_pos_t> {
        SJJunctions& sj = sj_obj.cast<SJJunctions&>();
        const size_t offset = offsetof(majiq::SJJunction, data.numpos);
        return ArrayFromVectorAndOffset<majiq::junction_pos_t, majiq::SJJunction>(
            sj.data(), offset, sj_obj);
        },
        "array[int] for total number of nonzero positions")
    .def("__len__", &SJJunctions::size, "Number of junctions")
    .def_static("from_netcdf", &JunctionsFromNetcdf,
        "Load junctions from netcdf",
        py::arg("netcdf_path"))
    .def("df",
        [](py::object& sj) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(sj, "jidx",
            {"contig_idx", "start", "end", "strand"},
            {"numreads", "numpos"});
        },
        "View on junction information as xarray Dataset");

  pySJJunctionsPositions
    .def_property_readonly("_junctions", &SJJunctionsPositions::junctions,
        "Underlying junctions and summarized counts")
    .def_property_readonly("num_positions",
        &SJJunctionsPositions::num_positions,
        "Number of valid positions possible (function of max read length)")
    .def("__len__", &SJJunctionsPositions::size,
        "Number of junction positions")
    .def_property_readonly("position_reads",
        [](py::object& sj_obj) {
        SJJunctionsPositions& sj = sj_obj.cast<SJJunctionsPositions&>();
        const size_t offset = offsetof(majiq::PositionReads, reads);
        return ArrayFromVectorAndOffset<majiq::junction_ct_t, majiq::PositionReads>(
            sj.reads(), offset, sj_obj);
        },
        "Number of reads for a junction/position")
    .def_property_readonly("position",
        [](py::object& sj_obj) {
        SJJunctionsPositions& sj = sj_obj.cast<SJJunctionsPositions&>();
        const size_t offset = offsetof(majiq::PositionReads, pos);
        return ArrayFromVectorAndOffset<majiq::junction_pos_t, majiq::PositionReads>(
            sj.reads(), offset, sj_obj);
        },
        "Position index for junction/position")
    .def_property_readonly("_offsets",
        [](py::object& sj_obj) {
        SJJunctionsPositions& sj = sj_obj.cast<SJJunctionsPositions&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            sj.offsets(), 0, sj_obj);
        },
        "Raw offsets for jpidx_start and jpidx_end")
    .def_property_readonly("jpidx_start",
        [](py::object& sj_obj) {
        SJJunctionsPositions& sj = sj_obj.cast<SJJunctionsPositions&>();
        return ArrayFromOffsetsVector<size_t>(sj.offsets(), true, sj_obj);
        },
        "First index into junctions/positions for each junction")
    .def_property_readonly("jpidx_end",
        [](py::object& sj_obj) {
        SJJunctionsPositions& sj = sj_obj.cast<SJJunctionsPositions&>();
        return ArrayFromOffsetsVector<size_t>(sj.offsets(), false, sj_obj);
        },
        "One after last index into junctions/positions for each junction")
    .def("df",
        [](py::object& sj) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(sj, "jpidx", {},
            {"position_reads", "position"});
        },
        "View on junction-position information as xarray Dataset")
    .def_property_readonly("junctions",
        [](py::object& sj) -> py::object {
        using namespace py::literals;
        auto base = sj.attr("_junctions").attr("df")();
        auto get_xr = [&sj](py::str key) {
          return py::module_::import("xarray").attr("DataArray")(
              sj.attr(key), "dims"_a = "jidx");
        };
        return base.attr("assign_coords")(
            "jpidx_start"_a = get_xr("jpidx_start"),
            "jpidx_end"_a = get_xr("jpidx_end"));
        },
        "View on junction information as xarray Dataset")
    .def("to_netcdf",
        [](py::object& sj, py::str output_path) {
        // don't write to existing file
        py::object Path = py::module_::import("pathlib").attr("Path");
        if (Path(output_path).attr("exists")().cast<bool>()) {
          std::ostringstream oss;
          oss << "Cannot save result to already existing file "
              << output_path.cast<std::string>();
          throw std::invalid_argument(oss.str());
        }
        using namespace py::literals;
        // save contigs
        sj.attr("_junctions").attr("_contigs").attr("to_netcdf")(
            output_path, "mode"_a = "w");
        // save junctions
        sj.attr("_junctions").attr("df")().attr("to_netcdf")(
            output_path,
            "group"_a = py::str(SJ_JUNCTIONS_NC_GROUP),
            "mode"_a = "a");
        // save junction positions
        auto xr_jp
          = sj.attr("df")()
          .attr("assign_coords")(
              "_offsets"_a = py::module_::import("xarray").attr("DataArray")(
                sj.attr("_offsets"), "dims"_a = "offset_idx"))
          .attr("assign_attrs")("num_positions"_a = sj.attr("num_positions"));
        xr_jp.attr("to_netcdf")(
            output_path,
            "group"_a = py::str(SJ_JUNCTIONS_RAW_NC_GROUP),
            "mode"_a = "a");
        return;
        },
        "Serialize junction counts to specified file",
        py::arg("output_path"))
    .def_static("from_netcdf", &JunctionsRawFromNetcdf,
        "Load junctions and per-position counts from netcdf",
        py::arg("netcdf_path"))
    .def_static("from_bam",
        [](const char* infile, ExperimentStrandness strand, int nthreads) {
        return majiq::JunctionsFromBam<BAM_MIN_OVERHANG>(
            infile, nthreads, strand);
        },
        R"pbdoc(
        Load junctions and per-position counts for an aligned BAM file

        Parameters
        ----------
        bam_path: str
            Path for input BAM fille
        experiment_strandness: ExperimentStrandness
            Strandness of RNA-seq library
        nthreads: int
            Number of threads to use when reading in BAM file
        )pbdoc",
        py::arg("bam_path"),
        py::arg("experiment_strandness") = ExperimentStrandness::NONE,
        py::arg("nthreads") = 1);
}
