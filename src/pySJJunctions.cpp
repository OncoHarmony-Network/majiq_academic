/**
 * pySJJunctions.cpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <array>
#include <memory>
#include <cstddef>

#include "pybind_utils.hpp"
#include "internals/SJJunctions.hpp"
#include "internals/MajiqTypes.hpp"


namespace py = pybind11;


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
    .def("jidx_to_jpidx",
        [](py::object& sj) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(
            sj, "jidx", {"jpidx_start", "jpidx_end"});
        },
        "Map from jidx to starts/ends for jpidx")
    .def_property_readonly("junctions",
        [](py::object& sj) -> py::object {
        auto base = sj.attr("_junctions").attr("df")();
        auto jpidx_map = sj.attr("jidx_to_jpidx")();
        return base.attr("combine_first")(jpidx_map);
        },
        "View on junction information as xarray Dataset")
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
