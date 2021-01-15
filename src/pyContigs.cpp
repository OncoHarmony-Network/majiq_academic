/**
 * pyContigs.cpp
 *
 * Set up python bindings to contigs
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>

#include "pybind_utils.hpp"
#include "internals/Contigs.hpp"

namespace py = pybind11;
using majiq::Contigs;

void init_Contigs(py::class_<Contigs, std::shared_ptr<Contigs>>& pyContigs) {
  using majiq_pybind::XarrayDatasetFromObject;
  using majiq_pybind::OpenXarrayDataset;
  using majiq_pybind::CONTIGS_NC_GROUP;
  using majiq::seqid_t;
  using majiq::Contig;
  using namespace py::literals;
  pyContigs
    .def_static("from_netcdf",
        [](py::str x) { return majiq_pybind::ContigsFromNetcdf(x); },
        "Load contigs from netcdf file", py::arg("netcdf_path"))
    .def("to_netcdf",
        [](py::object& self, py::str out, py::str mode) {
        self.attr("df")().attr("to_netcdf")(
            out, "mode"_a = mode, "group"_a = py::str(CONTIGS_NC_GROUP));
        return;
        },
        "Save contigs to netcdf file",
        py::arg("netcdf_path"), py::arg("mode"))
    .def_property_readonly("seqid", &Contigs::seqids,
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("df",
        [](py::object& contigs) -> py::object {
        return XarrayDatasetFromObject(contigs, "contig_idx", {"seqid"});
        },
        "View on contig information as xarray Dataset")
    .def("__repr__", [](const Contigs& self) -> std::string {
        std::ostringstream oss;
        oss << self;
        return oss.str();
        })
    .def("__len__", &Contigs::size)
    .def("__contains__",
        [](const Contigs& s, seqid_t x) -> bool { return s.contains(x); });
}
