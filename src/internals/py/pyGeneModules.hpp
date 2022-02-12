/**
 * pyGeneModules.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYGENEMODULES_HPP
#define MAJIQ_PYBIND_PYGENEMODULES_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>

#include "../GeneModules.hpp"

#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyGeneModules_t = pyClassShared_t<GeneModules>;

inline void init_GeneModules(pyGeneModules_t& pyGeneModules) {
  define_coordinates_properties(pyGeneModules);
  pyGeneModules
    .def("event_is_masked",
        [](const GeneModules& self,
        pybind11::array_t<size_t> exon_idx,
        pybind11::array_t<bool> is_source) -> pybind11::array_t<size_t> {
        auto f = [&self](size_t idx, bool is_src) -> size_t {
          if (idx >= self.exon_connections()->num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.event_is_masked(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        "Return array[bool] indicating if selected events were masked",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("event_module_idx",
        [](const GeneModules& self,
        pybind11::array_t<size_t> exon_idx,
        pybind11::array_t<bool> is_source) -> pybind11::array_t<size_t> {
        auto f = [&self](size_t idx, bool is_src) -> size_t {
          if (idx >= self.exon_connections()->num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.event_module_index(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        "Return array[int] indicating module to which input events belong",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def_property_readonly("introns_module_idx",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.introns_module_idx(), 0, self_obj);
        },
        "array[int] indicating the module to which each intron belongs")
    .def_property_readonly("junctions_module_idx",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.junctions_module_idx(), 0, self_obj);
        },
        "array[int] indicating the module to which each junction belongs")
    .def_property_readonly("start_exon_idx",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        const size_t offset = offsetof(GeneModule, data.start_exon_idx);
        return ArrayFromVectorAndOffset<size_t, GeneModule>(
            self.data(), offset, self_obj);
        },
        "array[int] of start exon idx of each module")
    .def_property_readonly("end_exon_idx",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        const size_t offset = offsetof(GeneModule, data.end_exon_idx);
        return ArrayFromVectorAndOffset<size_t, GeneModule>(
            self.data(), offset, self_obj);
        },
        "array[int] of end exon idx of each module")
    .def_property_readonly("start_intron_idx",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        const size_t offset = offsetof(GeneModule, data.start_intron_idx);
        return ArrayFromVectorAndOffset<size_t, GeneModule>(
            self.data(), offset, self_obj);
        },
        "array[int] of start intron idx of each module")
    .def_property_readonly("num_introns",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        const size_t offset = offsetof(GeneModule, data.num_introns);
        return ArrayFromVectorAndOffset<size_t, GeneModule>(
            self.data(), offset, self_obj);
        },
        "array[int] of number of introns for each module")
    .def_property_readonly("start_junction_idx",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        const size_t offset = offsetof(GeneModule, data.start_junction_idx);
        return ArrayFromVectorAndOffset<size_t, GeneModule>(
            self.data(), offset, self_obj);
        },
        "array[int] of start junction idx of each module")
    .def_property_readonly("num_junctions",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GeneModules& self = self_obj.cast<GeneModules&>();
        const size_t offset = offsetof(GeneModule, data.num_junctions);
        return ArrayFromVectorAndOffset<size_t, GeneModule>(
            self.data(), offset, self_obj);
        },
        "array[int] of number of junctions for each module")
    .def_property_readonly("exon_connections",
        &GeneModules::exon_connections,
        "underlying exon connections")
    .def_property_readonly("mask",
        &GeneModules::mask,
        "underlying mask over splicegraph introns and junctions")
    .def(
        pybind11::init<
          std::shared_ptr<ExonConnections>,
          std::shared_ptr<SpliceGraphMask>
        >(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Initialize GeneModules from exon connections and splicegraph mask",
        pybind11::arg("exon_connections"),
        pybind11::arg("sg_mask"));
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYGENEMODULES_HPP
