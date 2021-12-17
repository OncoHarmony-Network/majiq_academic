/**
 * pyExonConnections.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYEXONCONNECTIONS_HPP
#define MAJIQ_PYBIND_PYEXONCONNECTIONS_HPP

#include <pybind11/pybind11.h>

#include <memory>
#include <vector>

#include "../ExonConnections.hpp"

#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pyExonConnections_t = pyClassShared_t<majiq::ExonConnections>;

inline void init_ExonConnections(pyExonConnections_t& pyExonConnections) {
  pyExonConnections
    .def(
        pybind11::init<
          const std::shared_ptr<Exons>&,
          const std::shared_ptr<GeneIntrons>&,
          const std::shared_ptr<GeneJunctions>&
        >(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Track junctions/introns by exons/events assocciated with them",
        pybind11::arg("exons"),
        pybind11::arg("introns"),
        pybind11::arg("junctions"))
    .def_property_readonly("src_intron_idx",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.idx_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "intron_idx for exon_connections in src_exon sorted order")
    .def_property_readonly("src_intron_exon_offsets",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.exon_offsets_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "offsets into src_intron_idx for each exon")
    .def_property_readonly("dst_intron_idx",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.idx_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "intron_idx for exon_connections in dst_exon sorted order")
    .def_property_readonly("dst_intron_exon_offsets",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.exon_offsets_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "offsets into dst_intron_idx for each exon")
    .def_property_readonly("src_junction_idx",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.idx_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "junction_idx for exon_connections in src_exon sorted order")
    .def_property_readonly("src_junction_exon_offsets",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.exon_offsets_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "offsets into src_junction_idx for each exon")
    .def_property_readonly("dst_junction_idx",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.idx_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "junction_idx for exon_connections in dst_exon sorted order")
    .def_property_readonly("dst_junction_exon_offsets",
        [](pybind11::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(
                indexes.exon_offsets_, 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "offsets into dst_junction_idx for each exon")
    .def_property_readonly("_exons", &ExonConnections::exons,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "underlying exons")
    .def_property_readonly("_introns", &ExonConnections::introns,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "underlying introns")
    .def_property_readonly("_junctions", &ExonConnections::junctions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "underlying junctions")
    .def("strict_lsvs", &ExonConnections::StrictLSVs,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct strict LSV Events")
    .def("permissive_lsvs", &ExonConnections::PermissiveLSVs,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct permissive LSV Events")
    .def("source_lsvs", &ExonConnections::SourceLSVs,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct source LSV Events")
    .def("target_lsvs", &ExonConnections::TargetLSVs,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct target LSV Events")
    .def("constitutive", &ExonConnections::ConstitutiveEvents,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct Constitutive Events")
    .def("events_for",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) {
        if (exon_idx.ndim() != 1) {
          throw std::invalid_argument("exon_idx must be one-dimensional");
        } else if (is_source.ndim() != 1) {
          throw std::invalid_argument("is_source must be one-dimensional");
        } else if (is_source.shape(0) != exon_idx.shape(0)) {
          throw std::invalid_argument(
              "exon_idx and is_source must have same size");
        }
        auto _exon_idx = exon_idx.unchecked<1>();
        auto _is_source = is_source.unchecked<1>();
        std::vector<Event> events(_exon_idx.shape(0));
        for (pybind11::ssize_t i = 0; i < _exon_idx.shape(0); ++i) {
          if (_exon_idx(i) >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has out of range values");
          }
          events[i] = Event{
            _exon_idx(i),
            _is_source(i) ? EventType::SRC_EVENT : EventType::DST_EVENT};
        }
        return self.CreateEvents(std::move(events));
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct events for specified exons/directions",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("has_intron",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.has_intron(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if events have introns or not",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("event_size",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<size_t> {
        auto f = [&self](size_t idx, bool is_src) -> size_t {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.event_size(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate size of event",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("passed",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.passed(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event was passed",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("redundant",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.redundant(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event was redundant",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("is_target_LSV",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_target) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_target_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_target);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event is target LSV",
        pybind11::arg("exon_idx"), pybind11::arg("is_target"))
    .def("is_source_LSV",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_source_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event is source LSV",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("is_permissive_LSV",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_permissive_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event is permissive LSV",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("is_strict_LSV",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_strict_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event is strict LSV",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("is_constitutive",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> exon_idx,
          pybind11::array_t<bool> is_source) -> pybind11::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_constitutive(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return pybind11::vectorize(f)(exon_idx, is_source);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if event is constitutive",
        pybind11::arg("exon_idx"), pybind11::arg("is_source"))
    .def("event_id",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> _exon_idx,
          pybind11::array_t<std::array<char, 1>> _event_type) {
        if (_exon_idx.ndim() != 1 || _event_type.ndim() != 1) {
        throw std::invalid_argument("exon_idx and event_type must be 1D");
        } else if (_exon_idx.shape(0) != _event_type.shape(0)) {
        throw std::invalid_argument(
                "exon_idx and event_type must have same shape");
        }
        std::vector<std::string> result(_exon_idx.shape(0));
        {
          auto exon_idx = _exon_idx.unchecked<1>();
          auto event_type = _event_type.unchecked<1>();
          for (size_t i = 0; i < result.size(); ++i) {
          if (exon_idx(i) >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          result[i] = self.id(Event{
              exon_idx(i), static_cast<EventType>(event_type(i)[0])});
          }
        }
        return result;
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "List of event_id for specified events",
        pybind11::arg("exon_idx"), pybind11::arg("event_type"))
    .def("event_description",
        [](const ExonConnections& self,
          pybind11::array_t<size_t> _exon_idx,
          pybind11::array_t<std::array<char, 1>> _event_type) {
        if (_exon_idx.ndim() != 1 || _event_type.ndim() != 1) {
        throw std::invalid_argument("exon_idx and event_type must be 1D");
        } else if (_exon_idx.shape(0) != _event_type.shape(0)) {
        throw std::invalid_argument(
                "exon_idx and event_type must have same shape");
        }
        std::vector<std::string> result(_exon_idx.shape(0));
        {
          auto exon_idx = _exon_idx.unchecked<1>();
          auto event_type = _event_type.unchecked<1>();
          for (size_t i = 0; i < result.size(); ++i) {
          if (exon_idx(i) >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          result[i] = self.description(Event{
              exon_idx(i), static_cast<EventType>(event_type(i)[0])});
          }
        }
        return result;
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "List of description for specified events",
        pybind11::arg("exon_idx"), pybind11::arg("event_type"));
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYEXONCONNECTIONS_HPP
