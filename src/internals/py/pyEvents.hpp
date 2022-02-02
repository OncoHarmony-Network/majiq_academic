/**
 * pyEvents.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYEVENTS_HPP
#define MAJIQ_PYBIND_PYEVENTS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <utility>
#include <vector>

#include "../GeneIntrons.hpp"
#include "../GeneJunctions.hpp"
#include "../Events.hpp"

#include "constants.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pyEvents_t = pyClassShared_t<Events>;

inline void init_Events(pyEvents_t& pyEvents) {
  pyEvents
    .def(pybind11::init([](
            const std::shared_ptr<majiq::GeneIntrons>& introns,
            const std::shared_ptr<majiq::GeneJunctions>& junctions,
            // make events vector (ref_exon_idx, event_type)
            pybind11::array_t<size_t> _ref_exon_idx,
            pybind11::array_t<std::array<char, 1>> _event_type,
            // connection offsets (1 longer than events)
            pybind11::array_t<size_t> _offsets,
            // connections (is_intron, connection_idx)
            pybind11::array_t<bool> _is_intron,
            pybind11::array_t<size_t> _connection_idx) {
          if (_ref_exon_idx.ndim() != 1
              || _event_type.ndim() != 1
              || _offsets.ndim() != 1
              || _is_intron.ndim() != 1
              || _connection_idx.ndim() != 1) {
            throw std::invalid_argument(
                "Events::init input arrays must be 1D");
          }
          if (_ref_exon_idx.shape(0) != _event_type.shape(0)) {
            throw std::runtime_error(
                "ref_exon_idx and event_type must have same length");
          }
          std::vector<Event> event_vec(_ref_exon_idx.shape(0));
          {
            auto ref_exon_idx = _ref_exon_idx.unchecked<1>();
            auto event_type = _event_type.unchecked<1>();
            for (size_t i = 0; i < event_vec.size(); ++i) {
              event_vec[i] = Event{
                ref_exon_idx(i), static_cast<EventType>(event_type(i)[0])};
            }
          }
          std::vector<size_t> offsets_vec(_offsets.shape(0));
          {
            auto offsets = _offsets.unchecked<1>();
            for (size_t i = 0; i < offsets_vec.size(); ++i) {
              offsets_vec[i] = offsets(i);
            }
          }
          if (_is_intron.shape(0) != _connection_idx.shape(0)) {
            throw std::runtime_error(
                "is_intron and connection_idx must have same length");
          }
          std::vector<ConnectionIndex> connections_vec(_is_intron.shape(0));
          {
            auto is_intron = _is_intron.unchecked<1>();
            auto connection_idx = _connection_idx.unchecked<1>();
            for (size_t i = 0; i < connections_vec.size(); ++i) {
              connections_vec[i] = ConnectionIndex{
                is_intron(i), connection_idx(i)};
            }
          }
          pybind11::gil_scoped_release release;  // release GIL at this stage
          return Events{introns, junctions, std::move(event_vec),
            std::move(offsets_vec), std::move(connections_vec)};
          }),
        "Initialize events object from numpy arrays",
        pybind11::arg("introns"),
        pybind11::arg("junctions"),
        pybind11::arg("ref_exon_idx"),
        pybind11::arg("event_type"),
        pybind11::arg("offsets"),
        pybind11::arg("is_intron"),
        pybind11::arg("connection_idx"))
    .def_property_readonly("introns", &Events::introns,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "underlying introns")
    .def_property_readonly("junctions", &Events::junctions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "underlying junctions")
    .def_property_readonly("exons", &Events::exons,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "underlying exons")
    .def("index",
        [](const Events& self,
          pybind11::array_t<size_t> ref_exon_idx,
          pybind11::array_t<bool> is_source) {
        auto f = [&self](size_t idx, bool is_src) -> std::ptrdiff_t {
          if (idx >= self.exons()->size()) {
            throw std::invalid_argument("ref_exon_idx has values out of range");
          }
          Event query{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT};
          auto lb = self.events_lower_bound(query);
          return (lb == self.events_end() || *lb != query)
              ? -1 : lb - self.events_begin();
        };
        return pybind11::vectorize(f)(ref_exon_idx, is_source);
        },
        "Get indexes for specified events (-1 if not found)",
        pybind11::arg("ref_exon_idx"),
        pybind11::arg("is_source"))
    .def_property_readonly("ref_exon_idx",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(Event, ref_exon_idx_);
        return ArrayFromVectorAndOffset<size_t, Event>(
            self.events(), offset, self_obj); },
        "Event reference exon")
    .def_property_readonly("event_type",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(Event, type_);
        return ArrayFromVectorAndOffset<std::array<char, 1>, Event>(
            self.events(), offset, self_obj); },
        "Event type")
    .def("event_id",
        [](const Events& self, pybind11::array_t<size_t> _event_idx) {
        if (_event_idx.ndim() != 1) {
          throw std::invalid_argument("event_idx must be 1D");
        }
        pybind11::list result{_event_idx.shape(0)};
        auto event_idx = _event_idx.unchecked<1>();
        for (size_t i = 0; i < result.size(); ++i) {
          if (event_idx(i) >= self.num_events()) {
            throw std::invalid_argument("event_idx has values out of range");
          }
          result[i] = self.event_id(event_idx(i));
        }
        return result;
        },
        "List of event_id for specified events",
        pybind11::arg("event_idx"))
    .def_property_readonly("_offsets",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.connection_offsets(), 0, self_obj);
        },
        "Raw offsets for events into connections")
    .def_property_readonly("connection_idx_start",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.connection_offsets(), true, self_obj); },
        "First index into event connections for each event")
    .def_property_readonly("connection_idx_end",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.connection_offsets(), false, self_obj); },
        "One after last index into event connections for each event")
    .def_property_readonly("_gene_offsets",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.gene_idx_offsets(), 0, self_obj);
        },
        "Raw offsets for genes into events")
    .def_property_readonly("event_idx_start",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.gene_idx_offsets(), true, self_obj); },
        "First index into gene events for each gene")
    .def_property_readonly("event_idx_end",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.gene_idx_offsets(), false, self_obj); },
        "One after last index into gene events for each gene")
    .def_property_readonly("is_intron",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(ConnectionIndex, is_intron_);
        return ArrayFromVectorAndOffset<bool, ConnectionIndex>(
            self.connections(), offset, self_obj); },
        "Event connection is intron (false --> is junction)")
    .def_property_readonly("idx",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(ConnectionIndex, idx_);
        return ArrayFromVectorAndOffset<size_t, ConnectionIndex>(
            self.connections(), offset, self_obj); },
        "Event connection index into corresponding Introns or Junctions")
    .def_property_readonly("connection_event_idx",
        [](pybind11::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.connection_event_idx(), 0, self_obj); },
        "Event connection index back into events")
    .def_property_readonly("num_events", &Events::num_events,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def_property_readonly("num_connections", &Events::num_connections,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def_property_readonly("num_junctions", &Events::num_junctions,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def_property_readonly("num_introns", &Events::num_introns,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("connection_gene_idx",
        [](const Events& self, pybind11::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        if (idx >= self.num_connections()) {
          throw std::invalid_argument("connection_idx has values out of range");
        }
        return self.connection_gene(idx).idx_; };
        return pybind11::vectorize(f)(connection_idx);
        },
        "gene_idx for specified connection indexes",
        pybind11::arg("connection_idx"))
    .def("connection_start",
        [](const Events& self, pybind11::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        if (idx >= self.num_connections()) {
          throw std::invalid_argument("connection_idx has values out of range");
        }
        return self.connection_start(idx); };
        return pybind11::vectorize(f)(connection_idx);
        },
        "start for specified connection indexes",
        pybind11::arg("connection_idx"))
    .def("connection_end",
        [](const Events& self, pybind11::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        if (idx >= self.num_connections()) {
          throw std::invalid_argument("connection_idx has values out of range");
        }
        return self.connection_end(idx); };
        return pybind11::vectorize(f)(connection_idx);
        },
        "end for specified connection indexes",
        pybind11::arg("connection_idx"))
    .def("connection_denovo",
        [](const Events& self, pybind11::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        if (idx >= self.num_connections()) {
          throw std::invalid_argument("connection_idx has values out of range");
        }
        return self.connection_denovo(idx); };
        return pybind11::vectorize(f)(connection_idx);
        },
        "denovo status for specified connection indexes",
        pybind11::arg("connection_idx"))
    .def("connection_other_exon_idx",
        [](const Events& self, pybind11::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        if (idx >= self.num_connections()) {
          throw std::invalid_argument("connection_idx has values out of range");
        }
        return self.connection_other_exon_idx(idx); };
        return pybind11::vectorize(f)(connection_idx);
        },
        "index for other exon for specified connection indexes",
        pybind11::arg("connection_idx"))
    .def("__len__",
        &Events::size,
        pybind11::call_guard<pybind11::gil_scoped_release>());
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYEVENTS_HPP
