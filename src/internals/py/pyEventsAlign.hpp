/**
 * pyEventsAlign.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYEVENTSALIGN_HPP
#define MAJIQ_PYBIND_PYEVENTSALIGN_HPP

#include <pybind11/pybind11.h>

#include <memory>

#include "../EventsAlign.hpp"

#include "templateTypes.hpp"
#include "vectorArray.hpp"


namespace majiq {
namespace bindings {

using pyEventsAlign_t = pyClassShared_t<EventsAlign>;

inline void init_EventsAlign(pyEventsAlign_t& pyEventsAlign) {
  pyEventsAlign
    .def(pybind11::init<const Events&, const Events&>(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Obtain indexes of matching events in the two input Events containers",
        pybind11::arg("left_events"), pybind11::arg("right_events"))
    .def_static(
        "events_match",
        [](const Events& left_events, const Events& right_events,
          size_t left_idx, size_t right_idx) {
        if (left_idx >= left_events.num_events()) {
        throw std::runtime_error("left_idx is out of range for left_events");
        } else if (right_idx >= right_events.num_events()) {
        throw std::runtime_error("right_idx is out of range for right_events");
        }
        return EventsAlign::EventsMatch(
            left_events, right_events, left_idx, right_idx);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indicate if selected event indexes share the exact same connections",
        pybind11::arg("left_events"), pybind11::arg("right_events"),
        pybind11::arg("left_idx"), pybind11::arg("right_idx"))
    .def_property_readonly("left_event_idx",
        [](pybind11::object& self_obj) {
        EventsAlign& self = self_obj.cast<EventsAlign&>();
        const size_t offset = offsetof(EventsAlign::EventAligned, left_idx_);
        return ArrayFromVectorAndOffset<size_t, EventsAlign::EventAligned>(
            self.matched_, offset, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indexes for events in left_events used in constructor")
    .def_property_readonly("right_event_idx",
        [](pybind11::object& self_obj) {
        EventsAlign& self = self_obj.cast<EventsAlign&>();
        const size_t offset = offsetof(EventsAlign::EventAligned, right_idx_);
        return ArrayFromVectorAndOffset<size_t, EventsAlign::EventAligned>(
            self.matched_, offset, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Indexes for events in right_events used in constructor");
  return;
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYEVENTSALIGN_HPP
