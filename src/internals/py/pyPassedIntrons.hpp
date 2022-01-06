/**
 * pyPassedIntrons.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYPASSEDINTRONS_HPP
#define MAJIQ_PYBIND_PYPASSEDINTRONS_HPP

#include <pybind11/pybind11.h>

#include <memory>
#include <utility>
#include <vector>

#include "../PassedIntrons.hpp"

#include "constants.hpp"
#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyGroupIntronsGenerator_t = pyClassShared_t<GroupIntronsGenerator>;

inline void init_GroupIntronsGenerator(
    pyGroupIntronsGenerator_t& pyGroupIntronsGenerator) {
  pyGroupIntronsGenerator
    .def(
        pybind11::init<const std::shared_ptr<majiq::GeneIntrons>&>(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Accumulate SJIntronsBins for build group of introns",
        pybind11::arg("gene_introns"))
    .def_property_readonly("num_experiments",
        &GroupIntronsGenerator::num_experiments,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments in current group")
    .def_property_readonly("_introns", &GroupIntronsGenerator::introns,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def_property_readonly("num_passed",
        [](pybind11::object& self_obj) -> pybind11::array_t<size_t> {
        GroupIntronsGenerator& self = self_obj.cast<GroupIntronsGenerator&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.num_passed(), 0, self_obj);
        },
        "Number of experiments or which each intron has passed")
    .def("__len__", &GroupIntronsGenerator::size,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("add_experiment", &GroupIntronsGenerator::AddExperiment,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Add SJIntronsBins to build group",
        pybind11::arg("sj"),
        pybind11::arg("thresholds") = DEFAULT_THRESHOLDS)
    .def("update_introns", &GroupIntronsGenerator::UpdateInplace,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Pass introns that pass min-experiments in place, reset for next group",
        pybind11::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS);
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYPASSEDINTRONS_HPP
