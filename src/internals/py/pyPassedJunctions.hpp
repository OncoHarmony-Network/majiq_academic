/**
 * pyPassedJunctions.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYPASSEDJUNCTIONS_HPP
#define MAJIQ_PYBIND_PYPASSEDJUNCTIONS_HPP

#include <pybind11/pybind11.h>

#include <memory>
#include <utility>
#include <vector>

#include "../PassedJunctions.hpp"

#include "constants.hpp"
#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyGroupJunctionsGenerator_t = pyClassShared_t<GroupJunctionsGenerator>;
using pyPassedJunctionsGenerator_t = pyClassShared_t<PassedJunctionsGenerator>;

inline void init_GroupJunctionsGenerator(
    pyGroupJunctionsGenerator_t& pyGroupJunctionsGenerator) {
  pyGroupJunctionsGenerator
    .def(
      pybind11::init<const std::shared_ptr<majiq::GeneJunctions>&,
      const std::shared_ptr<majiq::Exons>&>(),
      pybind11::call_guard<pybind11::gil_scoped_release>(),
      "Accumulate SJJunctions for build group for input junctions/exons",
      pybind11::arg("junctions"),
      pybind11::arg("exons"))
    .def_property_readonly("_known", &GroupJunctionsGenerator::known,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get underlying known junctions")
    .def("add_experiment", &majiq::GroupJunctionsGenerator::AddExperiment,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Increment count of passed junctions from input experiment",
        pybind11::arg("sjp"),
        pybind11::arg("thresholds") = DEFAULT_THRESHOLDS,
        pybind11::arg("add_denovo") = DEFAULT_BUILD_DENOVO_JUNCTIONS)
    .def("pass_known_inplace",
        &majiq::GroupJunctionsGenerator::UpdateKnownInplace,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Update known junctions with their build status (ignore denovos)",
        pybind11::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def_property_readonly("num_experiments",
        &majiq::GroupJunctionsGenerator::num_experiments,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments that have been added to this group")
    .def_property_readonly("num_known",
        &majiq::GroupJunctionsGenerator::num_annotated,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of junctions known when constructed")
    .def_property_readonly("num_denovo",
        &majiq::GroupJunctionsGenerator::num_denovo,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of denovo junctions passing experiment-filters for 1+ inputs")
    .def("__len__", &majiq::GroupJunctionsGenerator::size,
        pybind11::call_guard<pybind11::gil_scoped_release>());
}

inline void init_PassedJunctionsGenerator(
    pyPassedJunctionsGenerator_t& pyPassedJunctionsGenerator) {
  pyPassedJunctionsGenerator
    .def(
      pybind11::init<const std::shared_ptr<majiq::GeneJunctions>&>(),
      pybind11::call_guard<pybind11::gil_scoped_release>(),
      R"pbdoc(
      Accumulator of GroupJunctionsGenerator for different build groups
      )pbdoc",
      pybind11::arg("junctions"))
    .def_property_readonly("_known", &PassedJunctionsGenerator::known,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get underlying known junctions")
    .def("add_group", &majiq::PassedJunctionsGenerator::AddGroup,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Combine passed junctions from build group of experiments",
        pybind11::arg("group"),
        pybind11::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def("add_junction",
        [](PassedJunctionsGenerator& self, size_t gene_idx,
          position_t start, position_t end) {
        self.AddJunction(gene_idx, start, end);
        return;
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Mark specific junction as passed",
        pybind11::arg("gene_idx"),
        pybind11::arg("start"),
        pybind11::arg("end"))
    .def("get_passed", &majiq::PassedJunctionsGenerator::PassedJunctions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get static, array-based representation of passed junctions",
        pybind11::arg("denovo_simplified"))
    .def_property_readonly("num_known",
        &majiq::PassedJunctionsGenerator::num_annotated,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of junctions known when constructed")
    .def_property_readonly("num_denovo",
        &majiq::PassedJunctionsGenerator::num_denovo,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of denovo junctions passing filters, previously not known")
    .def("__len__", &majiq::PassedJunctionsGenerator::size,
        pybind11::call_guard<pybind11::gil_scoped_release>());
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYPASSEDJUNCTIONS_HPP
