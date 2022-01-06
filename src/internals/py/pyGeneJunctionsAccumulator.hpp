/**
 * pyGeneJunctionsAccumulator.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYGENEJUNCTIONSACCUMULATOR_HPP
#define MAJIQ_PYBIND_PYGENEJUNCTIONSACCUMULATOR_HPP

#include <pybind11/pybind11.h>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../GeneJunctionsAccumulator.hpp"

#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyGeneJunctionsAccumulator_t
  = pyClassShared_t<GeneJunctionsAccumulator>;

inline void init_GeneJunctionsAccumulator(
    pyGeneJunctionsAccumulator_t& pyGeneJunctionsAccumulator) {
  pyGeneJunctionsAccumulator
    .def(pybind11::init<const std::shared_ptr<Genes>&>(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Initialize accumulator of GeneJunctions with Genes",
        pybind11::arg("genes"))
    .def_property_readonly(
        "_genes", &GeneJunctionsAccumulator::genes,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Genes used by GeneJunctionsAccumulator")
    .def(
        "add",
        &GeneJunctionsAccumulator::Add,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Add/update accumulated junctions",
        pybind11::arg("junctions"),
        pybind11::arg("make_annotated") = false)
    .def(
        "accumulated",
        &GeneJunctionsAccumulator::Accumulated,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Return GeneJunctions combining accumulated junctions");
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYGENEJUNCTIONSACCUMULATOR_HPP
