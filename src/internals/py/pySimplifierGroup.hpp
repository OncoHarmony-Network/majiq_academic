/**
 * pySimplifierGroup.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSIMPLIFIERGROUP_HPP
#define MAJIQ_PYBIND_PYSIMPLIFIERGROUP_HPP

#include <pybind11/pybind11.h>

#include <memory>

#include "../SimplifierGroup.hpp"

#include "constants.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pySimplifierGroup_t = pyClassShared_t<SimplifierGroup>;

inline void init_SimplifierGroup(pySimplifierGroup_t& pySimplifierGroup) {
  pySimplifierGroup
    .def_property_readonly("_exon_connections",
        &SimplifierGroup::exon_connections,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying exon connections being unsimplified")
    .def_property_readonly("introns_passed_src",
        [](pybind11::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, src_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.introns_passed(), offset, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as source for introns")
    .def_property_readonly("introns_passed_dst",
        [](pybind11::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, dst_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.introns_passed(), offset, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as target for introns")
    .def_property_readonly("junctions_passed_src",
        [](pybind11::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, src_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.junctions_passed(), offset, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as source for junctions")
    .def_property_readonly("junctions_passed_dst",
        [](pybind11::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, dst_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.junctions_passed(), offset, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as target for junctions")
    .def_property_readonly("num_experiments",
        &SimplifierGroup::num_experiments,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of experiments in current simplifier group")
    .def("add_experiment", &SimplifierGroup::AddExperiment,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Increment evidence to unsimplify connections from SpliceGraphReads",
        pybind11::arg("sg_reads"),
        pybind11::arg("simplify_min_psi") = DEFAULT_BUILD_SIMPL_MINPSI,
        pybind11::arg("simplify_minreads_annotated_junctions")
          = DEFAULT_BUILD_SIMPL_MINREADS_ANNOTATED_JUNCTION,
        pybind11::arg("simplify_minreads_denovo_junctions")
          = DEFAULT_BUILD_SIMPL_MINREADS_DENOVO_JUNCTION,
        pybind11::arg("simplify_minreads_introns")
          = DEFAULT_BUILD_SIMPL_MINREADS_INTRON)
    .def("update_connections", &SimplifierGroup::UpdateInplace,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Unsimplify introns/junctions with evidence, reset for next group",
        pybind11::arg("simplifier_min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def(pybind11::init<const std::shared_ptr<ExonConnections>&>(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Start simplification group for the specified exon connections",
        pybind11::arg("exon_connections"));
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSIMPLIFIERGROUP_HPP
