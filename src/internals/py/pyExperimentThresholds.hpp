/**
 * pyExperimentThresholds.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYEXPERIMENTTHRESHOLDS_HPP
#define MAJIQ_PYBIND_PYEXPERIMENTTHRESHOLDS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <utility>
#include <vector>

#include "../ExperimentThresholds.hpp"

#include "constants.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pyExperimentThresholds_t = pyClassShared_t<ExperimentThresholds>;
using pyIntronThresholdsGenerator_t
  = pyClassShared_t<majiq::IntronThresholdsGenerator>;

inline void init_ExperimentThresholds(
    pyExperimentThresholds_t& pyExperimentThresholds) {
  pyExperimentThresholds
    .def(
        pybind11::init<
          junction_ct_t, junction_ct_t, junction_pos_t, real_t, real_t, real_t
        >(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Thresholds on intron/junction coverage for inclusion in SpliceGraph",
        pybind11::arg("minreads") = DEFAULT_BUILD_MINREADS,
        pybind11::arg("mindenovo") = DEFAULT_BUILD_MINDENOVO,
        pybind11::arg("minpos") = DEFAULT_BUILD_MINPOS,
        pybind11::arg("max_pctbins") = DEFAULT_BUILD_MAX_PCTBINS,
        pybind11::arg("junction_acceptance_probability")
          = DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
        pybind11::arg("intron_acceptance_probability")
          = DEFAULT_BUILD_MATCH_INTRON_PROBABILITY)
    .def_property_readonly("minreads",
        [](const ExperimentThresholds& x) { return x.minreads_; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Minimum number of reads for an annotated junction to pass")
    .def_property_readonly("mindenovo",
        [](const ExperimentThresholds& x) { return x.mindenovo_; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Minimum number of reads for a denovo junction to pass")
    .def_property_readonly("minpos",
        [](const ExperimentThresholds& x) { return x.minpos_; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Minimum number of nonzero positions for a junction to pass")
    .def_property_readonly("max_pctbins",
        [](const ExperimentThresholds& x) { return x.max_pctbins_; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Maximum percentage of bins to require coverage in for intron to pass")
    .def_property_readonly("junction_acceptance_probability",
        [](const ExperimentThresholds& x) {
        return x.junction_acceptance_probability_; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Set intron thresholds to match junction readrate with probability")
    .def_property_readonly("intron_acceptance_probability",
        [](const ExperimentThresholds& x) {
        return x.intron_acceptance_probability_; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Intron thresholds pass per-position readrate with probability")
    .def("__repr__", [](const ExperimentThresholds& x) -> std::string {
        std::ostringstream oss;
        oss << "ExperimentThresholds(minreads=" << x.minreads_
            << ", mindenovo=" << x.mindenovo_
            << ", minpos=" << x.minpos_
            << ", max_pctbins=" << x.max_pctbins_
            << ", junction_acceptance_probability="
            << x.junction_acceptance_probability_
            << ", intron_acceptance_probability="
            << x.intron_acceptance_probability_
            << ")";
        return oss.str();
        },
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("intron_thresholds_generator",
        &ExperimentThresholds::intron_thresholds_generator,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        Create IntronThresholdsGenerator for intron thresholds by length

        Parameters
        ----------
        total_bins: int
            Maximum number of bins/positions for each junction/intron for
            quantification
        )pbdoc",
        pybind11::arg("total_bins"));
}

inline void init_IntronThresholdsGenerator(
    pyIntronThresholdsGenerator_t& pyIntronThresholdsGenerator) {
  PYBIND11_NUMPY_DTYPE(IntronThresholds, minreads_, minbins_, mincov_);
  pyIntronThresholdsGenerator
    .def("__call__",
        [](const IntronThresholdsGenerator& gen,
          const pybind11::array_t<junction_pos_t>& intron_lengths) {
        // function per-element of intron_lengths
        auto f = [&gen](junction_pos_t x) { return gen(x); };
        return pybind11::vectorize(f)(intron_lengths);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get intron thresholds for specified intron lengths",
        pybind11::arg("intron_lengths"));
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYEXPERIMENTTHRESHOLDS_HPP
