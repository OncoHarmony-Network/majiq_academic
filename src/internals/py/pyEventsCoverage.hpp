/**
 * pyEventsCoverage.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYEVENTSCOVERAGE_HPP
#define MAJIQ_PYBIND_PYEVENTSCOVERAGE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <utility>
#include <vector>

#include "../EventsCoverage.hpp"

#include "globalRNG.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"


namespace majiq {
namespace bindings {

using pyEventsCoverage_t = pyClassShared_t<EventsCoverage>;

inline void init_EventsCoverage(pyEventsCoverage_t& pyEventsCoverage) {
  pyEventsCoverage
    .def(pybind11::init([](
            const std::shared_ptr<Events>& events,
            pybind11::array_t<majiq::real_t> _numreads,
            pybind11::array_t<majiq::real_t> _numbins,
            pybind11::array_t<majiq::real_t> _bootstraps) {
          if (_numreads.ndim() != 1 || _numbins.ndim() != 1) {
            throw std::invalid_argument("numreads and numbins most both be 1D");
          }
          if (_bootstraps.ndim() != 2) {
            throw std::invalid_argument("bootstraps must be 2D");
          }
          if (_numreads.shape(0) != _numbins.shape(0)
              || _numbins.shape(0) != _bootstraps.shape(0)) {
            throw std::invalid_argument(
                "EventsCoverage arrays do not agree on first dimension");
          }
          std::vector<CoverageSummary> summaries_vec(_numreads.shape(0));
          {
            auto numreads = _numreads.unchecked<1>();
            auto numbins = _numbins.unchecked<1>();
            for (size_t i = 0; i < summaries_vec.size(); ++i) {
              summaries_vec[i] = CoverageSummary{numreads(i), numbins(i)};
            }
          }
          CoverageBootstraps coverage{
              static_cast<size_t>(_bootstraps.shape(0)),
              static_cast<size_t>(_bootstraps.shape(1))};
          {
            auto bootstraps = _bootstraps.unchecked<2>();
            for (size_t i = 0; i < coverage.num_connections(); ++i) {
              for (size_t j = 0; j < coverage.num_bootstraps(); ++j) {
                coverage(i, j) = bootstraps(i, j);
              }
            }
          }
          pybind11::gil_scoped_release release;  // release GIL at this stage
          return EventsCoverage{events,
              std::move(summaries_vec), std::move(coverage)};
          }),
        "construct events coverage from numpy arrays",
        pybind11::arg("events"),
        pybind11::arg("numreads"),
        pybind11::arg("numbins"),
        pybind11::arg("bootstraps"))
    .def_static("from_sj",
        [](
          const std::shared_ptr<Events>& events,
          const majiq::SJJunctionsBins& sj_junctions,
          const majiq::SJIntronsBins& sj_introns,
          size_t num_bootstraps,
          majiq::real_t pvalue_threshold) {
        // acquire ownership of random number generator
        // NOTE: need to keep pointer in scope to maintain ownership
        // (i.e. don't replace with *global_rng_pool.acquire())
        auto gen_ptr = global_rng_pool.acquire();
        auto& gen = *gen_ptr;
        return EventsCoverage::FromSJ(events, sj_junctions, sj_introns,
            num_bootstraps, gen, pvalue_threshold);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Obtain coverage for events from SJ junctions and introns",
        pybind11::arg("events"),
        pybind11::arg("sj_junctions"),
        pybind11::arg("sj_introns"),
        pybind11::arg("num_bootstraps") = DEFAULT_BUILD_NUM_BOOTSTRAPS,
        pybind11::arg("pvalue_threshold") = DEFAULT_BUILD_STACK_PVALUE)
    .def_property_readonly("numreads",
        [](pybind11::object& self_obj) -> pybind11::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const size_t offset = offsetof(CoverageSummary, numreads);
        return ArrayFromVectorAndOffset<majiq::real_t, CoverageSummary>(
            self.summaries(), offset, self_obj);
        },
        "Readrate of each event connection after stacks removed")
    .def_property_readonly("numbins",
        [](pybind11::object& self_obj) -> pybind11::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const size_t offset = offsetof(CoverageSummary, numbins);
        return ArrayFromVectorAndOffset<majiq::real_t, CoverageSummary>(
            self.summaries(), offset, self_obj);
        },
        "Number of nonzero bins for each connection after stacks removed")
    .def_property_readonly("bootstraps",
        [](pybind11::object& self_obj) -> pybind11::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const CoverageBootstraps& x = self.bootstraps();
        pybind11::array_t<majiq::real_t> result = pybind11::array_t(
            // shape
            {x.num_connections(), x.num_bootstraps()},
            // strides
            {sizeof(majiq::real_t) * x.num_bootstraps(), sizeof(majiq::real_t)},
            // pointer to first element
            x.data().data(),
            // handle for object
            self_obj);
        return result;
        },
        "Bootstrapped read coverage or each connection after stacks removed")
    .def_property_readonly("_events", &EventsCoverage::events,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Events for which the coverage information is defined")
    .def("__len__", &EventsCoverage::num_connections,
        pybind11::call_guard<pybind11::gil_scoped_release>());
  return;
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYEVENTSCOVERAGE_HPP
