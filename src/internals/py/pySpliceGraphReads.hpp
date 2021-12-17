/**
 * pySpliceGraphReads.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSPLICEGRAPHREADS_HPP
#define MAJIQ_PYBIND_PYSPLICEGRAPHREADS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <utility>
#include <vector>

#include "../SpliceGraphReads.hpp"

#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pySpliceGraphReads_t = pyClassShared_t<majiq::SpliceGraphReads>;

inline void init_SpliceGraphReads(pySpliceGraphReads_t& pySpliceGraphReads) {
  pySpliceGraphReads
    .def_property_readonly("_introns",
        &SpliceGraphReads::introns,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying introns")
    .def_property_readonly("_junctions",
        &SpliceGraphReads::junctions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying junctions")
    .def_property_readonly("introns_reads",
        [](pybind11::object& self_obj) {
        SpliceGraphReads& self = self_obj.cast<SpliceGraphReads&>();
        return ArrayFromVectorAndOffset<majiq::real_t, majiq::real_t>(
            self.introns_reads(), 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Raw readrates for each intron")
    .def_property_readonly("junctions_reads",
        [](pybind11::object& self_obj) {
        SpliceGraphReads& self = self_obj.cast<SpliceGraphReads&>();
        return ArrayFromVectorAndOffset<majiq::real_t, majiq::real_t>(
            self.junctions_reads(), 0, self_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Raw readrates for each junction")
    .def(pybind11::init([](
            const std::shared_ptr<GeneIntrons>& introns,
            const std::shared_ptr<GeneJunctions>& junctions,
            pybind11::array_t<majiq::real_t> _introns_reads,
            pybind11::array_t<majiq::real_t> _junctions_reads) {
          if (_introns_reads.ndim() != 1 || _junctions_reads.ndim() != 1) {
            throw std::invalid_argument(
                "introns/junctions reads must both be 1D");
          }
          std::vector<majiq::real_t> ireads_vec(_introns_reads.shape(0));
          {
            auto introns_reads = _introns_reads.unchecked<1>();
            for (size_t i = 0; i < ireads_vec.size(); ++i) {
            ireads_vec[i] = introns_reads(i);
            }
          }
          std::vector<majiq::real_t> jreads_vec(_junctions_reads.shape(0));
          {
            auto junctions_reads = _junctions_reads.unchecked<1>();
            for (size_t j = 0; j < jreads_vec.size(); ++j) {
            jreads_vec[j] = junctions_reads(j);
            }
          }
          return SpliceGraphReads{introns, junctions,
              std::move(ireads_vec), std::move(jreads_vec)};
          }),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Initialize SpliceGraphReads from numpy arrays",
        pybind11::arg("introns"),
        pybind11::arg("junctions"),
        pybind11::arg("introns_reads"),
        pybind11::arg("junctions_reads"))
    .def_static("from_sj", &SpliceGraphReads::FromSJ,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Obtain raw readrates for introns/junctions from experiment SJ",
        pybind11::arg("introns"),
        pybind11::arg("junctions"),
        pybind11::arg("sj_introns"),
        pybind11::arg("sj_junctions"));
}


}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSPLICEGRAPHREADS_HPP
