/**
 * pySJJunctions.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSJJUNCTIONS_HPP
#define MAJIQ_PYBIND_PYSJJUNCTIONS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../SJJunctions.hpp"

#include "constants.hpp"
#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pySJJunctions_t = pyClassShared_t<SJJunctions>;

inline void init_SJJunctions(pySJJunctions_t& pySJJunctions) {
  define_coordinates_properties(pySJJunctions);
  pySJJunctions
    .def(pybind11::init([](
            std::shared_ptr<Contigs> contigs,
            pybind11::array_t<size_t> _contig_idx,
            pybind11::array_t<position_t> _start,
            pybind11::array_t<position_t> _end,
            pybind11::array_t<std::array<char, 1>> _strand) {
          if (_contig_idx.ndim() != 1
              || _start.ndim() != 1
              || _end.ndim() != 1
              || _strand.ndim() != 1) {
            throw std::invalid_argument(
                "SJJunctions::init input arrays must be 1D");
          }
          if (_contig_idx.shape(0) != _start.shape(0)
              || _contig_idx.shape(0) != _end.shape(0)
              || _contig_idx.shape(0) != _strand.shape(0)) {
            throw std::invalid_argument(
                "SJJunctions::init input arrays have conflicting sizes");
          }
          auto contig_idx = _contig_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto strand = _strand.unchecked<1>();
          // fill in junctions
          std::vector<SJJunction> sj_vec(start.shape(0));
          for (size_t i = 0; i < sj_vec.size(); ++i) {
            auto strand_i = static_cast<GeneStrandness>(strand(i)[0]);
            if (strand_i != GeneStrandness::AMBIGUOUS
                && strand_i != GeneStrandness::FORWARD
                && strand_i != GeneStrandness::REVERSE) {
              throw std::invalid_argument(
                  "SJJunctions::init input arrays have invalid strand");
            }
            sj_vec[i] = SJJunction{
              KnownContig{contig_idx(i), contigs},
              OpenInterval{start(i), end(i)},
              strand_i};
          }
          return std::make_shared<SJJunctions>(
              contigs, std::move(sj_vec));
        }),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Create SJJunctions object from contigs and arrays",
        pybind11::arg("contigs"),
        pybind11::arg("contig_idx"),
        pybind11::arg("start"), pybind11::arg("end"),
        pybind11::arg("strand"))
    .def("to_unstranded", &SJJunctions::ToUnstranded,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Create unstranded SJJunctions from self")
    .def("flip_strand", &SJJunctions::FlipStrand,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Create SJJunctions with strands flipped in sorted order")
    .def_property_readonly("_contigs", &SJJunctions::parents,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying contigs corresponding to contig_idx");
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSJJUNCTIONS_HPP
