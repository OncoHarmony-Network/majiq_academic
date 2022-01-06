/**
 * pySJIntrons.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSJINTRONS_HPP
#define MAJIQ_PYBIND_PYSJINTRONS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../SJIntrons.hpp"

#include "constants.hpp"
#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pySJIntrons_t = pyClassShared_t<SJIntrons>;

inline void init_SJIntrons(pySJIntrons_t& pySJIntrons) {
  define_coordinates_properties(pySJIntrons);
  pySJIntrons
    .def(pybind11::init([](
            std::shared_ptr<Contigs> contigs,
            pybind11::array_t<size_t> _contig_idx,
            pybind11::array_t<position_t> _start,
            pybind11::array_t<position_t> _end,
            pybind11::array_t<std::array<char, 1>> _strand,
            pybind11::array_t<bool> _annotated) {
          if (_contig_idx.ndim() != 1
              || _start.ndim() != 1
              || _end.ndim() != 1
              || _strand.ndim() != 1
              || _annotated.ndim() != 1) {
            throw std::invalid_argument(
                "SJIntrons::init input arrays must be 1D");
          }
          if (_contig_idx.shape(0) != _start.shape(0)
              || _contig_idx.shape(0) != _end.shape(0)
              || _contig_idx.shape(0) != _strand.shape(0)
              || _contig_idx.shape(0) != _annotated.shape(0)) {
            throw std::invalid_argument(
                "SJIntrons::init input arrays have conflicting sizes");
          }
          auto contig_idx = _contig_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto strand = _strand.unchecked<1>();
          auto annotated = _annotated.unchecked<1>();
          std::vector<SJIntron> result(start.shape(0));
          for (size_t i = 0; i < result.size(); ++i) {
            auto strand_i = static_cast<GeneStrandness>(strand(i)[0]);
            if (strand_i != GeneStrandness::AMBIGUOUS
                && strand_i != GeneStrandness::FORWARD
                && strand_i != GeneStrandness::REVERSE) {
              throw std::invalid_argument(
                  "SJJunctions::init input arrays have invalid strand");
            }
            result[i] = SJIntron{
              majiq::KnownContig{contig_idx(i), contigs},
              majiq::ClosedInterval{start(i), end(i)},
              strand_i,
              annotated(i)};
          }
          pybind11::gil_scoped_release release;  // release GIL at this stage
          return std::make_shared<SJIntrons>(contigs, std::move(result));
        }),
        "Create SJIntrons object from contigs and arrays",
        pybind11::arg("contigs"),
        pybind11::arg("contig_idx"),
        pybind11::arg("start"), pybind11::arg("end"),
        pybind11::arg("strand"),
        pybind11::arg("annotated"))
    .def_static("from_exons_and_introns", &SJIntrons::FromGeneExonsAndIntrons,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Construct sj introns for input exons/introns",
        pybind11::arg("exons"), pybind11::arg("introns"),
        pybind11::arg("stranded"))
    .def_property_readonly("annotated",
        [](pybind11::object& introns_obj) -> pybind11::array_t<bool> {
        SJIntrons& introns = introns_obj.cast<SJIntrons&>();
        const size_t offset = offsetof(majiq::SJIntron, data.annotated_);
        return ArrayFromVectorAndOffset<bool, majiq::SJIntron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if intron is annotated (exon in annotation)");
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSJINTRONS_HPP
