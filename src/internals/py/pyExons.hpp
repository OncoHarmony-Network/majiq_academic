/**
 * pyExons.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYEXONS_HPP
#define MAJIQ_PYBIND_PYEXONS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../GeneIntrons.hpp"
#include "../Exons.hpp"

#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyExons_t = pyClassShared_t<majiq::Exons>;

inline void init_Exons(pyExons_t& pyExons) {
  define_coordinates_properties(pyExons);
  pyExons
    .def(pybind11::init([](
            std::shared_ptr<majiq::Genes> genes,
            pybind11::array_t<size_t> _gene_idx,
            pybind11::array_t<position_t> _start,
            pybind11::array_t<position_t> _end,
            pybind11::array_t<position_t> _ann_start,
            pybind11::array_t<position_t> _ann_end) {
          if (_gene_idx.ndim() != 1
              || _start.ndim() != 1
              || _end.ndim() != 1
              || _ann_start.ndim() != 1
              || _ann_end.ndim() != 1) {
            throw std::invalid_argument("Exons::init input arrays must be 1D");
          }
          if (_gene_idx.shape(0) != _start.shape(0)
              || _gene_idx.shape(0) != _end.shape(0)
              || _gene_idx.shape(0) != _ann_start.shape(0)
              || _gene_idx.shape(0) != _ann_end.shape(0)) {
            throw std::invalid_argument(
                "Exons::init input arrays have conflicting sizes");
          }
          // unchecked accesses to numpy array
          auto gene_idx = _gene_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto ann_start = _ann_start.unchecked<1>();
          auto ann_end = _ann_end.unchecked<1>();
          // create vector of genes matching input arrays
          std::vector<majiq::Exon> exon_vec{};
          exon_vec.reserve(gene_idx.shape(0));
          for (pybind11::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
            exon_vec.push_back(majiq::Exon{
                majiq::KnownGene{gene_idx(i), genes},
                majiq::ExonIntervalT{start(i), end(i)},
                majiq::ExonIntervalT{ann_start(i), ann_end(i)}});
          }
          pybind11::gil_scoped_release release;  // release GIL at this stage
          return std::make_shared<Exons>(genes, std::move(exon_vec));
        }),
        "Create Exons object using Genes and info about each exon",
        pybind11::arg("genes"), pybind11::arg("gene_idx"),
        pybind11::arg("start"), pybind11::arg("end"),
        pybind11::arg("annotated_start"), pybind11::arg("annotated_end"))
    .def("checksum",
        [](const Exons& self) { return majiq::checksum(self); },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "checksum of exons")
    .def("get_annotated", &Exons::get_annotated,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Return original annotated exons")
    .def("potential_introns",
        [](const std::shared_ptr<Exons>& exons_ptr, bool make_simplified) {
        return majiq::GeneIntrons::PotentialIntrons(exons_ptr, make_simplified);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Return denovo, nonpassed introns corresponding to these exons",
        pybind11::arg("make_simplified"))
    .def_property_readonly("annotated_start",
        [](pybind11::object& exons_obj) -> pybind11::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, data.start);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of annotated exon starts")
    .def_property_readonly("annotated_end",
        [](pybind11::object& exons_obj) -> pybind11::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, data.end);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of annotated exon ends")
    .def("is_denovo",
        [](const Exons& self,
          pybind11::array_t<size_t> exon_idx) -> pybind11::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_denovo(); };
        return pybind11::vectorize(f)(exon_idx);
        },
        "Indicate if selected exons are denovo",
        pybind11::arg("exon_idx"))
    .def("is_exon_extension",
        [](const Exons& self,
          pybind11::array_t<size_t> exon_idx) -> pybind11::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_exon_extension(); };
        return pybind11::vectorize(f)(exon_idx);
        },
        "Indicate if selected exons have exon extension",
        pybind11::arg("exon_idx"))
    .def("is_full_exon",
        [](const Exons& self,
          pybind11::array_t<size_t> exon_idx) -> pybind11::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_full_exon(); };
        return pybind11::vectorize(f)(exon_idx);
        },
        "Indicate if selected exons are full exons",
        pybind11::arg("exon_idx"))
    .def("is_half_exon",
        [](const Exons& self,
          pybind11::array_t<size_t> exon_idx) -> pybind11::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_half_exon(); };
        return pybind11::vectorize(f)(exon_idx);
        },
        "Indicate if selected exons are half exons",
        pybind11::arg("exon_idx"));
}


}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYEXONS_HPP
