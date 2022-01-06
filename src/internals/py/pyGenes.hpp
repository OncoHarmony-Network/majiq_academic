/**
 * pyGenes.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYGENES_HPP
#define MAJIQ_PYBIND_PYGENES_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../Genes.hpp"

#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyGenes_t = pyClassShared_t<majiq::Genes>;

inline void init_Genes(pyGenes_t& pyGenes) {
  define_coordinates_properties(pyGenes);
  pyGenes
    .def("checksum",
        [](const Genes& self) { return majiq::checksum(self); },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "checksum of genes")
    .def(
        pybind11::init([](
            std::shared_ptr<Contigs> contigs,
            pybind11::array_t<size_t> _contig_idx,
            pybind11::array_t<position_t> _start,
            pybind11::array_t<position_t> _end,
            pybind11::array_t<std::array<char, 1>> _strand,
            pybind11::list geneid,
            pybind11::list genename) {
          if (_contig_idx.ndim() != 1
              || _start.ndim() != 1
              || _end.ndim() != 1
              || _strand.ndim() != 1) {
            throw std::invalid_argument("Genes::init input arrays must be 1D");
          }
          if (_contig_idx.shape(0) != _start.shape(0)
              || _contig_idx.shape(0) != _end.shape(0)
              || _contig_idx.shape(0) != _strand.shape(0)
              || static_cast<size_t>(_contig_idx.shape(0)) != geneid.size()
              || static_cast<size_t>(_contig_idx.shape(0)) != genename.size()) {
            throw std::invalid_argument(
                "Genes::init input arrays have conflicting sizes");
          }
          auto contig_idx = _contig_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto strand = _strand.unchecked<1>();
          std::vector<majiq::Gene> gene_vec{};
          gene_vec.reserve(geneid.size());
          for (size_t i = 0; i < geneid.size(); ++i) {
            auto strand_i = static_cast<GeneStrandness>(strand(i)[0]);
            if (strand_i != GeneStrandness::FORWARD
                && strand_i != GeneStrandness::REVERSE) {
              throw std::invalid_argument(
                  "Genes::init input arrays have ambiguous/invalid strand");
            }
            gene_vec.push_back(majiq::Gene{
                majiq::KnownContig{contig_idx(i), contigs},
                majiq::ClosedInterval{start(i), end(i)},
                strand_i,
                geneid[i].cast<majiq::geneid_t>(),
                genename[i].cast<majiq::genename_t>()});
          }
          pybind11::gil_scoped_release release;  // release GIL at this stage
          return Genes::create(contigs, std::move(gene_vec));
        }),
        "Create Genes object using Contigs object and arrays defining genes",
        pybind11::arg("contigs"), pybind11::arg("contig_idx"),
        pybind11::arg("start"), pybind11::arg("end"), pybind11::arg("strand"),
        pybind11::arg("gene_id"), pybind11::arg("gene_name"))
    .def_property_readonly("gene_id", &majiq::Genes::geneids,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Sequence[str] of gene ids in order matching gene_idx")
    .def_property_readonly("gene_name", &majiq::Genes::genenames,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Sequence[str] of gene names in order matching gene_idx")
    .def("__contains__",
        [](const Genes& self, geneid_t x) { return self.count(x) > 0; },
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("__getitem__",
        [](const Genes& self, geneid_t x) { return self.get_idx(x); },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "gene_idx for specified gene_id", pybind11::arg("gene_id"));
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYGENES_HPP
