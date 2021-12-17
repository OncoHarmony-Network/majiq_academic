/**
 * coordinatesProperties.hpp
 *
 * Add coordinates properties to pyClassShared_t
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PYBIND_COORDINATESPROPERTIES_HPP
#define MAJIQ_PYBIND_COORDINATESPROPERTIES_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <array>
#include <memory>

#include "../Genes.hpp"
#include "../Meta.hpp"

#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

template <
  typename RegionsT,
  typename RegionT = base_t<decltype(std::declval<RegionsT>().data()[0])>,
  typename ParentT = base_t<decltype(std::declval<RegionT>().parent())>,
  typename IntervalT = decltype(std::declval<RegionT>().coordinates)>
inline void define_coordinates_properties(pyClassShared_t<RegionsT>& pyRegions) {
  // general properties
  pyRegions
    .def("__len__",
        &RegionsT::size, pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("__eq__", [](const RegionsT& x, const RegionsT& y) { return x == y; },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        pybind11::is_operator())
    .def_property_readonly("_parents", &RegionsT::parents,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get parents object on which regions are defined (e.g. contigs, genes)")
    .def_property_readonly("_parent_idx_start",
        [](pybind11::object& regions_obj) {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        return ArrayFromOffsetsVector<size_t>(
            regions.parent_idx_offsets(), true, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "First index into regions corresponding to associated parent")
    .def_property_readonly("_parent_idx_end",
        [](pybind11::object& regions_obj) {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        return ArrayFromOffsetsVector<size_t>(
            regions.parent_idx_offsets(), false, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "One after last index into regions corresponding to associated parent")
    .def_property_readonly("start",
        [](pybind11::object& regions_obj) -> pybind11::array_t<position_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of starts for each feature")
    .def_property_readonly("end",
        [](pybind11::object& regions_obj) -> pybind11::array_t<position_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of ends for each feature");
  // if has contigs field
  if constexpr(majiq::detail::has_contig_field<RegionT>::value) {
    pyRegions.def_property_readonly("contig_idx",
        [](pybind11::object& regions_obj) -> pybind11::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, contig.idx_);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of indexes indicating contig feature belongs to");
  }
  // if has strand field
  if constexpr(majiq::detail::has_strand_field<RegionT>::value) {
    pyRegions.def_property_readonly("strand",
        [](pybind11::object& regions_obj)
        -> pybind11::array_t<std::array<char, 1>> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, strand);
        return ArrayFromVectorAndOffset<std::array<char, 1>, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[char] of characters indicating strand of each feature");
  }
  // if it has the gene field
  if constexpr(majiq::detail::has_gene_field<RegionT>::value) {
    pyRegions.def_property_readonly("gene_idx",
        [](pybind11::object& regions_obj) -> pybind11::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, gene.idx_);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[int] of indexes indicating gene feature belongs to");
  }
  // has genes as parent, enable quick lookup of specific intervals
  if constexpr(std::is_same_v<ParentT, majiq::KnownGene>) {
    pyRegions
      .def("index",
          [](const RegionsT& self,
            pybind11::array_t<size_t> gene_idx,
            pybind11::array_t<position_t> start,
            pybind11::array_t<position_t> end) {
          auto f = [&self](
              size_t g, position_t s, position_t e) -> std::ptrdiff_t {
            if (g >= self.parents_->size()) { return -1; }
            IntervalT iv;
            try {
              iv = IntervalT{s, e};
            } catch (std::invalid_argument& e) {
              return -1;
            }
            auto it = self.find(RegionT{(*self.parents_)[g], iv});
            return it == self.end() ? -1 : it - self.begin(); };
          return pybind11::vectorize(f)(gene_idx, start, end);
          },
          pybind11::call_guard<pybind11::gil_scoped_release>(),
          "Get indexes for specified regions (or -1 if it doesn't exist)",
          pybind11::arg("gene_idx"),
          pybind11::arg("start"),
          pybind11::arg("end"));
  }
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_COORDINATESPROPERTIES_HPP
