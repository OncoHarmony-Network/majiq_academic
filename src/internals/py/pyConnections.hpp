/**
 * pyConnections.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYCONNECTIONS_HPP
#define MAJIQ_PYBIND_PYCONNECTIONS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../Exons.hpp"
#include "../GeneIntrons.hpp"
#include "../GeneJunctions.hpp"
#include "../PassedIntrons.hpp"

#include "constants.hpp"
#include "coordinatesProperties.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pyGeneIntrons_t = pyClassShared_t<GeneIntrons>;
using pyGeneJunctions_t = pyClassShared_t<GeneJunctions>;

template <
  typename RegionsT,
  typename RegionT = base_t<decltype(std::declval<RegionsT>().data()[0])>,
  typename IntervalT = decltype(std::declval<RegionT>().coordinates)>
inline void define_connections_properties(
        pyClassShared_t<RegionsT>& pyRegions) {
  pyRegions
    .def(
        pybind11::init([](
            std::shared_ptr<Genes> genes,
            pybind11::array_t<size_t> _gene_idx,
            pybind11::array_t<position_t> _start,
            pybind11::array_t<position_t> _end,
            pybind11::array_t<bool> _denovo,
            pybind11::array_t<bool> _passed_build,
            pybind11::array_t<bool> _simplified) {
          if (_gene_idx.ndim() != 1
              || _start.ndim() != 1
              || _end.ndim() != 1
              || _denovo.ndim() != 1
              || _passed_build.ndim() != 1
              || _simplified.ndim() != 1) {
            throw std::invalid_argument(
                "Connections init input arrays must be 1D");
          }
          if (_gene_idx.shape(0) != _start.shape(0)
              || _gene_idx.shape(0) != _end.shape(0)
              || _gene_idx.shape(0) != _denovo.shape(0)
              || _gene_idx.shape(0) != _passed_build.shape(0)
              || _gene_idx.shape(0) != _simplified.shape(0)) {
            throw std::invalid_argument(
                "Connections init input arrays have conflicting sizes");
          }
          auto gene_idx = _gene_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto denovo = _denovo.unchecked<1>();
          auto passed_build = _passed_build.unchecked<1>();
          auto simplified = _simplified.unchecked<1>();
          // Create vector of the regions to add matching input arrays
          std::vector<RegionT> connection_vec{};
          connection_vec.reserve(gene_idx.shape(0));
          for (pybind11::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
            connection_vec.emplace_back(
                    KnownGene{gene_idx(i), genes},
                    IntervalT{start(i), end(i)},
                    denovo(i), passed_build(i), simplified(i));
          }
          return std::make_shared<RegionsT>(genes, std::move(connection_vec));
        }),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Create connections using Genes and arrays defining each connection",
        pybind11::arg("genes"),
        pybind11::arg("gene_idx"), pybind11::arg("start"), pybind11::arg("end"),
        pybind11::arg("denovo"), pybind11::arg("passed_build"),
        pybind11::arg("simplified"))
    .def("checksum",
        [](const RegionsT& self) {
        return self.template checksum<true>();
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "checksum of connections including data (connections, passed, etc)")
    .def("checksum_nodata",
        [](const RegionsT& self) {
        return self.template checksum<false>();
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "checksum of connections gene/coordinates")
    .def("_pass_all", &RegionsT::pass_all,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "pass all connections")
    .def("_simplify_all", &RegionsT::simplify_all,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "simplify all connections")
    .def("_unsimplify_all", &RegionsT::unsimplify_all,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "unsimplify all connections")
    .def_property_readonly("denovo",
        [](pybind11::object& regions_obj) -> pybind11::array_t<bool> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, data.denovo);
        return ArrayFromVectorAndOffset<bool, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[bool] indicating if connection was not found in annotations")
    .def_property_readonly("passed_build",
        [](pybind11::object& regions_obj) -> pybind11::array_t<bool> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, data.passed_build);
        return ArrayFromVectorAndOffset<bool, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[bool] indicating if passed build criteria to be in LSV")
    .def_property_readonly("simplified",
        [](pybind11::object& regions_obj) -> pybind11::array_t<bool> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, data.simplified);
        return ArrayFromVectorAndOffset<bool, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "array[bool] indicating if the connection is simplified")
    .def("connect_exons", &RegionsT::connect_exons,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Connect regions to specified exons, updataing {start,end}_exon_idx",
        pybind11::arg("exons"))
    .def_property_readonly("connected_exons",
        [](RegionsT& self) -> std::optional<std::shared_ptr<Exons>> {
        using return_t = std::optional<std::shared_ptr<Exons>>;
        return self.is_connected()
          ? return_t{self.connected_exons()} : return_t{};
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Exons connected to (or None if not connected to any exons)")
    .def("src_exon_idx",
        [](const RegionsT& self, pybind11::array_t<size_t> region_idx)
        -> pybind11::array_t<size_t> {
        auto f = [&self](size_t i) -> size_t {
          if (i >= self.size()) {
            throw std::invalid_argument("region_idx has values out of range");
          }
          return self[i].src_exon_idx(); };
        return pybind11::vectorize(f)(region_idx);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        array[int] indicating exon_idx for connection source exon

        Note: Uninitialized values default to all 0.
        )pbdoc",
        pybind11::arg("region_idx"))
    .def("dst_exon_idx",
        [](const RegionsT& self, pybind11::array_t<size_t> region_idx)
        -> pybind11::array_t<size_t> {
        auto f = [&self](size_t i) -> size_t {
          if (i >= self.size()) {
            throw std::invalid_argument("region_idx has values out of range");
          }
          return self[i].dst_exon_idx(); };
        return pybind11::vectorize(f)(region_idx);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        array[int] indicating exon_idx for connection target exon

        Note: Uninitialized values default to all 0.
        )pbdoc",
        pybind11::arg("region_idx"))
    .def_property_readonly("start_exon_idx",
        [](pybind11::object& regions_obj) -> pybind11::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, data.start_exon_idx);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        array[int] indicating exon_idx for connection start

        Note: Uninitialized values default to all 0.
        )pbdoc")
    .def_property_readonly("end_exon_idx",
        [](pybind11::object& regions_obj) -> pybind11::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, data.end_exon_idx);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        array[int] indicating exon_idx for connection end

        Note: Uninitialized values default to all 0.
        )pbdoc");
}

inline void init_GeneIntrons(pyGeneIntrons_t& pyGeneIntrons) {
  define_coordinates_properties(pyGeneIntrons);
  define_connections_properties(pyGeneIntrons);
  pyGeneIntrons
    .def("update_flags_from", &GeneIntrons::UpdateFlagsFrom,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Update intron flags using introns that overlap from input",
        pybind11::arg("donor_introns"))
    .def("build_group", [](std::shared_ptr<GeneIntrons>& gene_introns) {
        return GroupIntronsGenerator(gene_introns);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Create build group to update passed introns in place")
    .def("filter_passed", &GeneIntrons::FilterPassed,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        Get subset of introns that passed build filters

        Parameters
        ----------
        keep_annotated: bool
            Keep all annotated introns regardless of whether they passed
        discard_denovo: bool
            Discard all denovo introns regardless of whether they passed
        )pbdoc",
        pybind11::arg("keep_annotated") = DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        pybind11::arg("discard_denovo") = !DEFAULT_BUILD_DENOVO_IR)
    .def("__repr__", [](const GeneIntrons& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneIntrons<" << self.size() << " total>";
        return oss.str();
        },
        pybind11::call_guard<pybind11::gil_scoped_release>());
}

inline void init_GeneJunctions(pyGeneJunctions_t& pyGeneJunctions) {
  define_coordinates_properties(pyGeneJunctions);
  define_connections_properties(pyGeneJunctions);
  pyGeneJunctions
    .def("__repr__", [](const GeneJunctions& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneJunctions<" << self.size() << " total>";
        return oss.str();
        },
        pybind11::call_guard<pybind11::gil_scoped_release>());
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYCONNECTIONS_HPP
