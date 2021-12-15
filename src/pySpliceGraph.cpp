/**
 * pySpliceGraph.cpp
 *
 * Set up python bindings to contigs
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <memory>
#include <optional>
#include <cstddef>
#include <stdexcept>

#include "pybind_utils.hpp"
#include "internals/Contigs.hpp"
#include "internals/Genes.hpp"
#include "internals/Exons.hpp"
#include "internals/GeneJunctions.hpp"
#include "internals/GeneIntrons.hpp"
#include "internals/SpliceGraph.hpp"
#include "internals/GFF3.hpp"
#include "internals/SJJunctions.hpp"
#include "internals/SJIntrons.hpp"
#include "internals/SJBinsReads.hpp"
#include "internals/PassedJunctions.hpp"
#include "internals/PassedIntrons.hpp"
#include "internals/ExonConnections.hpp"
#include "internals/Events.hpp"
#include "internals/EventsCoverage.hpp"
#include "internals/EventsAlign.hpp"
#include "internals/SpliceGraphReads.hpp"
#include "internals/SimplifierGroup.hpp"
#include "internals/GeneJunctionsAccumulator.hpp"
#include "internals/Meta.hpp"

#include "internals/ExperimentThresholds.hpp"

#include <majiqinclude/ResourcePool.hpp>

// NOTE: right now we have a global PRNG, which is not threadsafe, but
// multithreading is really only used by htslib internally at this time
static MajiqInclude::ResourcePool<majiq::rng_t> global_rng_pool{};

// default value of ExperimentThresholds
static const auto DEFAULT_THRESHOLDS = majiq::ExperimentThresholds(
    DEFAULT_BUILD_MINREADS, DEFAULT_BUILD_MINDENOVO, DEFAULT_BUILD_MINPOS,
    DEFAULT_BUILD_MAX_PCTBINS, DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
    DEFAULT_BUILD_MATCH_INTRON_PROBABILITY);

namespace py = pybind11;
template <typename T>
using pyClassShared_t = py::class_<T, std::shared_ptr<T>>;
using pyContigs_t = pyClassShared_t<majiq::Contigs>;
using pyGenes_t = pyClassShared_t<majiq::Genes>;
using pyExons_t = pyClassShared_t<majiq::Exons>;
using pyGeneIntrons_t = pyClassShared_t<majiq::GeneIntrons>;
using pyGeneJunctions_t = pyClassShared_t<majiq::GeneJunctions>;
using pySJIntrons_t = pyClassShared_t<majiq::SJIntrons>;
using pySJJunctions_t = pyClassShared_t<majiq::SJJunctions>;
using pySJJunctionsBins_t = pyClassShared_t<majiq::SJJunctionsBins>;
using pyGroupJunctionsGen_t = pyClassShared_t<majiq::GroupJunctionsGenerator>;
using pyPassedJunctionsGen_t = pyClassShared_t<majiq::PassedJunctionsGenerator>;
using pyGroupIntronsGen_t = pyClassShared_t<majiq::GroupIntronsGenerator>;
using pySJIntronsBins_t = pyClassShared_t<majiq::SJIntronsBins>;
using pyExonConnections_t = pyClassShared_t<majiq::ExonConnections>;
using pySimplifierGroup_t = pyClassShared_t<majiq::SimplifierGroup>;
using pyEvents_t = pyClassShared_t<majiq::Events>;
using pyEventsCoverage_t = pyClassShared_t<majiq::EventsCoverage>;
using pyEventsAlign_t = pyClassShared_t<majiq::EventsAlign>;
using pySpliceGraphReads_t = pyClassShared_t<majiq::SpliceGraphReads>;
using pyGeneJunctionsAccumulator_t = pyClassShared_t<majiq::GeneJunctionsAccumulator>;

using pyExperimentThresholds_t = pyClassShared_t<majiq::ExperimentThresholds>;
using pyIntronThresholdsGenerator_t
  = pyClassShared_t<majiq::IntronThresholdsGenerator>;

template <typename T>
using base_t = std::remove_cv_t<std::remove_reference_t<T>>;

template <typename SJBinsT,
         typename RegionsT
           = base_t<decltype(*(std::declval<SJBinsT>().regions()))>,
         typename BinReads
           = base_t<decltype(std::declval<SJBinsT>().reads()[0])>,
         typename CountT = decltype(std::declval<BinReads>().bin_reads)>
void define_sjbins_properties(pyClassShared_t<SJBinsT>& pySJBins) {
  using majiq_pybind::ArrayFromVectorAndOffset;
  using majiq_pybind::ArrayFromOffsetsVector;
  pySJBins
    .def_property_readonly("_regions", &SJBinsT::regions,
        py::call_guard<py::gil_scoped_release>(),
        "Underlying regions bin reads are defined over")
    .def_property_readonly("total_bins", &SJBinsT::total_bins,
        py::call_guard<py::gil_scoped_release>(),
        "the total number of bins possible (positions for junctions, too)")
    .def_property_readonly("bin_reads",
        [](py::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        const size_t offset = offsetof(BinReads, bin_reads);
        return ArrayFromVectorAndOffset<CountT, BinReads>(
            sj.reads(), offset, sj_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of reads for each bin")
    .def_property_readonly("bin_idx",
        [](py::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        const size_t offset = offsetof(BinReads, bin_idx);
        return ArrayFromVectorAndOffset<majiq::junction_pos_t, BinReads>(
            sj.reads(), offset, sj_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Bin index for the each bin reads")
    .def_property_readonly("_offsets",
        [](py::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            sj.offsets(), 0, sj_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Raw offsets for regions into bin reads")
    .def("numstacks",
        [](const SJBinsT& self,
          py::array_t<size_t> idx, py::array_t<majiq::real_t> pvalue) {
        auto f = [&self](size_t i, majiq::real_t p) {
          if (i >= self.num_regions()) {
            throw std::invalid_argument("idx has out of range values");
          }
          return self.numstacks(i, p); };
        return py::vectorize(f)(idx, pvalue);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Get number of stacks for the specified regions given threshold",
        py::arg("region_idx"),
        py::arg("pvalue_threshold") = DEFAULT_BUILD_STACK_PVALUE)
    .def("numbins",
        [](const SJBinsT& self,
          py::array_t<size_t> idx,
          py::array_t<CountT> minreads) {
        auto f = [&self](size_t i, CountT r) {
          if (i >= self.num_regions()) {
            throw std::invalid_argument("idx has out of range values");
          }
          return self.numbins_minreads(i, r); };
        return py::vectorize(f)(idx, minreads);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of bins for regions with more than specified number of reads",
        py::arg("region_idx"), py::arg("minreads"))
    .def("numreads",
        [](const SJBinsT& self,
          py::array_t<size_t> idx,
          py::array_t<majiq::junction_pos_t> num_stacks) {
        auto f = [&self](size_t i, majiq::junction_pos_t n) {
          if (i >= self.num_regions()) {
            throw std::invalid_argument("idx has out of range values");
          }
          return self.numreads(i, n); };
        return py::vectorize(f)(idx, num_stacks);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of reads for regions given known number of stacks",
        py::arg("region_idx"), py::arg("num_stacks"))
    .def("__len__", &SJBinsT::size,
        py::call_guard<py::gil_scoped_release>(),
        "Total number of bin reads")
    .def(py::init([](
            std::shared_ptr<RegionsT> regions,
            py::array_t<CountT> _bin_reads,
            py::array_t<majiq::junction_pos_t> _bin_idx,
            py::array_t<size_t> _offsets,
            majiq::junction_pos_t total_bins) {
          auto check_1d = [](const auto& x) {
            if (x.ndim() != 1)
            throw std::runtime_error("Bins arrays must be 1D"); };
          check_1d(_bin_reads);
          check_1d(_bin_idx);
          check_1d(_offsets);
          std::vector<size_t> offsets_vec(_offsets.shape(0));
          {
            auto offsets = _offsets.unchecked<1>();
            for (size_t i = 0; i < offsets_vec.size(); ++i) {
              offsets_vec[i] = offsets(i);
            }
          }
          if (_bin_reads.shape(0) != _bin_idx.shape(0)) {
            throw std::runtime_error(
                "bin_reads and bin_idx should be same length");
          }
          std::vector<BinReads> br_vec(_bin_reads.shape(0));
          {
            auto bin_reads = _bin_reads.template unchecked<1>();
            auto bin_idx = _bin_idx.unchecked<1>();
            for (size_t i = 0; i < br_vec.size(); ++i) {
              br_vec[i] = BinReads{bin_idx(i), bin_reads(i)};
            }
          }
          return SJBinsT{regions,
            std::move(br_vec), std::move(offsets_vec), total_bins};
          }),
        py::call_guard<py::gil_scoped_release>(),
        "Initialize bins over specified regions with per-bin coverage",
        py::arg("sj"),
        py::arg("bin_reads"), py::arg("bin_idx"), py::arg("_offsets"),
        py::arg("total_bins"));
}

template <typename RegionsT,
         typename RegionT
           = base_t<decltype(std::declval<RegionsT>().data()[0])>,
         typename ParentT = base_t<decltype(std::declval<RegionT>().parent())>,
         typename IntervalT = decltype(std::declval<RegionT>().coordinates)>
void define_coordinates_properties(pyClassShared_t<RegionsT>& pyRegions) {
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  using majiq_pybind::ArrayFromOffsetsVector;
  if constexpr(std::is_same_v<ParentT, majiq::KnownGene>) {
    pyRegions
      .def("index",
          [](const RegionsT& self,
            py::array_t<size_t> gene_idx,
            py::array_t<position_t> start,
            py::array_t<position_t> end) {
          auto f = [&self](size_t g, position_t s, position_t e) -> std::ptrdiff_t {
            if (g >= self.parents_->size()) { return -1; }
            IntervalT iv;
            try {
              iv = IntervalT{s, e};
            } catch (std::invalid_argument& e) {
              return -1;
            }
            auto it = self.find(RegionT{(*self.parents_)[g], iv});
            return it == self.end() ? -1 : it - self.begin(); };
          return py::vectorize(f)(gene_idx, start, end);
          },
          py::call_guard<py::gil_scoped_release>(),
          "Get indexes for specified regions (or -1 if it doesn't exist)",
          py::arg("gene_idx"),
          py::arg("start"),
          py::arg("end"));
  }
  pyRegions
    .def("__len__", &RegionsT::size, py::call_guard<py::gil_scoped_release>())
    .def("__eq__", [](const RegionsT& x, const RegionsT& y) { return x == y; },
        py::call_guard<py::gil_scoped_release>(),
        py::is_operator())
    .def_property_readonly("_parents", &RegionsT::parents,
        py::call_guard<py::gil_scoped_release>(),
        "Get parents object on which regions are defined (e.g. contigs, genes)")
    .def_property_readonly("_parent_idx_start",
        [](py::object& regions_obj) {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        return ArrayFromOffsetsVector<size_t>(
            regions.parent_idx_offsets(), true, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "First index into regions corresponding to associated parent")
    .def_property_readonly("_parent_idx_end",
        [](py::object& regions_obj) {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        return ArrayFromOffsetsVector<size_t>(
            regions.parent_idx_offsets(), false, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "One after last index into regions corresponding to associated parent")
    .def_property_readonly("start",
        [](py::object& regions_obj) -> py::array_t<position_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[int] of starts for each feature")
    .def_property_readonly("end",
        [](py::object& regions_obj) -> py::array_t<position_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[int] of ends for each feature");
  if constexpr(majiq::detail::has_contig_field<RegionT>::value) {
    pyRegions.def_property_readonly("contig_idx",
        [](py::object& regions_obj) -> py::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, contig.idx_);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[int] of indexes indicating contig feature belongs to");
  }
  if constexpr(majiq::detail::has_strand_field<RegionT>::value) {
    pyRegions.def_property_readonly("strand",
        [](py::object& regions_obj) -> py::array_t<std::array<char, 1>> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, strand);
        return ArrayFromVectorAndOffset<std::array<char, 1>, RegionT>(
            regions.data(), offset, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[char] of characters indicating strand of each feature");
  }
  if constexpr(majiq::detail::has_gene_field<RegionT>::value) {
    pyRegions.def_property_readonly("gene_idx",
        [](py::object& regions_obj) -> py::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, gene.idx_);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[int] of indexes indicating gene feature belongs to");
  }
  if constexpr(
      std::is_same_v<decltype(std::declval<RegionT>().data),
                     majiq::detail::ConnectionData>) {
    pyRegions
      .def("checksum",
          [](const RegionsT& self) {
          return self.template checksum<true>();
          },
          py::call_guard<py::gil_scoped_release>(),
          "checksum of connections including data (connections, passed, etc)")
      .def("checksum_nodata",
          [](const RegionsT& self) {
          return self.template checksum<false>();
          },
          py::call_guard<py::gil_scoped_release>(),
          "checksum of connections gene/coordinates")
      .def("_pass_all", &RegionsT::pass_all,
          py::call_guard<py::gil_scoped_release>(), "pass all connections")
      .def("_simplify_all", &RegionsT::simplify_all,
          py::call_guard<py::gil_scoped_release>(), "simplify all connections")
      .def("_unsimplify_all", &RegionsT::unsimplify_all,
          py::call_guard<py::gil_scoped_release>(),
          "unsimplify all connections")
      .def_property_readonly("denovo",
          [](py::object& regions_obj) -> py::array_t<bool> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.denovo);
          return ArrayFromVectorAndOffset<bool, RegionT>(
              regions.data(), offset, regions_obj);
          },
          py::call_guard<py::gil_scoped_release>(),
          "array[bool] indicating if connection was not found in annotations")
      .def_property_readonly("passed_build",
          [](py::object& regions_obj) -> py::array_t<bool> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.passed_build);
          return ArrayFromVectorAndOffset<bool, RegionT>(
              regions.data(), offset, regions_obj);
          },
          py::call_guard<py::gil_scoped_release>(),
          "array[bool] indicating if passed build criteria to be in LSV")
      .def_property_readonly("simplified",
          [](py::object& regions_obj) -> py::array_t<bool> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.simplified);
          return ArrayFromVectorAndOffset<bool, RegionT>(
              regions.data(), offset, regions_obj);
          },
          py::call_guard<py::gil_scoped_release>(),
          "array[bool] indicating if the connection is simplified")
      .def("connect_exons", &RegionsT::connect_exons,
          py::call_guard<py::gil_scoped_release>(),
          "Connect regions to specified exons, updataing {start,end}_exon_idx",
          py::arg("exons"))
      .def_property_readonly("connected_exons",
          [](RegionsT& self) -> std::optional<std::shared_ptr<majiq::Exons>> {
          using return_t = std::optional<std::shared_ptr<majiq::Exons>>;
          return self.is_connected()
            ? return_t{self.connected_exons()} : return_t{};
          },
          py::call_guard<py::gil_scoped_release>(),
          "Exons connected to (or None if not connected to any exons)")
      .def("src_exon_idx",
          [](const RegionsT& self, py::array_t<size_t> region_idx) -> py::array_t<size_t> {
          auto f = [&self](size_t i) -> size_t {
            if (i >= self.size()) {
              throw std::invalid_argument("region_idx has values out of range");
            }
            return self[i].src_exon_idx(); };
          return py::vectorize(f)(region_idx);
          },
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
          array[int] indicating exon_idx for connection source exon

          Note: Uninitialized values default to all 0.
          )pbdoc",
          py::arg("region_idx"))
      .def("dst_exon_idx",
          [](const RegionsT& self, py::array_t<size_t> region_idx) -> py::array_t<size_t> {
          auto f = [&self](size_t i) -> size_t {
            if (i >= self.size()) {
              throw std::invalid_argument("region_idx has values out of range");
            }
            return self[i].dst_exon_idx(); };
          return py::vectorize(f)(region_idx);
          },
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
          array[int] indicating exon_idx for connection target exon

          Note: Uninitialized values default to all 0.
          )pbdoc",
          py::arg("region_idx"))
      .def_property_readonly("start_exon_idx",
          [](py::object& regions_obj) -> py::array_t<size_t> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.start_exon_idx);
          return ArrayFromVectorAndOffset<size_t, RegionT>(
              regions.data(), offset, regions_obj);
          },
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
          array[int] indicating exon_idx for connection start

          Note: Uninitialized values default to all 0.
          )pbdoc")
      .def_property_readonly("end_exon_idx",
          [](py::object& regions_obj) -> py::array_t<size_t> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.end_exon_idx);
          return ArrayFromVectorAndOffset<size_t, RegionT>(
              regions.data(), offset, regions_obj);
          },
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
          array[int] indicating exon_idx for connection end

          Note: Uninitialized values default to all 0.
          )pbdoc");
  }
  return;
}

struct ConnectionsArrays {
  py::array_t<size_t> gene_idx;
  py::array_t<majiq::position_t> start;
  py::array_t<majiq::position_t> end;
  py::array_t<bool> denovo;
  py::array_t<bool> passed_build;
  py::array_t<bool> simplified;
};

template <class Connections>
std::shared_ptr<Connections> MakeConnections(
    std::shared_ptr<majiq::Genes> genes, const ConnectionsArrays& values) {
  // unchecked accesses to numpy array
  auto gene_idx = values.gene_idx.unchecked<1>();
  auto start = values.start.unchecked<1>();
  auto end = values.end.unchecked<1>();
  auto denovo = values.denovo.unchecked<1>();
  auto passed_build = values.passed_build.unchecked<1>();
  auto simplified = values.simplified.unchecked<1>();
  // create vector of genes matching input arrays for specified connections
  using IntervalT = typename Connections::IntervalT;
  using RegionT = typename Connections::value_type;
  std::vector<RegionT> connection_vec{};
  connection_vec.reserve(gene_idx.shape(0));
  for (py::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
    if (gene_idx(i) >= genes->size()) {
      throw std::invalid_argument("gene_idx has values out of range");
    }
    connection_vec.push_back(RegionT{
        majiq::KnownGene{gene_idx(i), genes}, IntervalT{start(i), end(i)},
        denovo(i), passed_build(i), simplified(i)});
  }
  // create shared pointer to connections object
  return std::make_shared<Connections>(genes, std::move(connection_vec));
}

void init_Contigs(pyContigs_t& pyContigs) {
  using majiq::seqid_t;
  using majiq::Contigs;
  using majiq::Contig;
  pyContigs
    .def("checksum",
        [](const Contigs& self) { return majiq::checksum(self); },
        py::call_guard<py::gil_scoped_release>(),
        "checksum of contigs")
    .def(
        py::init([](py::list seqids) {
          auto result = Contigs::create();
          for (auto seqid : seqids) {
            result->add(Contig{seqid.cast<seqid_t>()});
          }
          return result;
        }),
        py::call_guard<py::gil_scoped_release>(),
        "Set up Contigs object using specified identifiers",
        py::arg("seqids"))
    .def_property_readonly("seqid", &Contigs::seqids,
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("__repr__", [](const Contigs& self) -> std::string {
        std::ostringstream oss;
        oss << self;
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>())
    .def("__len__", &Contigs::size, py::call_guard<py::gil_scoped_release>())
    .def("__contains__",
        [](const Contigs& s, seqid_t x) -> bool { return s.count(x) > 0; },
        py::call_guard<py::gil_scoped_release>())
    .def("__getitem__",
        [](const Contigs& self, seqid_t x) { return self.get_idx(x); },
        py::call_guard<py::gil_scoped_release>(),
        "contig_idx for specified seqid", py::arg("seqid"));
}

void init_Genes(pyGenes_t& pyGenes) {
  using majiq_pybind::ArrayFromVectorAndOffset;
  using majiq::Contigs;
  using majiq::Genes;
  using majiq::position_t;
  using majiq::geneid_t;
  define_coordinates_properties(pyGenes);
  pyGenes
    .def("checksum",
        [](const Genes& self) { return majiq::checksum(self); },
        py::call_guard<py::gil_scoped_release>(),
        "checksum of genes")
    .def(
        py::init([](
            std::shared_ptr<Contigs> contigs,
            py::array_t<size_t> _contig_idx,
            py::array_t<position_t> _start,
            py::array_t<position_t> _end,
            py::array_t<std::array<char, 1>> _strand,
            py::list geneid,
            py::list genename) {
          auto contig_idx = _contig_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto strand = _strand.unchecked<1>();
          std::vector<majiq::Gene> gene_vec{};
          gene_vec.reserve(geneid.size());
          for (size_t i = 0; i < geneid.size(); ++i) {
            if (contig_idx(i) >= contigs->size()) {
              throw std::invalid_argument("contig_idx has values out of range");
            }
            gene_vec.push_back(majiq::Gene{
                majiq::KnownContig{contig_idx(i), contigs},
                majiq::ClosedInterval{start(i), end(i)},
                static_cast<majiq::GeneStrandness>(strand(i)[0]),
                geneid[i].cast<majiq::geneid_t>(),
                genename[i].cast<majiq::genename_t>()});
          }
          return Genes::create(contigs, std::move(gene_vec));
        }),
        py::call_guard<py::gil_scoped_release>(),
        "Create Genes object using Contigs object and arrays defining genes",
        py::arg("contigs"), py::arg("contig_idx"),
        py::arg("start"), py::arg("end"), py::arg("strand"),
        py::arg("gene_id"), py::arg("gene_name"))
    .def_property_readonly("gene_id", &majiq::Genes::geneids,
        py::call_guard<py::gil_scoped_release>(),
        "Sequence[str] of gene ids in order matching gene_idx")
    .def_property_readonly("gene_name", &majiq::Genes::genenames,
        py::call_guard<py::gil_scoped_release>(),
        "Sequence[str] of gene names in order matching gene_idx")
    .def("__repr__", [](const majiq::Genes& self) -> std::string {
        std::ostringstream oss;
        oss << "Genes<" << self.size() << " total>";
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>())
    .def("__contains__",
        [](const Genes& self, geneid_t x) { return self.count(x) > 0; },
        py::call_guard<py::gil_scoped_release>())
    .def("__getitem__",
        [](const Genes& self, geneid_t x) { return self.get_idx(x); },
        py::call_guard<py::gil_scoped_release>(),
        "gene_idx for specified gene_id", py::arg("gene_id"));
}

void init_Exons(pyExons_t& pyExons) {
  using majiq::Exons;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties(pyExons);
  pyExons
    .def(py::init([](
            std::shared_ptr<majiq::Genes> genes,
            py::array_t<size_t> _gene_idx,
            py::array_t<position_t> _start,
            py::array_t<position_t> _end,
            py::array_t<position_t> _ann_start,
            py::array_t<position_t> _ann_end) {
          // unchecked accesses to numpy array
          auto gene_idx = _gene_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto ann_start = _ann_start.unchecked<1>();
          auto ann_end = _ann_end.unchecked<1>();
          // create vector of genes matching input arrays
          std::vector<majiq::Exon> exon_vec{};
          exon_vec.reserve(gene_idx.shape(0));
          for (py::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
            if (gene_idx(i) >= genes->size()) {
              throw std::invalid_argument("gene_idx has values out of range");
            }
            exon_vec.push_back(majiq::Exon{
                majiq::KnownGene{gene_idx(i), genes},
                majiq::ClosedInterval{start(i), end(i)},
                majiq::ClosedInterval{ann_start(i), ann_end(i)}});
          }
          return std::make_shared<Exons>(genes, std::move(exon_vec));
        }),
        py::call_guard<py::gil_scoped_release>(),
        "Create Exons object using Genes and info about each exon",
        py::arg("genes"), py::arg("gene_idx"),
        py::arg("start"), py::arg("end"),
        py::arg("annotated_start"), py::arg("annotated_end"))
    .def("checksum",
        [](const Exons& self) { return majiq::checksum(self); },
        py::call_guard<py::gil_scoped_release>(),
        "checksum of exons")
    .def("potential_introns",
        [](const std::shared_ptr<Exons>& exons_ptr, bool make_simplified) {
        return majiq::GeneIntrons::PotentialIntrons(exons_ptr, make_simplified);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Return denovo, nonpassed introns corresponding to these exons",
        py::arg("make_simplified"))
    .def_property_readonly("annotated_start",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, data.start);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[int] of annotated exon starts")
    .def_property_readonly("annotated_end",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, data.end);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[int] of annotated exon ends")
    .def("is_denovo",
        [](const Exons& self,
          py::array_t<size_t> exon_idx) -> py::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_denovo(); };
        return py::vectorize(f)(exon_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if selected exons are denovo",
        py::arg("exon_idx"))
    .def("is_exon_extension",
        [](const Exons& self,
          py::array_t<size_t> exon_idx) -> py::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_exon_extension(); };
        return py::vectorize(f)(exon_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if selected exons have exon extension",
        py::arg("exon_idx"))
    .def("is_full_exon",
        [](const Exons& self,
          py::array_t<size_t> exon_idx) -> py::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_full_exon(); };
        return py::vectorize(f)(exon_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if selected exons are full exons",
        py::arg("exon_idx"))
    .def("is_half_exon",
        [](const Exons& self,
          py::array_t<size_t> exon_idx) -> py::array_t<bool> {
        auto f = [&self](size_t i) {
          if (i >= self.size()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self[i].is_half_exon(); };
        return py::vectorize(f)(exon_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if selected exons are half exons",
        py::arg("exon_idx"))
    .def("__repr__", [](const Exons& self) -> std::string {
        std::ostringstream oss;
        oss << "Exons<" << self.size() << " total>";
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>());
}

void init_GeneJunctionsAccumulator(
    pyGeneJunctionsAccumulator_t& pyGeneJunctionsAccumulator) {
  using majiq::Genes;
  using majiq::GeneJunctionsAccumulator;
  using majiq::GeneJunctions;
  pyGeneJunctionsAccumulator
    .def(py::init<const std::shared_ptr<Genes>&>(),
        py::call_guard<py::gil_scoped_release>(),
        "Initialize accumulator of GeneJunctions with Genes",
        py::arg("genes"))
    .def_property_readonly(
        "_genes", &GeneJunctionsAccumulator::genes,
        "Genes used by GeneJunctionsAccumulator")
    .def(
        "add",
        &GeneJunctionsAccumulator::Add,
        "Add/update accumulated junctions",
        py::arg("junctions"),
        py::arg("make_annotated") = false)
    .def(
        "accumulated",
        &GeneJunctionsAccumulator::Accumulated,
        "Return GeneJunctions combining accumulated junctions");
}

void init_GeneJunctions(pyGeneJunctions_t& pyGeneJunctions) {
  using majiq::GeneJunctions;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties(pyGeneJunctions);
  pyGeneJunctions
    .def(py::init([](
            std::shared_ptr<majiq::Genes> genes,
            py::array_t<size_t> gene_idx,
            py::array_t<position_t> start,
            py::array_t<position_t> end,
            py::array_t<bool> denovo,
            py::array_t<bool> passed_build,
            py::array_t<bool> simplified) {
          return MakeConnections<GeneJunctions>(genes, ConnectionsArrays{
              gene_idx, start, end, denovo, passed_build, simplified});
        }),
        py::call_guard<py::gil_scoped_release>(),
        "Create GeneJunctions using Genes and arrays defining each junction",
        py::arg("genes"),
        py::arg("gene_idx"), py::arg("start"), py::arg("end"),
        py::arg("denovo"), py::arg("passed_build"), py::arg("simplified"))
    .def("__repr__", [](const GeneJunctions& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneJunctions<" << self.size() << " total>";
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>());
}

void init_SJIntrons(pySJIntrons_t& pySJIntrons) {
  using majiq::position_t;
  using majiq::Contigs;
  using majiq::SJIntron;
  using majiq::SJIntrons;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties(pySJIntrons);
  pySJIntrons
    .def(py::init([](
            std::shared_ptr<Contigs> contigs,
            py::array_t<size_t> _contig_idx,
            py::array_t<position_t> _start,
            py::array_t<position_t> _end,
            py::array_t<std::array<char, 1>> _strand,
            py::array_t<bool> _annotated) {
          auto contig_idx = _contig_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto strand = _strand.unchecked<1>();
          auto annotated = _annotated.unchecked<1>();
          std::vector<SJIntron> result(start.shape(0));
          for (size_t i = 0; i < result.size(); ++i) {
            if (contig_idx(i) >= contigs->size()) {
              throw std::invalid_argument("contig_idx has values out of range");
            }
            result[i] = SJIntron{
              majiq::KnownContig{contig_idx(i), contigs},
              majiq::ClosedInterval{start(i), end(i)},
              static_cast<majiq::GeneStrandness>(strand(i)[0]),
              annotated(i)};
          }
          return std::make_shared<SJIntrons>(contigs, std::move(result));
          }),
        py::call_guard<py::gil_scoped_release>(),
        "Initialize SJIntrons from contigs and numpy arrays",
        py::arg("contigs"),
        py::arg("contig_idx"),
        py::arg("start"),
        py::arg("end"),
        py::arg("strand"),
        py::arg("annotated"))
    .def_static("from_exons_and_introns", &SJIntrons::FromGeneExonsAndIntrons,
        py::call_guard<py::gil_scoped_release>(),
        "Construct sj introns for input exons/introns",
        py::arg("exons"), py::arg("introns"), py::arg("stranded"))
    .def_property_readonly("annotated",
        [](py::object& introns_obj) -> py::array_t<bool> {
        SJIntrons& introns = introns_obj.cast<SJIntrons&>();
        const size_t offset = offsetof(majiq::SJIntron, data.annotated_);
        return ArrayFromVectorAndOffset<bool, majiq::SJIntron>(
            introns.data(), offset, introns_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "array[bool] indicating if intron is annotated (exon in annotation)")
    .def("__repr__", [](const majiq::SJIntrons& self) -> std::string {
        std::ostringstream oss;
        oss << "SJIntrons<" << self.size() << " total>";
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>());
}

void init_PySpliceGraphReads(pySpliceGraphReads_t& pySpliceGraphReads) {
  using majiq::SpliceGraphReads;
  using majiq::GeneIntrons;
  using majiq::GeneJunctions;
  using majiq::SJIntronsBins;
  using majiq::SJJunctionsBins;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pySpliceGraphReads
    .def_property_readonly("_introns",
        &SpliceGraphReads::introns,
        py::call_guard<py::gil_scoped_release>(),
        "Underlying introns")
    .def_property_readonly("_junctions",
        &SpliceGraphReads::junctions,
        py::call_guard<py::gil_scoped_release>(),
        "Underlying junctions")
    .def_property_readonly("introns_reads",
        [](py::object& self_obj) {
        SpliceGraphReads& self = self_obj.cast<SpliceGraphReads&>();
        return ArrayFromVectorAndOffset<majiq::real_t, majiq::real_t>(
            self.introns_reads(), 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Raw readrates for each intron")
    .def_property_readonly("junctions_reads",
        [](py::object& self_obj) {
        SpliceGraphReads& self = self_obj.cast<SpliceGraphReads&>();
        return ArrayFromVectorAndOffset<majiq::real_t, majiq::real_t>(
            self.junctions_reads(), 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Raw readrates for each junction")
    .def(py::init([](
            const std::shared_ptr<GeneIntrons>& introns,
            const std::shared_ptr<GeneJunctions>& junctions,
            py::array_t<majiq::real_t> _introns_reads,
            py::array_t<majiq::real_t> _junctions_reads) {
          if (_introns_reads.ndim() != 1) {
            throw std::runtime_error("introns_reads must be 1D");
          } else if (_junctions_reads.ndim() != 1) {
            throw std::runtime_error("junctions_reads must be 1D");
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
        py::call_guard<py::gil_scoped_release>(),
        "Initialize SpliceGraphReads from numpy arrays",
        py::arg("introns"),
        py::arg("junctions"),
        py::arg("introns_reads"),
        py::arg("junctions_reads"))
    .def_static("from_sj", &SpliceGraphReads::FromSJ,
        py::call_guard<py::gil_scoped_release>(),
        "Obtain raw readrates for introns/junctions from experiment SJ",
        py::arg("introns"),
        py::arg("junctions"),
        py::arg("sj_introns"),
        py::arg("sj_junctions"));
}

void init_PySimplifierGroup(pySimplifierGroup_t& pySimplifierGroup) {
  using majiq::SimplifierCount;
  using majiq::SimplifierGroup;
  using majiq::ExonConnections;
  using majiq::SpliceGraphReads;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pySimplifierGroup
    .def_property_readonly("_exon_connections",
        &SimplifierGroup::exon_connections,
        py::call_guard<py::gil_scoped_release>(),
        "Underlying exon connections being unsimplified")
    .def_property_readonly("introns_passed_src",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, src_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.introns_passed(), offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as source for introns")
    .def_property_readonly("introns_passed_dst",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, dst_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.introns_passed(), offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as target for introns")
    .def_property_readonly("junctions_passed_src",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, src_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.junctions_passed(), offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as source for junctions")
    .def_property_readonly("junctions_passed_dst",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, dst_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.junctions_passed(), offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments with evidence for unsimplification as target for junctions")
    .def_property_readonly("num_experiments",
        &SimplifierGroup::num_experiments,
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments in current simplifier group")
    .def("add_experiment", &SimplifierGroup::AddExperiment,
        py::call_guard<py::gil_scoped_release>(),
        "Increment evidence to unsimplify connections from SpliceGraphReads",
        py::arg("sg_reads"),
        py::arg("simplify_min_psi") = DEFAULT_BUILD_SIMPL_MINPSI,
        py::arg("simplify_minreads_annotated_junctions")
          = DEFAULT_BUILD_SIMPL_MINREADS_ANNOTATED_JUNCTION,
        py::arg("simplify_minreads_denovo_junctions")
          = DEFAULT_BUILD_SIMPL_MINREADS_DENOVO_JUNCTION,
        py::arg("simplify_minreads_introns")
          = DEFAULT_BUILD_SIMPL_MINREADS_INTRON)
    .def("update_connections", &SimplifierGroup::UpdateInplace,
        py::call_guard<py::gil_scoped_release>(),
        "Unsimplify introns/junctions with evidence, reset for next group",
        py::arg("simplifier_min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def(py::init<const std::shared_ptr<ExonConnections>&>(),
        py::call_guard<py::gil_scoped_release>(),
        "Start simplification group for the specified exon connections",
        py::arg("exon_connections"));
}

void init_PyEventsCoverage(pyEventsCoverage_t& pyEventsCoverage) {
  using majiq::Events;
  using majiq::EventsCoverage;
  using majiq::CoverageSummary;
  using majiq::CoverageBootstraps;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyEventsCoverage
    .def(py::init([](
            const std::shared_ptr<Events>& events,
            py::array_t<majiq::real_t> _numreads,
            py::array_t<majiq::real_t> _numbins,
            py::array_t<majiq::real_t> _bootstraps) {
          if (_numreads.ndim() != 1) {
            throw std::runtime_error("numreads must be 1D");
          } else if (_numbins.ndim() != 1) {
            throw std::runtime_error("numbins must be 1D");
          } else if (_bootstraps.ndim() != 2) {
            throw std::runtime_error("bootstraps must be 2D");
          }
          if (_numreads.shape(0) != _numbins.shape(0)
              || _numbins.shape(0) != _bootstraps.shape(0)) {
            throw std::runtime_error(
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
          return EventsCoverage{events,
              std::move(summaries_vec), std::move(coverage)};
          }),
        py::call_guard<py::gil_scoped_release>(),
        "construct events coverage from numpy arrays",
        py::arg("events"),
        py::arg("numreads"),
        py::arg("numbins"),
        py::arg("bootstraps"))
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
        py::call_guard<py::gil_scoped_release>(),
        "Obtain coverage for events from SJ junctions and introns",
        py::arg("events"),
        py::arg("sj_junctions"),
        py::arg("sj_introns"),
        py::arg("num_bootstraps") = DEFAULT_BUILD_NUM_BOOTSTRAPS,
        py::arg("pvalue_threshold") = DEFAULT_BUILD_STACK_PVALUE)
    .def_property_readonly("numreads",
        [](py::object& self_obj) -> py::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const size_t offset = offsetof(CoverageSummary, numreads);
        return ArrayFromVectorAndOffset<majiq::real_t, CoverageSummary>(
            self.summaries(), offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Readrate of each event connection after stacks removed")
    .def_property_readonly("numbins",
        [](py::object& self_obj) -> py::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const size_t offset = offsetof(CoverageSummary, numbins);
        return ArrayFromVectorAndOffset<majiq::real_t, CoverageSummary>(
            self.summaries(), offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of nonzero bins for each connection after stacks removed")
    .def_property_readonly("bootstraps",
        [](py::object& self_obj) -> py::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const CoverageBootstraps& x = self.bootstraps();
        py::array_t<majiq::real_t> result = py::array_t(
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
        py::call_guard<py::gil_scoped_release>(),
        "Bootstrapped read coverage or each connection after stacks removed")
    .def_property_readonly("_events", &EventsCoverage::events,
        py::call_guard<py::gil_scoped_release>(),
        "Events for which the coverage information is defined")
    .def("__len__", &EventsCoverage::num_connections,
        py::call_guard<py::gil_scoped_release>());
}

void init_pyEventsAlign(pyEventsAlign_t& pyEventsAlign) {
  using majiq::EventsAlign;
  using majiq::Events;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyEventsAlign
    .def(py::init<const Events&, const Events&>(),
        py::call_guard<py::gil_scoped_release>(),
        "Obtain indexes of matching events in the two input Events containers",
        py::arg("left_events"), py::arg("right_events"))
    .def_static(
        "events_match",
        [](const Events& left_events, const Events& right_events,
          size_t left_idx, size_t right_idx) {
        if (left_idx >= left_events.num_events()) {
        throw std::runtime_error("left_idx is out of range for left_events");
        } else if (right_idx >= right_events.num_events()) {
        throw std::runtime_error("right_idx is out of range for right_events");
        }
        return EventsAlign::EventsMatch(
            left_events, right_events, left_idx, right_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if selected event indexes share the exact same connections",
        py::arg("left_events"), py::arg("right_events"),
        py::arg("left_idx"), py::arg("right_idx"))
    .def_property_readonly("left_event_idx",
        [](py::object& self_obj) {
        EventsAlign& self = self_obj.cast<EventsAlign&>();
        const size_t offset = offsetof(EventsAlign::EventAligned, left_idx_);
        return ArrayFromVectorAndOffset<size_t, EventsAlign::EventAligned>(
            self.matched_, offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indexes for events in left_events used in constructor")
    .def_property_readonly("right_event_idx",
        [](py::object& self_obj) {
        EventsAlign& self = self_obj.cast<EventsAlign&>();
        const size_t offset = offsetof(EventsAlign::EventAligned, right_idx_);
        return ArrayFromVectorAndOffset<size_t, EventsAlign::EventAligned>(
            self.matched_, offset, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indexes for events in right_events used in constructor");
  return;
}

void init_PyEvents(pyEvents_t& pyEvents) {
  using majiq::Event;
  using majiq::Events;
  using majiq::ConnectionIndex;
  using majiq_pybind::ArrayFromOffsetsVector;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyEvents
    .def(py::init([](
            const std::shared_ptr<majiq::GeneIntrons>& introns,
            const std::shared_ptr<majiq::GeneJunctions>& junctions,
            // make events vector (ref_exon_idx, event_type)
            py::array_t<size_t> _ref_exon_idx,
            py::array_t<std::array<char, 1>> _event_type,
            // connection offsets (1 longer than events)
            py::array_t<size_t> _offsets,
            // connections (is_intron, connection_idx)
            py::array_t<bool> _is_intron,
            py::array_t<size_t> _connection_idx) {
          using majiq::EventType;
          using majiq::ConnectionIndex;
          auto check_1d = [](const auto& x) {
            if (x.ndim() != 1)
            throw std::runtime_error("Events arrays must be 1D"); };
          check_1d(_ref_exon_idx);
          check_1d(_event_type);
          check_1d(_offsets);
          check_1d(_is_intron);
          check_1d(_connection_idx);
          if (_ref_exon_idx.shape(0) != _event_type.shape(0)) {
            throw std::runtime_error(
                "ref_exon_idx and event_type must have same length");
          }
          std::vector<Event> event_vec(_ref_exon_idx.shape(0));
          {
            auto ref_exon_idx = _ref_exon_idx.unchecked<1>();
            auto event_type = _event_type.unchecked<1>();
            for (size_t i = 0; i < event_vec.size(); ++i) {
              event_vec[i] = Event{
                ref_exon_idx(i), static_cast<EventType>(event_type(i)[0])};
            }
          }
          std::vector<size_t> offsets_vec(_offsets.shape(0));
          {
            auto offsets = _offsets.unchecked<1>();
            for (size_t i = 0; i < offsets_vec.size(); ++i) {
              offsets_vec[i] = offsets(i);
            }
          }
          if (_is_intron.shape(0) != _connection_idx.shape(0)) {
            throw std::runtime_error(
                "is_intron and connection_idx must have same length");
          }
          std::vector<ConnectionIndex> connections_vec(_is_intron.shape(0));
          {
            auto is_intron = _is_intron.unchecked<1>();
            auto connection_idx = _connection_idx.unchecked<1>();
            for (size_t i = 0; i < connections_vec.size(); ++i) {
              if (connection_idx(i)
                  >= (is_intron(i) ? introns->size() : junctions->size())) {
                throw std::invalid_argument("connection_idx has out of range values");
              }
              connections_vec[i] = ConnectionIndex{
                is_intron(i), connection_idx(i)};
            }
          }
          return Events{introns, junctions, std::move(event_vec),
            std::move(offsets_vec), std::move(connections_vec)};
          }),
        py::call_guard<py::gil_scoped_release>(),
        "Initialize events object from numpy arrays",
        py::arg("introns"),
        py::arg("junctions"),
        py::arg("ref_exon_idx"),
        py::arg("event_type"),
        py::arg("offsets"),
        py::arg("is_intron"),
        py::arg("connection_idx"))
    .def_property_readonly("introns", &Events::introns,
        py::call_guard<py::gil_scoped_release>(),
        "underlying introns")
    .def_property_readonly("junctions", &Events::junctions,
        py::call_guard<py::gil_scoped_release>(),
        "underlying junctions")
    .def_property_readonly("ref_exon_idx",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(Event, ref_exon_idx_);
        return ArrayFromVectorAndOffset<size_t, Event>(
            self.events(), offset, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "Event reference exon")
    .def_property_readonly("event_type",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(Event, type_);
        return ArrayFromVectorAndOffset<std::array<char, 1>, Event>(
            self.events(), offset, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "Event type")
    .def_property_readonly("_offsets",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.connection_offsets(), 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Raw offsets for events into connections")
    .def_property_readonly("connection_idx_start",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.connection_offsets(), true, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "First index into event connections for each event")
    .def_property_readonly("connection_idx_end",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.connection_offsets(), false, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "One after last index into event connections for each event")
    .def_property_readonly("is_intron",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(ConnectionIndex, is_intron_);
        return ArrayFromVectorAndOffset<bool, ConnectionIndex>(
            self.connections(), offset, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "Event connection is intron (false --> is junction)")
    .def_property_readonly("idx",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(ConnectionIndex, idx_);
        return ArrayFromVectorAndOffset<size_t, ConnectionIndex>(
            self.connections(), offset, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "Event connection index into corresponding Introns or Junctions")
    .def_property_readonly("connection_event_idx",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.connection_event_idx(), 0, self_obj); },
        py::call_guard<py::gil_scoped_release>(),
        "Event connection index back into events")
    .def_property_readonly("num_events", &Events::num_events,
        py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("num_connections", &Events::num_connections,
        py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("num_junctions", &Events::num_junctions,
        py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("num_introns", &Events::num_introns,
        py::call_guard<py::gil_scoped_release>())
    .def("connection_gene_idx",
        [](const Events& self, py::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        return self.connection_gene(idx).idx_; };
        return py::vectorize(f)(connection_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "gene_idx for specified connection indexes",
        py::arg("connection_idx"))
    .def("connection_start",
        [](const Events& self, py::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        return self.connection_start(idx); };
        return py::vectorize(f)(connection_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "start for specified connection indexes",
        py::arg("connection_idx"))
    .def("connection_end",
        [](const Events& self, py::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        return self.connection_end(idx); };
        return py::vectorize(f)(connection_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "end for specified connection indexes",
        py::arg("connection_idx"))
    .def("connection_denovo",
        [](const Events& self, py::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        return self.connection_denovo(idx); };
        return py::vectorize(f)(connection_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "denovo status for specified connection indexes",
        py::arg("connection_idx"))
    .def("connection_other_exon_idx",
        [](const Events& self, py::array_t<size_t> connection_idx) {
        auto f = [&self](size_t idx) {
        return self.connection_other_exon_idx(idx); };
        return py::vectorize(f)(connection_idx);
        },
        py::call_guard<py::gil_scoped_release>(),
        "index for other exon for specified connection indexes",
        py::arg("connection_idx"))
    .def("__len__", &Events::size, py::call_guard<py::gil_scoped_release>());
}

void init_pyExonConnections(pyExonConnections_t& pyExonConnections) {
  using majiq::Event;
  using majiq::EventType;
  using majiq::ExonConnections;
  using ExonsPtrT = std::shared_ptr<majiq::Exons>;
  using IntronsPtrT = std::shared_ptr<majiq::GeneIntrons>;
  using JunctionsPtrT = std::shared_ptr<majiq::GeneJunctions>;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyExonConnections
    .def(
        py::init<const ExonsPtrT&, const IntronsPtrT&, const JunctionsPtrT&>(),
        py::call_guard<py::gil_scoped_release>(),
        "Track junctions/introns by exons/events assocciated with them",
        py::arg("exons"),
        py::arg("introns"),
        py::arg("junctions"))
    .def_property_readonly("src_intron_idx",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.idx_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "intron_idx for exon_connections in src_exon sorted order")
    .def_property_readonly("src_intron_exon_offsets",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.exon_offsets_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "offsets into src_intron_idx for each exon")
    .def_property_readonly("dst_intron_idx",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.idx_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "intron_idx for exon_connections in dst_exon sorted order")
    .def_property_readonly("dst_intron_exon_offsets",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_introns();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.exon_offsets_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "offsets into dst_intron_idx for each exon")
    .def_property_readonly("src_junction_idx",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.idx_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "junction_idx for exon_connections in src_exon sorted order")
    .def_property_readonly("src_junction_exon_offsets",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.src_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.exon_offsets_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "offsets into src_junction_idx for each exon")
    .def_property_readonly("dst_junction_idx",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.idx_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "junction_idx for exon_connections in dst_exon sorted order")
    .def_property_readonly("dst_junction_exon_offsets",
        [](py::object& self_obj) {
        ExonConnections& self = self_obj.cast<ExonConnections&>();
        const auto& indexes = self.dst_junctions();
        return ArrayFromVectorAndOffset<size_t, size_t>(indexes.exon_offsets_, 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "offsets into dst_junction_idx for each exon")
    .def_property_readonly("_exons", &ExonConnections::exons,
        py::call_guard<py::gil_scoped_release>(), "underlying exons")
    .def_property_readonly("_introns", &ExonConnections::introns,
        py::call_guard<py::gil_scoped_release>(), "underlying introns")
    .def_property_readonly("_junctions", &ExonConnections::junctions,
        py::call_guard<py::gil_scoped_release>(), "underlying junctions")
    .def("strict_lsvs", &ExonConnections::StrictLSVs,
        py::call_guard<py::gil_scoped_release>(), "Construct strict LSV Events")
    .def("permissive_lsvs", &ExonConnections::PermissiveLSVs,
        py::call_guard<py::gil_scoped_release>(), "Construct permissive LSV Events")
    .def("source_lsvs", &ExonConnections::SourceLSVs,
        py::call_guard<py::gil_scoped_release>(), "Construct source LSV Events")
    .def("target_lsvs", &ExonConnections::TargetLSVs,
        py::call_guard<py::gil_scoped_release>(), "Construct target LSV Events")
    .def("constitutive", &ExonConnections::ConstitutiveEvents,
        py::call_guard<py::gil_scoped_release>(), "Construct Constitutive Events")
    .def("events_for",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) {
        if (exon_idx.ndim() != 1) {
          throw std::invalid_argument("exon_idx must be one-dimensional");
        } else if (is_source.ndim() != 1) {
          throw std::invalid_argument("is_source must be one-dimensional");
        } else if (is_source.shape(0) != exon_idx.shape(0)) {
          throw std::invalid_argument(
              "exon_idx and is_source must have same size");
        }
        auto _exon_idx = exon_idx.unchecked<1>();
        auto _is_source = is_source.unchecked<1>();
        std::vector<Event> events(_exon_idx.shape(0));
        for (py::ssize_t i = 0; i < _exon_idx.shape(0); ++i) {
          if (_exon_idx(i) >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has out of range values");
          }
          events[i] = Event{_exon_idx(i), _is_source(i) ? EventType::SRC_EVENT : EventType::DST_EVENT};
        }
        return self.CreateEvents(std::move(events));
        },
        py::call_guard<py::gil_scoped_release>(),
        "Construct events for specified exons/directions",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("has_intron",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.has_intron(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if events have introns or not",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("event_size",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<size_t> {
        auto f = [&self](size_t idx, bool is_src) -> size_t {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.event_size(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate size of event",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("passed",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.passed(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event was passed",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("redundant",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.redundant(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event was redundant",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("is_target_LSV",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_target) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_target_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_target);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event is target LSV",
        py::arg("exon_idx"), py::arg("is_target"))
    .def("is_source_LSV",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_source_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event is source LSV",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("is_permissive_LSV",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_permissive_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event is permissive LSV",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("is_strict_LSV",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_strict_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event is strict LSV",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("is_constitutive",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          if (idx >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          return self.is_constitutive(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Indicate if event is constitutive",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("event_id",
        [](const ExonConnections& self,
          py::array_t<size_t> _exon_idx,
          py::array_t<std::array<char, 1>> _event_type) {
        if (_exon_idx.ndim() != 1 || _event_type.ndim() != 1) {
        throw std::runtime_error("exon_idx and event_type must be 1D");
        } else if (_exon_idx.shape(0) != _event_type.shape(0)) {
        throw std::runtime_error("exon_idx and event_type must have same shape");
        }
        std::vector<std::string> result(_exon_idx.shape(0));
        {
          auto exon_idx = _exon_idx.unchecked<1>();
          auto event_type = _event_type.unchecked<1>();
          for (size_t i = 0; i < result.size(); ++i) {
          if (exon_idx(i) >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          result[i] = self.id(Event{
              exon_idx(i), static_cast<EventType>(event_type(i)[0])});
          }
        }
        return result;
        },
        py::call_guard<py::gil_scoped_release>(),
        "List of event_id for specified events",
        py::arg("exon_idx"), py::arg("event_type"))
    .def("event_description",
        [](const ExonConnections& self,
          py::array_t<size_t> _exon_idx,
          py::array_t<std::array<char, 1>> _event_type) {
        if (_exon_idx.ndim() != 1 || _event_type.ndim() != 1) {
        throw std::runtime_error("exon_idx and event_type must be 1D");
        } else if (_exon_idx.shape(0) != _event_type.shape(0)) {
        throw std::runtime_error("exon_idx and event_type must have same shape");
        }
        std::vector<std::string> result(_exon_idx.shape(0));
        {
          auto exon_idx = _exon_idx.unchecked<1>();
          auto event_type = _event_type.unchecked<1>();
          for (size_t i = 0; i < result.size(); ++i) {
          if (exon_idx(i) >= self.num_exons()) {
            throw std::invalid_argument("exon_idx has values out of range");
          }
          result[i] = self.description(Event{
              exon_idx(i), static_cast<EventType>(event_type(i)[0])});
          }
        }
        return result;
        },
        py::call_guard<py::gil_scoped_release>(),
        "List of description for specified events",
        py::arg("exon_idx"), py::arg("event_type"));
}

void init_SJIntronsBins(pySJIntronsBins_t& pySJIntronsBins) {
  using majiq::SJIntronsBins;
  define_sjbins_properties(pySJIntronsBins);
  pySJIntronsBins
    .def_static("from_bam", &SJIntronsBins::FromBam,
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Load introns and per-bin counts for an aligned BAM file

        Parameters
        ----------
        bam_path: str
            Path for input BAM fille
        num_bins: int
            Number of bins to split coverage. Typically set to num_positions
            from junctions
        exons: Exons
            Gene exons defining potential introns for coverage
        gene_introns: GeneIntrons
            Gene introns indicating annotated introns for coverage
        experiment_strandness: ExperimentStrandness
            Strandness of RNA-seq library
        nthreads: int
            Number of threads to use when reading in BAM file
        )pbdoc",
        py::arg("bam_path"),
        py::arg("num_bins"),
        py::arg("exons"),
        py::arg("gene_introns"),
        py::arg("experiment_strandness") = DEFAULT_BAM_STRANDNESS,
        py::arg("nthreads") = DEFAULT_BAM_NTHREADS);
}

void init_GeneIntrons(pyGeneIntrons_t& pyGeneIntrons) {
  using majiq::GeneIntrons;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties(pyGeneIntrons);
  pyGeneIntrons
    .def(py::init([](
            std::shared_ptr<majiq::Genes> genes,
            py::array_t<size_t> gene_idx,
            py::array_t<position_t> start,
            py::array_t<position_t> end,
            py::array_t<bool> denovo,
            py::array_t<bool> passed_build,
            py::array_t<bool> simplified) {
          return MakeConnections<GeneIntrons>(genes, ConnectionsArrays{
              gene_idx, start, end, denovo, passed_build, simplified});
        }),
        py::call_guard<py::gil_scoped_release>(),
        "Create GeneIntrons using Genes and arrays defining each intron",
        py::arg("genes"),
        py::arg("gene_idx"), py::arg("start"), py::arg("end"),
        py::arg("denovo"), py::arg("passed_build"), py::arg("simplified"))
    .def("update_flags_from", &GeneIntrons::UpdateFlagsFrom,
        py::call_guard<py::gil_scoped_release>(),
        "Update intron flags using introns that overlap from input",
        py::arg("donor_introns"))
    .def("build_group", [](std::shared_ptr<GeneIntrons>& gene_introns) {
        return majiq::GroupIntronsGenerator(gene_introns);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Create build group to update passed introns in place")
    .def("filter_passed", &GeneIntrons::FilterPassed,
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Get subset of introns that passed build filters

        Parameters
        ----------
        keep_annotated: bool
            Keep all annotated introns regardless of whether they passed
        discard_denovo: bool
            Discard all denovo introns regardless of whether they passed
        )pbdoc",
        py::arg("keep_annotated") = DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        py::arg("discard_denovo") = !DEFAULT_BUILD_DENOVO_IR)
    .def("__repr__", [](const GeneIntrons& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneIntrons<" << self.size() << " total>";
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>());
}

void init_SJJunctions(pySJJunctions_t& pySJJunctions) {
  using majiq::position_t;
  using majiq::SJJunctions;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties(pySJJunctions);
  pySJJunctions
    .def(py::init([](
            std::shared_ptr<majiq::Contigs> contigs,
            py::array_t<size_t> _contig_idx,
            py::array_t<position_t> _start,
            py::array_t<position_t> _end,
            py::array_t<std::array<char, 1>> _strand) {
          auto contig_idx = _contig_idx.unchecked<1>();
          auto start = _start.unchecked<1>();
          auto end = _end.unchecked<1>();
          auto strand = _strand.unchecked<1>();
          // fill in junctions
          std::vector<majiq::SJJunction> sj_vec(start.shape(0));
          for (size_t i = 0; i < sj_vec.size(); ++i) {
            if (contig_idx(i) >= contigs->size()) {
              throw std::invalid_argument("contig_idx has values out of range");
            }
            sj_vec[i] = majiq::SJJunction{
              majiq::KnownContig{contig_idx(i), contigs},
              majiq::OpenInterval{start(i), end(i)},
              static_cast<majiq::GeneStrandness>(strand(i)[0])};
          }
          return std::make_shared<majiq::SJJunctions>(
              contigs, std::move(sj_vec));
        }),
        py::call_guard<py::gil_scoped_release>(),
        "Create SJJunctions object from contigs and arrays",
        py::arg("contigs"),
        py::arg("contig_idx"), py::arg("start"), py::arg("end"),
        py::arg("strand"))
    .def("to_unstranded", &SJJunctions::ToUnstranded,
        py::call_guard<py::gil_scoped_release>(),
        "Create unstranded SJJunctions from self")
    .def("flip_strand", &SJJunctions::FlipStrand,
        py::call_guard<py::gil_scoped_release>(),
        "Create SJJunctions with strands flipped in sorted order")
    .def_property_readonly("_contigs", &SJJunctions::parents,
        py::call_guard<py::gil_scoped_release>(),
        "Underlying contigs corresponding to contig_idx");
}
void init_SJJunctionsBins(pySJJunctionsBins_t& pySJJunctionsBins) {
  using majiq::SJJunctionsBins;
  define_sjbins_properties(pySJJunctionsBins);
  pySJJunctionsBins
    .def("project_reads", &SJJunctionsBins::ProjectReads,
        py::call_guard<py::gil_scoped_release>(),
        "Project reads to a different set of SJJunctions")
    .def_static("from_bam", &SJJunctionsBins::FromBam,
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Load junctions and per-position counts for an aligned BAM file

        Parameters
        ----------
        bam_path: str
            Path for input BAM file
        experiment_strandness: ExperimentStrandness
            Strandness of RNA-seq library
        nthreads: int
            Number of threads to use when reading in BAM file
        )pbdoc",
        py::arg("bam_path"),
        py::arg("experiment_strandness") = DEFAULT_BAM_STRANDNESS,
        py::arg("nthreads") = DEFAULT_BAM_NTHREADS);
}

void init_pyGroupJunctionsGen(pyGroupJunctionsGen_t& pyGroupJunctionsGen) {
  pyGroupJunctionsGen
    .def(
      py::init<const std::shared_ptr<majiq::GeneJunctions>&,
      const std::shared_ptr<majiq::Exons>&>(),
      py::call_guard<py::gil_scoped_release>(),
      "Accumulate SJJunctions for build group for input junctions/exons",
      py::arg("junctions"),
      py::arg("exons"))
    .def("add_experiment", &majiq::GroupJunctionsGenerator::AddExperiment,
        py::call_guard<py::gil_scoped_release>(),
        "Increment count of passed junctions from input experiment",
        py::arg("sjp"),
        py::arg("thresholds") = DEFAULT_THRESHOLDS,
        py::arg("add_denovo") = DEFAULT_BUILD_DENOVO_JUNCTIONS)
    .def("pass_known_inplace",
        &majiq::GroupJunctionsGenerator::UpdateKnownInplace,
        py::call_guard<py::gil_scoped_release>(),
        "Update known junctions with their build status (ignore denovos)",
        py::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def_property_readonly("num_experiments",
        &majiq::GroupJunctionsGenerator::num_experiments,
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments that have been added to this group")
    .def_property_readonly("num_known",
        &majiq::GroupJunctionsGenerator::num_annotated,
        py::call_guard<py::gil_scoped_release>(),
        "Number of junctions known when constructed")
    .def_property_readonly("num_denovo",
        &majiq::GroupJunctionsGenerator::num_denovo,
        py::call_guard<py::gil_scoped_release>(),
        "Number of denovo junctions passing experiment-filters for 1+ inputs")
    .def("__len__", &majiq::GroupJunctionsGenerator::size,
        py::call_guard<py::gil_scoped_release>());
}

void init_pyGroupIntronsGen(pyGroupIntronsGen_t& pyGroupIntronsGen) {
  using majiq::GroupIntronsGenerator;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyGroupIntronsGen
    .def(
        py::init<const std::shared_ptr<majiq::GeneIntrons>&>(),
        py::call_guard<py::gil_scoped_release>(),
        "Accumulate SJIntronsBins for build group of introns",
        py::arg("gene_introns"))
    .def_property_readonly("num_experiments",
        &GroupIntronsGenerator::num_experiments,
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments in current group")
    .def_property_readonly("_introns", &GroupIntronsGenerator::introns,
        py::call_guard<py::gil_scoped_release>())
    .def_property_readonly("num_passed",
        [](py::object& self_obj) -> py::array_t<size_t> {
        GroupIntronsGenerator& self = self_obj.cast<GroupIntronsGenerator&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.num_passed(), 0, self_obj);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Number of experiments or which each intron has passed")
    .def("__len__", &GroupIntronsGenerator::size,
        py::call_guard<py::gil_scoped_release>())
    .def("add_experiment", &GroupIntronsGenerator::AddExperiment,
        py::call_guard<py::gil_scoped_release>(),
        "Add SJIntronsBins to build group",
        py::arg("sj"),
        py::arg("thresholds") = DEFAULT_THRESHOLDS)
    .def("update_introns", &GroupIntronsGenerator::UpdateInplace,
        py::call_guard<py::gil_scoped_release>(),
        "Pass introns that pass min-experiments in place, reset for next group",
        py::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS);
}

void init_pyPassedJunctionsGen(pyPassedJunctionsGen_t& pyPassedJunctionsGen) {
  pyPassedJunctionsGen
    .def(
      py::init<const std::shared_ptr<majiq::GeneJunctions>&>(),
      py::call_guard<py::gil_scoped_release>(),
      R"pbdoc(
      Accumulator of GroupJunctionsGenerator for different build groups
      )pbdoc",
      py::arg("junctions"))
    .def("add_group", &majiq::PassedJunctionsGenerator::AddGroup,
        py::call_guard<py::gil_scoped_release>(),
        "Combine passed junctions from build group of experiments",
        py::arg("group"),
        py::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def("get_passed", &majiq::PassedJunctionsGenerator::PassedJunctions,
        py::call_guard<py::gil_scoped_release>(),
        "Get static, array-based representation of passed junctions",
        py::arg("denovo_simplified"))
    .def_property_readonly("num_known",
        &majiq::PassedJunctionsGenerator::num_annotated,
        py::call_guard<py::gil_scoped_release>(),
        "Number of junctions known when constructed")
    .def_property_readonly("num_denovo",
        &majiq::PassedJunctionsGenerator::num_denovo,
        py::call_guard<py::gil_scoped_release>(),
        "Number of denovo junctions passing filters, previously not known")
    .def("__len__", &majiq::PassedJunctionsGenerator::size,
        py::call_guard<py::gil_scoped_release>());
}

void init_SpliceGraph(py::class_<majiq::SpliceGraph>& pySpliceGraph) {
  using majiq::SpliceGraph;
  using majiq::position_t;
  pySpliceGraph
    // expose constructor from individual components
    .def(py::init<const std::shared_ptr<majiq::Contigs>&,
                  const std::shared_ptr<majiq::Genes>&,
                  const std::shared_ptr<majiq::Exons>&,
                  const std::shared_ptr<majiq::GeneJunctions>&,
                  const std::shared_ptr<majiq::GeneIntrons>&>(),
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Initialize splicegraph from components

        Initialize splicegraph from components. Typically will want to use the
        factory methods `from_gff3` to create all components
        )pbdoc",
        py::arg("contigs"), py::arg("genes"), py::arg("exons"),
        py::arg("junctions"), py::arg("introns"))
    // constructors from gff3
    .def_static("from_gff3",
        [](std::string gff3_path, bool process_ir,
            majiq::gff3::featuretype_map_t gff3_types) {
          using majiq::gff3::GFF3ExonHierarchy;
          using majiq::gff3::GFF3TranscriptModels;
          using majiq::gff3::ToTranscriptModels;
          // load gff3 exon hierarchy, convert to MAJIQ gene/transcript/exons
          auto gff3_models
            = ToTranscriptModels(GFF3ExonHierarchy{gff3_path, gff3_types});
          // TODO(jaicher) print info about skipped types more Python-friendly
          for (const auto& [tx_type, tx_ct]
              : gff3_models.skipped_transcript_type_ct_) {
            std::cerr << "Skipped exons for "
              << tx_ct
              << " potential transcripts with unaccepted parent type '"
              << tx_type << "'\n";
          }
          for (const auto& [gene_type, gene_ct]
              : gff3_models.skipped_gene_type_ct_) {
            std::cerr << "Skipped exons for "
              << gene_ct
              << " potential genes with unaccepted top-level type '"
              << gene_type << "'\n";
          }
          // convert to splicegraph
          return gff3_models.models_.ToSpliceGraph(process_ir);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Create splicegraph from input GFF3 file",
        py::arg("gff3_path"),
        py::arg("process_ir") = DEFAULT_BUILD_PROCESS_IR,
        py::arg_v("gff3_types", majiq::gff3::default_gff3_types,
            "new_majiq._default_gff3_types()"))
    // XXX debug
    .def_static("infer_exons", &SpliceGraph::InferExons,
        py::call_guard<py::gil_scoped_release>(),
        "Infer exons from base annotated exons and junctions",
        py::arg("base_exons"), py::arg("junctions"))
    .def("make_group_junctions", &SpliceGraph::MakeGroupGenerator,
        py::call_guard<py::gil_scoped_release>(),
        "Create GroupJunctionsGenerator for the splicegraph junctions/exons")
    .def("make_build_junctions", &SpliceGraph::MakePassedGenerator,
        py::call_guard<py::gil_scoped_release>(),
        "Create PassedJunctionsGenerator for splicegraph junctions")
    .def("close_to_annotated_exon",
        [](SpliceGraph& sg, size_t gene_idx, position_t x, bool to_following) {
        if (gene_idx >= sg.genes()->size()) {
          throw std::invalid_argument("gene_idx is out of range");
        }
        majiq::KnownGene g = (*sg.genes())[gene_idx];
        const majiq::Exons& exons = *sg.exons();
        return to_following
          ? majiq::detail::CloseToFollowingAnnotatedExon(exons, g, x)
          : majiq::detail::CloseToPrecedingAnnotatedExon(exons, g, x);
        },
        py::call_guard<py::gil_scoped_release>(),
        "True if position close to following/preceding annotated exon in gene",
        py::arg("gene_idx"),
        py::arg("x"),
        py::arg("to_following") = true)
    // access underlying data
    .def_property_readonly("_exons", &SpliceGraph::exons,
        py::call_guard<py::gil_scoped_release>(),
        "Access the splicegraph's exons")
    .def_property_readonly("_introns", &SpliceGraph::introns,
        py::call_guard<py::gil_scoped_release>(),
        "Access the splicegraph's introns")
    .def_property_readonly("_junctions", &SpliceGraph::junctions,
        py::call_guard<py::gil_scoped_release>(),
        "Access the splicegraph's junctions")
    .def_property_readonly("_genes", &SpliceGraph::genes,
        py::call_guard<py::gil_scoped_release>(),
        "Access the splicegraph's genes")
    .def_property_readonly("_contigs", &SpliceGraph::contigs,
        py::call_guard<py::gil_scoped_release>(),
        "Access the splicegraph's contigs")
    .def_property_readonly("_exon_connections", &SpliceGraph::exon_connections,
        py::call_guard<py::gil_scoped_release>(),
        "Access the splicegraph's exon connections")
    // get sj introns on which coverage may be read
    .def("sj_introns", [](SpliceGraph& sg, bool stranded) {
        return majiq::SJIntrons::FromGeneExonsAndIntrons(
            *sg.exons(), *sg.introns(), stranded);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Get contig introns (by strand or not) for splicegraph",
        py::arg("stranded"))
    // that's really for debugging because it depends on experiment. Instead,
    // just create them while loading BAM
    .def("sj_introns_from_bam",
        [](SpliceGraph& sg, const char* infile, majiq::junction_pos_t num_bins,
          majiq::ExperimentStrandness exp_strandness, int nthreads) {
        return majiq::SJIntronsBins::FromBam(infile, num_bins, *sg.exons(),
            *sg.introns(), exp_strandness, nthreads);
        },
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Load introns and per-bin counts for an aligned BAM file

        Parameters
        ----------
        bam_path: str
            Path for input BAM fille
        num_bins: int
            Number of bins to split coverage. Typically set to num_positions
            from junctions
        experiment_strandness: ExperimentStrandness
            Strandness of RNA-seq library
        nthreads: int
            Number of threads to use when reading in BAM file
        )pbdoc",
        py::arg("bam_path"),
        py::arg("num_bins"),
        py::arg("experiment_strandness") = DEFAULT_BAM_STRANDNESS,
        py::arg("nthreads") = DEFAULT_BAM_NTHREADS)
    // string representation of splicegraph
    .def("__repr__", [](const SpliceGraph& sg) -> std::string {
        std::ostringstream oss;
        oss << sg;
        return oss.str();
        },
        py::call_guard<py::gil_scoped_release>());
}

void init_pyExperimentThresholds(
    pyExperimentThresholds_t& pyExperimentThresholds) {
  using majiq::ExperimentThresholds;
  using majiq::junction_ct_t;
  using majiq::junction_pos_t;
  using majiq::real_t;
  pyExperimentThresholds
    .def(py::init<junction_ct_t, junction_ct_t, junction_pos_t, real_t, real_t, real_t>(),
        py::call_guard<py::gil_scoped_release>(),
        "Thresholds on intron and junction coverage for inclusion in SpliceGraph",
        py::arg("minreads") = DEFAULT_BUILD_MINREADS,
        py::arg("mindenovo") = DEFAULT_BUILD_MINDENOVO,
        py::arg("minpos") = DEFAULT_BUILD_MINPOS,
        py::arg("max_pctbins") = DEFAULT_BUILD_MAX_PCTBINS,
        py::arg("junction_acceptance_probability")
          = DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
        py::arg("intron_acceptance_probability")
          = DEFAULT_BUILD_MATCH_INTRON_PROBABILITY)
    .def_property_readonly("minreads",
        [](const ExperimentThresholds& x) { return x.minreads_; },
        py::call_guard<py::gil_scoped_release>(),
        "Minimum number of reads for an annotated junction to pass")
    .def_property_readonly("mindenovo",
        [](const ExperimentThresholds& x) { return x.mindenovo_; },
        py::call_guard<py::gil_scoped_release>(),
        "Minimum number of reads for a denovo junction to pass")
    .def_property_readonly("minpos",
        [](const ExperimentThresholds& x) { return x.minpos_; },
        py::call_guard<py::gil_scoped_release>(),
        "Minimum number of nonzero positions for a junction to pass")
    .def_property_readonly("max_pctbins",
        [](const ExperimentThresholds& x) { return x.max_pctbins_; },
        py::call_guard<py::gil_scoped_release>(),
        "Maximum percentage of bins to require coverage in for intron to pass")
    .def_property_readonly("junction_acceptance_probability",
        [](const ExperimentThresholds& x) {
        return x.junction_acceptance_probability_; },
        py::call_guard<py::gil_scoped_release>(),
        "Set intron thresholds to match junction readrate with probability")
    .def_property_readonly("intron_acceptance_probability",
        [](const ExperimentThresholds& x) {
        return x.intron_acceptance_probability_; },
        py::call_guard<py::gil_scoped_release>(),
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
        py::call_guard<py::gil_scoped_release>())
    .def("intron_thresholds_generator",
        &ExperimentThresholds::intron_thresholds_generator,
        py::call_guard<py::gil_scoped_release>(),
        R"pbdoc(
        Create IntronThresholdsGenerator for intron thresholds by length

        Parameters
        ----------
        total_bins: int
            Maximum number of bins/positions for each junction/intron for
            quantification
        )pbdoc",
        py::arg("total_bins"));
}

void init_pyIntronThresholdsGenerator(
    pyIntronThresholdsGenerator_t& pyIntronThresholdsGenerator) {
  using majiq::IntronThresholdsGenerator;
  using majiq::IntronThresholds;
  using majiq::junction_ct_t;
  using majiq::junction_pos_t;
  PYBIND11_NUMPY_DTYPE(IntronThresholds, minreads_, minbins_, mincov_);
  pyIntronThresholdsGenerator
    .def("__call__",
        [](const IntronThresholdsGenerator& gen,
          const py::array_t<junction_pos_t>& intron_lengths) {
        // function per-element of intron_lengths
        auto f = [&gen](junction_pos_t x) { return gen(x); };
        return py::vectorize(f)(intron_lengths);
        },
        py::call_guard<py::gil_scoped_release>(),
        "Get intron thresholds for specified intron lengths",
        py::arg("intron_lengths"));
}

void init_SpliceGraphAll(py::module_& m) {
  m.def("rng_seed",
      [](int64_t x) { global_rng_pool.seed(x); },
      "Set seed for pool of random number generators in new_majiq.internals",
      py::arg("seed"));
  m.def("rng_resize",
      [](int64_t n) { global_rng_pool.resize(n); },
      "Resize pool of random number generators for at least n simultaneous threads",
      py::arg("n"));

  using majiq::Contigs;
  using majiq::Genes;
  using majiq::Exons;
  using majiq::GeneIntrons;
  using majiq::GeneJunctions;
  using majiq::SpliceGraph;
  using majiq::SJJunctions;
  using majiq::SJJunctionsBins;
  using majiq::ExperimentStrandness;
  using majiq::GeneStrandness;
  using majiq::SJIntrons;
  auto pyContigs = pyContigs_t(m, "Contigs", "Splicegraph contigs");
  auto pyGenes = pyGenes_t(m, "Genes",
      "Splicegraph genes");
  auto pyExons = pyExons_t(m, "Exons",
      "Splicegraph exons");
  auto pyGeneIntrons = pyGeneIntrons_t(m, "GeneIntrons",
      "Splicegraph introns");
  auto pyGeneJunctions = pyGeneJunctions_t(
      m, "GeneJunctions", "Splicegraph junctions");
  auto pyGeneJunctionsAccumulator = pyGeneJunctionsAccumulator_t(
      m, "GeneJunctionsAccumulator", "Accumulator of splicegraph junctions");
  auto pyEvents = pyEvents_t(
      m, "Events", "Events from reference exon with junctions/introns");
  auto pyEventsAlign = pyEventsAlign_t(
      m, "EventsAlign", "Indexes to shared events between two Events containers");
  auto pyEventsCoverage = pyEventsCoverage_t(
      m, "EventsCoverage", "Coverage over events for a single experiment");
  auto pyExonConnections = pyExonConnections_t(
      m, "ExonConnections", "Connections from exons to junctions, introns");
  auto pySJIntrons = pySJIntrons_t(m, "SJIntrons");
  auto pySJIntronsBins = pySJIntronsBins_t(m, "SJIntronsBins",
      "Summarized and per-bin counts for introns from an experiment");
  auto pySpliceGraphReads = pySpliceGraphReads_t(
      m, "SpliceGraphReads", "Raw readrates for each intron and junction");
  auto pySJJunctions = pySJJunctions_t(
      m, "SJJunctions", "Summarized junction counts for an experiment");
  auto pySJJunctionsBins = pySJJunctionsBins_t(
      m, "SJJunctionsBins",
      "Summarized and per-position counts for an experiment");
  auto pyGroupJunctionsGen = pyGroupJunctionsGen_t(
      m, "GroupJunctionsGenerator",
      "Accumulator of SJJunctions in the same build group");
  auto pyPassedJunctionsGen = pyPassedJunctionsGen_t(
      m, "PassedJunctionsGenerator",
      "Accumulator of GroupJunctionsGenerator from multiple build groups");
  auto pyGroupIntronsGen = pyGroupIntronsGen_t(
      m, "GroupIntronsGenerator",
      "Accumulator of SJIntronsBins in the same build group");
  auto pySimplifierGroup = pySimplifierGroup_t(
      m, "SimplifierGroup",
      "Accumulator of SpliceGraphReads to unsimplify introns and junctions");
  auto pyGeneStrandness = py::enum_<GeneStrandness>(
      m, "GeneStrandness")
    .value("forward", GeneStrandness::FORWARD)
    .value("reverse", GeneStrandness::REVERSE)
    .value("ambiguous", GeneStrandness::AMBIGUOUS);
  auto pyExperimentStrandness = py::enum_<ExperimentStrandness>(
      m, "ExperimentStrandness")
    .value("FORWARD", ExperimentStrandness::FORWARD)
    .value("REVERSE", ExperimentStrandness::REVERSE)
    .value("NONE", ExperimentStrandness::NONE);
  auto pyGFFFeatureTypes = py::enum_<majiq::gff3::FeatureType>(
      m, "GFF3FeatureType")
    .value("EXON", majiq::gff3::FeatureType::EXON,
        "Indicates feature type defining transcript exons")
    .value("ACCEPT_GENE", majiq::gff3::FeatureType::ACCEPT_GENE,
        R"pbdoc(
        Indicates feature type that would be accepted as a gene

        Also accepted as a transcript if has exons as direct children
        )pbdoc")
    .value("ACCEPT_TRANSCRIPT", majiq::gff3::FeatureType::ACCEPT_TRANSCRIPT,
        R"pbdoc(
        Indicates feature type accepted as transcript if has exon children

        Unlike ACCEPT_GENE, requires parent that maps to ACCEPT_GENE, otherwise
        any exons belonging to it will be ignored
        )pbdoc")
    .value("REJECT_SILENT", majiq::gff3::FeatureType::REJECT_SILENT,
        R"pbdoc(
        Known feature type that is not accepted as transcript or gene

        Can be part of chain of ancestors of an exon for determining transcript
        gene, but:
        + cannot serve as transcript: direct exon children are ignored
        + cannot serve as gene: if top-level feature, all descendants
          are ignored
        Unlike REJECT_OTHER, these rejections as potential transcripts/genes
        will not be kept track of (rejected silently)
        )pbdoc")
    .value("REJECT_OTHER", majiq::gff3::FeatureType::REJECT_OTHER,
        R"pbdoc(
        Potentially unknown feature type not accepted as transcript or gene

        Can be part of chain of ancestors of an exon for determining transcript
        gene, but:
        + cannot serve as transcript: direct exon children are ignored
        + cannot serve as gene: if top-level feature, all descendants
          are ignored
        When these features are parents of exons or top-level features for
        descendant exons, will be accounted for
        )pbdoc")
    .value("HARD_SKIP", majiq::gff3::FeatureType::HARD_SKIP,
        R"pbdoc(
        Ignore this feature when building exon hierarchy

        If it is a parent of an exon or an ancestor of an accepted transcript,
        will raise exception when parsing annotation. To be used with care
        to ignore records that are completely unrelated to exons
        )pbdoc");
  py::bind_map<majiq::gff3::featuretype_map_t>(m, "GFF3Types");
  m.def("_default_gff3_types",
      []() { return majiq::gff3::default_gff3_types; });
  auto pySpliceGraph = py::class_<SpliceGraph>(m, "SpliceGraph",
      "Splicegraph managing exons, junctions, and introns within genes");

  auto pyExperimentThresholds
    = pyExperimentThresholds_t(m, "ExperimentThresholds",
        "Thresholds on intron and junction coverage for inclusion in SpliceGraph");
  auto pyIntronThresholdsGenerator
    = pyIntronThresholdsGenerator_t(m, "IntronThresholdsGenerator");

  init_pyExperimentThresholds(pyExperimentThresholds);
  init_pyIntronThresholdsGenerator(pyIntronThresholdsGenerator);
  init_Contigs(pyContigs);
  init_Genes(pyGenes);
  init_Exons(pyExons);
  init_GeneJunctions(pyGeneJunctions);
  init_GeneJunctionsAccumulator(pyGeneJunctionsAccumulator);
  init_GeneIntrons(pyGeneIntrons);
  init_SJJunctions(pySJJunctions);
  init_SJJunctionsBins(pySJJunctionsBins);
  init_SJIntrons(pySJIntrons);
  init_SJIntronsBins(pySJIntronsBins);
  init_pyGroupJunctionsGen(pyGroupJunctionsGen);
  init_pyPassedJunctionsGen(pyPassedJunctionsGen);
  init_pyGroupIntronsGen(pyGroupIntronsGen);
  init_PyEvents(pyEvents);
  init_pyEventsAlign(pyEventsAlign);
  init_PyEventsCoverage(pyEventsCoverage);
  init_PySpliceGraphReads(pySpliceGraphReads);
  init_pyExonConnections(pyExonConnections);
  init_PySimplifierGroup(pySimplifierGroup);
  init_SpliceGraph(pySpliceGraph);
}
