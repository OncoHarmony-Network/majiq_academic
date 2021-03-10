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
#include "internals/SpliceGraphReads.hpp"
#include "internals/SimplifierGroup.hpp"
#include "internals/Meta.hpp"

#include "internals/ExperimentThresholds.hpp"

// NOTE: right now we have a global PRNG, which is not threadsafe, but
// multithreading is really only used by htslib internally at this time
static majiq::rng_t rng;  // PRNG for new_majiq

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
using pySpliceGraphReads_t = pyClassShared_t<majiq::SpliceGraphReads>;

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
        "Underlying regions bin reads are defined over")
    .def_property_readonly("total_bins", &SJBinsT::total_bins,
        "the total number of bins possible (positions for junctions, too)")
    .def_property_readonly("bin_reads",
        [](py::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        const size_t offset = offsetof(BinReads, bin_reads);
        return ArrayFromVectorAndOffset<CountT, BinReads>(
            sj.reads(), offset, sj_obj);
        },
        "Number of reads for each bin")
    .def_property_readonly("bin_idx",
        [](py::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        const size_t offset = offsetof(BinReads, bin_idx);
        return ArrayFromVectorAndOffset<majiq::junction_pos_t, BinReads>(
            sj.reads(), offset, sj_obj);
        },
        "Bin index for the each bin reads")
    .def_property_readonly("_offsets",
        [](py::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            sj.offsets(), 0, sj_obj);
        },
        "Raw offsets for regions into bin reads")
    .def("numstacks",
        [](const SJBinsT& self,
          py::array_t<size_t> idx, py::array_t<majiq::real_t> pvalue) {
        auto f = [&self](size_t i, majiq::real_t p) {
          return self.numstacks(i, p); };
        return py::vectorize(f)(idx, pvalue);
        },
        "Get number of stacks for the specified regions given threshold",
        py::arg("region_idx"),
        py::arg("pvalue_threshold") = DEFAULT_BUILD_STACK_PVALUE)
    .def("numbins",
        [](const SJBinsT& self,
          py::array_t<size_t> idx,
          py::array_t<CountT> minreads) {
        auto f = [&self](size_t i, CountT r) {
          return self.numbins_minreads(i, r); };
        return py::vectorize(f)(idx, minreads);
        },
        "Number of bins for regions with more than specified number of reads",
        py::arg("region_idx"), py::arg("minreads"))
    .def("numreads",
        [](const SJBinsT& self,
          py::array_t<size_t> idx,
          py::array_t<majiq::junction_pos_t> num_stacks) {
        auto f = [&self](size_t i, majiq::junction_pos_t n) {
          return self.numreads(i, n); };
        return py::vectorize(f)(idx, num_stacks);
        },
        "Number of reads for regions given known number of stacks",
        py::arg("region_idx"), py::arg("num_stacks"))
    .def("__len__", &SJBinsT::size, "Total number of bin reads")
    .def(py::init([](
            std::shared_ptr<RegionsT> regions,
            py::array_t<majiq::junction_ct_t> _position_reads,
            py::array_t<majiq::junction_pos_t> _position,
            py::array_t<size_t> _offsets,
            majiq::junction_pos_t num_bins) {
          auto position_reads = _position_reads.unchecked<1>();
          auto position = _position.unchecked<1>();
          auto offsets = _offsets.unchecked<1>();
          std::vector<size_t> offsets_vec(offsets.shape(0));
          for (size_t i = 0; i < offsets_vec.size(); ++i) {
            offsets_vec[i] = offsets(i);
          }
          std::vector<BinReads> pr_vec(position_reads.shape(0));
          for (size_t i = 0; i < pr_vec.size(); ++i) {
            pr_vec[i] = BinReads{position(i), position_reads(i)};
          }
          return SJBinsT{regions, std::move(pr_vec), std::move(offsets_vec), num_bins};
          }),
        "Initialize bins over specified regions with per-bin coverage",
        py::arg("sj"),
        py::arg("bin_reads"), py::arg("bin_idx"), py::arg("_offsets"),
        py::arg("total_bins"));
}

template <const char* REGIONS_NC_GROUP,
         typename RegionsT,
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
      .def("find",
          [](const RegionsT& self,
            size_t gene_idx, position_t start, position_t end) {
          auto it = self.find(
              RegionT{(*self.parents_)[gene_idx], IntervalT{start, end}});
          return it == self.end() ? -1 : it - self.begin();
          },
          "Get index for specified region (or -1 if it doesn't exist)",
          py::arg("gene_idx"),
          py::arg("start"),
          py::arg("end"));
  }
  pyRegions
    .def("__len__", &RegionsT::size)
    .def_property_readonly("_parents", &RegionsT::parents,
        "Get parents object on which regions are defined (e.g. contigs, genes)")
    .def_property_readonly("_parent_idx_start",
        [](py::object& regions_obj) {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        return ArrayFromOffsetsVector<size_t>(
            regions.parent_idx_offsets(), true, regions_obj);
        },
        "First index into regions corresponding to associated parent")
    .def_property_readonly("_parent_idx_end",
        [](py::object& regions_obj) {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        return ArrayFromOffsetsVector<size_t>(
            regions.parent_idx_offsets(), false, regions_obj);
        },
        "One after last index into regions corresponding to associated parent")
    .def_property_readonly("start",
        [](py::object& regions_obj) -> py::array_t<position_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        "array[int] of starts for each feature")
    .def_property_readonly("end",
        [](py::object& regions_obj) -> py::array_t<position_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
        "array[int] of ends for each feature");
  if constexpr(majiq::detail::has_contig_field<RegionT>::value) {
    pyRegions.def_property_readonly("contig_idx",
        [](py::object& regions_obj) -> py::array_t<size_t> {
        RegionsT& regions = regions_obj.cast<RegionsT&>();
        const size_t offset = offsetof(RegionT, contig.idx_);
        return ArrayFromVectorAndOffset<size_t, RegionT>(
            regions.data(), offset, regions_obj);
        },
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
        "array[int] of indexes indicating gene feature belongs to");
  }
  if constexpr(
      std::is_same_v<decltype(std::declval<RegionT>().data),
                     majiq::detail::ConnectionData>) {
    pyRegions
      .def("_pass_all", &RegionsT::pass_all, "pass all connections")
      .def("_simplify_all",
          &RegionsT::simplify_all, "simplify all connections")
      .def("_unsimplify_all",
          &RegionsT::unsimplify_all, "unsimplify all connections")
      .def_property_readonly("denovo",
          [](py::object& regions_obj) -> py::array_t<bool> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.denovo);
          return ArrayFromVectorAndOffset<bool, RegionT>(
              regions.data(), offset, regions_obj);
          },
          "array[bool] indicating if connection was not found in annotations")
      .def_property_readonly("passed_build",
          [](py::object& regions_obj) -> py::array_t<bool> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.passed_build);
          return ArrayFromVectorAndOffset<bool, RegionT>(
              regions.data(), offset, regions_obj);
          },
          "array[bool] indicating if passed build criteria to be in LSV")
      .def_property_readonly("simplified",
          [](py::object& regions_obj) -> py::array_t<bool> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.simplified);
          return ArrayFromVectorAndOffset<bool, RegionT>(
              regions.data(), offset, regions_obj);
          },
          "array[bool] indicating if the connection is simplified")
      .def("connect_exons", &RegionsT::connect_exons,
          "Connect regions to specified exons, updataing {start,end}_exon_idx",
          py::arg("exons"))
      .def_property_readonly("connected_exons",
          [](RegionsT& self) -> std::optional<std::shared_ptr<majiq::Exons>> {
          using return_t = std::optional<std::shared_ptr<majiq::Exons>>;
          return self.is_connected()
            ? return_t{self.connected_exons()} : return_t{};
          },
          "Exons connected to (or None if not connected to any exons)")
      .def_property_readonly("start_exon_idx",
          [](py::object& regions_obj) -> py::array_t<size_t> {
          RegionsT& regions = regions_obj.cast<RegionsT&>();
          const size_t offset = offsetof(RegionT, data.start_exon_idx);
          return ArrayFromVectorAndOffset<size_t, RegionT>(
              regions.data(), offset, regions_obj);
          },
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
    .def(
        py::init([](py::list seqids) {
          auto result = Contigs::create();
          for (auto seqid : seqids) {
            result->add(Contig{seqid.cast<seqid_t>()});
          }
          return result;
        }),
        "Set up Contigs object using specified identifiers",
        py::arg("seqids"))
    .def_property_readonly("seqid", &Contigs::seqids,
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("__repr__", [](const Contigs& self) -> std::string {
        std::ostringstream oss;
        oss << self;
        return oss.str();
        })
    .def("__len__", &Contigs::size)
    .def("__contains__",
        [](const Contigs& s, seqid_t x) -> bool { return s.count(x) > 0; })
    .def("__getitem__",
        [](const Contigs& self, seqid_t x) { return self.get_idx(x); },
        "contig_idx for specified seqid", py::arg("seqid"));
}

void init_Genes(pyGenes_t& pyGenes) {
  using majiq_pybind::ArrayFromVectorAndOffset;
  using majiq::Contigs;
  using majiq::Genes;
  using majiq::position_t;
  using majiq::geneid_t;
  define_coordinates_properties<GENES_NC_GROUP>(pyGenes);
  pyGenes
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
            gene_vec.push_back(majiq::Gene{
                majiq::KnownContig{contig_idx(i), contigs},
                majiq::ClosedInterval{start(i), end(i)},
                static_cast<majiq::GeneStrandness>(strand(i)[0]),
                geneid[i].cast<majiq::geneid_t>(),
                genename[i].cast<majiq::genename_t>()});
          }
          return Genes::create(contigs, std::move(gene_vec));
        }),
        "Create Genes object using Contigs object and arrays defining genes",
        py::arg("contigs"), py::arg("contig_idx"),
        py::arg("start"), py::arg("end"), py::arg("strand"),
        py::arg("gene_id"), py::arg("gene_name"))
    .def_property_readonly("gene_id", &majiq::Genes::geneids,
        "Sequence[str] of gene ids in order matching gene_idx")
    .def_property_readonly("gene_name", &majiq::Genes::genenames,
        "Sequence[str] of gene names in order matching gene_idx")
    .def("__repr__", [](const majiq::Genes& self) -> std::string {
        std::ostringstream oss;
        oss << "Genes<" << self.size() << " total>";
        return oss.str();
        })
    .def("__contains__",
        [](const Genes& self, geneid_t x) { return self.count(x) > 0; })
    .def("__getitem__",
        [](const Genes& self, geneid_t x) { return self.get_idx(x); },
        "gene_idx for specified gene_id", py::arg("gene_id"));
}

void init_Exons(pyExons_t& pyExons) {
  using majiq::Exons;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties<EXONS_NC_GROUP>(pyExons);
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
            exon_vec.push_back(majiq::Exon{
                majiq::KnownGene{gene_idx(i), genes},
                majiq::ClosedInterval{start(i), end(i)},
                majiq::ClosedInterval{ann_start(i), ann_end(i)}});
          }
          return std::make_shared<Exons>(genes, std::move(exon_vec));
        }),
        "Create Exons object using Genes and info about each exon",
        py::arg("genes"), py::arg("gene_idx"),
        py::arg("start"), py::arg("end"),
        py::arg("annotated_start"), py::arg("annotated_end"))
    .def_property_readonly("annotated_start",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, data.start);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of annotated exon starts")
    .def_property_readonly("annotated_end",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, data.end);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of annotated exon ends")
    .def("__repr__", [](const Exons& self) -> std::string {
        std::ostringstream oss;
        oss << "Exons<" << self.size() << " total>";
        return oss.str();
        });
}

void init_GeneJunctions(pyGeneJunctions_t& pyGeneJunctions) {
  using majiq::GeneJunctions;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties<JUNCTIONS_NC_GROUP>(pyGeneJunctions);
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
        "Create GeneJunctions using Genes and arrays defining each junction",
        py::arg("genes"),
        py::arg("gene_idx"), py::arg("start"), py::arg("end"),
        py::arg("denovo"), py::arg("passed_build"), py::arg("simplified"))
    .def("__repr__", [](const GeneJunctions& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneJunctions<" << self.size() << " total>";
        return oss.str();
        });
}

void init_SJIntrons(pySJIntrons_t& pySJIntrons) {
  using majiq::position_t;
  using majiq::Contigs;
  using majiq::SJIntron;
  using majiq::SJIntrons;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties<SJ_INTRONS_NC_GROUP>(pySJIntrons);
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
            result[i] = SJIntron{
              majiq::KnownContig{contig_idx(i), contigs},
              majiq::ClosedInterval{start(i), end(i)},
              static_cast<majiq::GeneStrandness>(strand(i)[0]),
              annotated(i)};
          }
          return std::make_shared<SJIntrons>(contigs, std::move(result));
          }),
        "Initialize SJIntrons from contigs and numpy arrays",
        py::arg("contigs"),
        py::arg("contig_idx"),
        py::arg("start"),
        py::arg("end"),
        py::arg("strand"),
        py::arg("annotated"))
    .def_static("from_exons_and_introns", &SJIntrons::FromGeneExonsAndIntrons,
        "Construct sj introns for input exons/introns",
        py::arg("exons"), py::arg("introns"), py::arg("stranded"))
    .def_property_readonly("annotated",
        [](py::object& introns_obj) -> py::array_t<bool> {
        SJIntrons& introns = introns_obj.cast<SJIntrons&>();
        const size_t offset = offsetof(majiq::SJIntron, data.annotated_);
        return ArrayFromVectorAndOffset<bool, majiq::SJIntron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if intron is annotated (exon in annotation)")
    .def("__repr__", [](const majiq::SJIntrons& self) -> std::string {
        std::ostringstream oss;
        oss << "SJIntrons<" << self.size() << " total>";
        return oss.str();
        });
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
        "Underlying introns")
    .def_property_readonly("_junctions",
        &SpliceGraphReads::junctions,
        "Underlying junctions")
    .def_property_readonly("introns_reads",
        [](py::object& self_obj) {
        SpliceGraphReads& self = self_obj.cast<SpliceGraphReads&>();
        return ArrayFromVectorAndOffset<majiq::real_t, majiq::real_t>(
            self.introns_reads(), 0, self_obj);
        },
        "Raw readrates for each intron")
    .def_property_readonly("junctions_reads",
        [](py::object& self_obj) {
        SpliceGraphReads& self = self_obj.cast<SpliceGraphReads&>();
        return ArrayFromVectorAndOffset<majiq::real_t, majiq::real_t>(
            self.junctions_reads(), 0, self_obj);
        },
        "Raw readrates for each junction")
    .def_static("from_sj", &SpliceGraphReads::FromSJ,
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
        "Underlying exon connections being unsimplified")
    .def_property_readonly("introns_passed_src",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, src_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.introns_passed(), offset, self_obj);
        },
        "Number of experiments with evidence for unsimplification as source for introns")
    .def_property_readonly("introns_passed_dst",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, dst_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.introns_passed(), offset, self_obj);
        },
        "Number of experiments with evidence for unsimplification as target for introns")
    .def_property_readonly("junctions_passed_src",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, src_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.junctions_passed(), offset, self_obj);
        },
        "Number of experiments with evidence for unsimplification as source for junctions")
    .def_property_readonly("junctions_passed_dst",
        [](py::object& self_obj) {
        SimplifierGroup& self = self_obj.cast<SimplifierGroup&>();
        size_t offset = offsetof(SimplifierCount, dst_ct_unsimplify);
        return ArrayFromVectorAndOffset<size_t, SimplifierCount>(
            self.junctions_passed(), offset, self_obj);
        },
        "Number of experiments with evidence for unsimplification as target for junctions")
    .def_property_readonly("num_experiments",
        &SimplifierGroup::num_experiments,
        "Number of experiments in current simplifier group")
    .def("add_experiment", &SimplifierGroup::AddExperiment,
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
        "Unsimplify introns/junctions with evidence, reset for next group",
        py::arg("simplifier_min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def(py::init<const std::shared_ptr<ExonConnections>&>(),
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
    .def_static("from_sj",
        [](
          const std::shared_ptr<Events>& events,
          const majiq::SJJunctionsBins& sj_junctions,
          const majiq::SJIntronsBins& sj_introns,
          size_t num_bootstraps,
          majiq::real_t pvalue_threshold) {
        return EventsCoverage::FromSJ(events, sj_junctions, sj_introns,
            num_bootstraps, rng, pvalue_threshold);
        },
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
        "Readrate of each event connection after stacks removed")
    .def_property_readonly("numbins",
        [](py::object& self_obj) -> py::array_t<majiq::real_t> {
        EventsCoverage& self = self_obj.cast<EventsCoverage&>();
        const size_t offset = offsetof(CoverageSummary, numbins);
        return ArrayFromVectorAndOffset<majiq::real_t, CoverageSummary>(
            self.summaries(), offset, self_obj);
        },
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
        "Bootstrapped read coverage or each connection after stacks removed")
    .def_property_readonly("_events", &EventsCoverage::events,
        "Events for which the coverage information is defined")
    .def("__len__", &EventsCoverage::num_connections);
}

void init_PyEvents(pyEvents_t& pyEvents) {
  using majiq::Event;
  using majiq::Events;
  using majiq::ConnectionIndex;
  using majiq_pybind::ArrayFromOffsetsVector;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyEvents
    .def_property_readonly("ref_exon_idx",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(Event, ref_exon_idx_);
        return ArrayFromVectorAndOffset<size_t, Event>(
            self.events(), offset, self_obj); },
        "Event reference exon")
    .def_property_readonly("event_type",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(Event, type_);
        return ArrayFromVectorAndOffset<std::array<char, 1>, Event>(
            self.events(), offset, self_obj); },
        "Event type")
    .def_property_readonly("connection_idx_start",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.connection_offsets(), true, self_obj); },
        "First index into event connections for each event")
    .def_property_readonly("connection_idx_end",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        return ArrayFromOffsetsVector<size_t>(
            self.connection_offsets(), false, self_obj); },
        "One after last index into event connections for each event")
    .def("df",
        [](py::object& self) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(self, "event_idx",
            {"ref_exon_idx", "event_type", "connection_idx_start",
            "connection_idx_end"}); },
        "View on event information as xarray Dataset")
    .def_property_readonly("is_intron",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(ConnectionIndex, is_intron_);
        return ArrayFromVectorAndOffset<bool, ConnectionIndex>(
            self.connections(), offset, self_obj); },
        "Event connection is intron (false --> is junction)")
    .def_property_readonly("idx",
        [](py::object& self_obj) {
        Events& self = self_obj.cast<Events&>();
        const size_t offset = offsetof(ConnectionIndex, idx_);
        return ArrayFromVectorAndOffset<size_t, ConnectionIndex>(
            self.connections(), offset, self_obj); },
        "Event connection index into corresponding Introns or Junctions")
    .def("connections_df",
        [](py::object& self) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(self, "connection_idx",
            {"is_intron", "idx"}); },
        "View on event connection information as xarray Dataset")
    .def_property_readonly("num_events", &Events::num_events)
    .def_property_readonly("num_connections", &Events::num_connections)
    .def_property_readonly("num_junctions", &Events::num_junctions)
    .def_property_readonly("num_introns", &Events::num_introns)
    .def("__len__", &Events::size);
}

void init_pyExonConnections(pyExonConnections_t& pyExonConnections) {
  using majiq::Event;
  using majiq::EventType;
  using majiq::ExonConnections;
  using ExonsPtrT = std::shared_ptr<majiq::Exons>;
  using IntronsPtrT = std::shared_ptr<majiq::GeneIntrons>;
  using JunctionsPtrT = std::shared_ptr<majiq::GeneJunctions>;
  pyExonConnections
    .def(
        py::init<const ExonsPtrT&, const IntronsPtrT&, const JunctionsPtrT&>(),
        "Track junctions/introns by exons/events assocciated with them",
        py::arg("exons"),
        py::arg("introns"),
        py::arg("junctions"))
    .def("lsvs", &ExonConnections::LSVEvents, "Construct LSV Events")
    .def("constitutive", &ExonConnections::ConstitutiveEvents,
        "Construct Constitutive Events")
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
          events.emplace_back(
              _exon_idx(i),
              _is_source(i) ? EventType::SRC_EVENT : EventType::DST_EVENT);
        }
        return self.CreateEvents(std::move(events));
        },
        "Construct events for specified exons/directions",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("has_intron",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) {
          return self.has_intron(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate if events have introns or not",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("has_intron",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          return self.has_intron(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate if events have introns or not",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("event_size",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<size_t> {
        auto f = [&self](size_t idx, bool is_src) -> size_t {
          return self.event_size(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate size of event",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("passed",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          return self.passed(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate if event was passed",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("redundant",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          return self.redundant(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate if event was redundant",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("is_LSV",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          return self.is_LSV(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate if event is LSV",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("is_constitutive",
        [](const ExonConnections& self,
          py::array_t<size_t> exon_idx,
          py::array_t<bool> is_source) -> py::array_t<bool> {
        auto f = [&self](size_t idx, bool is_src) -> bool {
          return self.is_constitutive(Event{
              idx, is_src ? EventType::SRC_EVENT : EventType::DST_EVENT});
        };
        return py::vectorize(f)(exon_idx, is_source);
        },
        "Indicate if event is constitutive",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("event_id",
        [](const ExonConnections& self, size_t exon_idx, bool is_source) {
        return self.id(Event{
            exon_idx, is_source ? EventType::SRC_EVENT : EventType::DST_EVENT});
        },
        "Return identifier for event",
        py::arg("exon_idx"), py::arg("is_source"))
    .def("event_description",
        [](const ExonConnections& self, size_t exon_idx, bool is_source) {
        return self.description(Event{
            exon_idx, is_source ? EventType::SRC_EVENT : EventType::DST_EVENT});
        },
        "Return description for event",
        py::arg("exon_idx"), py::arg("is_source"));
}

void init_SJIntronsBins(pySJIntronsBins_t& pySJIntronsBins) {
  using majiq::SJIntronsBins;
  define_sjbins_properties(pySJIntronsBins);
  pySJIntronsBins
    .def_static("from_bam", &SJIntronsBins::FromBam,
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
  define_coordinates_properties<INTRONS_NC_GROUP>(pyGeneIntrons);
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
        "Create GeneIntrons using Genes and arrays defining each intron",
        py::arg("genes"),
        py::arg("gene_idx"), py::arg("start"), py::arg("end"),
        py::arg("denovo"), py::arg("passed_build"), py::arg("simplified"))
    .def("build_group", [](std::shared_ptr<GeneIntrons>& gene_introns) {
        return majiq::GroupIntronsGenerator(gene_introns);
        },
        "Create build group to update passed introns in place")
    .def("filter_passed", &GeneIntrons::FilterPassed,
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
    .def("potential_introns", &GeneIntrons::PotentialIntrons,
        "Get potential gene introns from exons keeping annotations from self",
        py::arg("exons"))
    .def("__repr__", [](const GeneIntrons& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneIntrons<" << self.size() << " total>";
        return oss.str();
        });
}

void init_SJJunctions(pySJJunctions_t& pySJJunctions) {
  using majiq::position_t;
  using majiq::SJJunctions;
  using majiq_pybind::ArrayFromVectorAndOffset;
  define_coordinates_properties<SJ_JUNCTIONS_NC_GROUP>(pySJJunctions);
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
            sj_vec[i] = majiq::SJJunction{
              majiq::KnownContig{contig_idx(i), contigs},
              majiq::OpenInterval{start(i), end(i)},
              static_cast<majiq::GeneStrandness>(strand(i)[0])};
          }
          return std::make_shared<majiq::SJJunctions>(
              contigs, std::move(sj_vec));
        }),
        "Create SJJunctions object from contigs and arrays",
        py::arg("contigs"),
        py::arg("contig_idx"), py::arg("start"), py::arg("end"),
        py::arg("strand"))
    .def_property_readonly("_contigs", &SJJunctions::parents,
        "Underlying contigs corresponding to contig_idx");
}
void init_SJJunctionsBins(pySJJunctionsBins_t& pySJJunctionsBins) {
  using majiq::SJJunctionsBins;
  define_sjbins_properties(pySJJunctionsBins);
  pySJJunctionsBins
    .def_static("from_bam", &SJJunctionsBins::FromBam,
        R"pbdoc(
        Load junctions and per-position counts for an aligned BAM file

        Parameters
        ----------
        bam_path: str
            Path for input BAM fille
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
      "Accumulate SJJunctions for build group for input junctions/exons",
      py::arg("junctions"),
      py::arg("exons"))
    .def("add_experiment", &majiq::GroupJunctionsGenerator::AddExperiment,
        "Increment count of passed junctions from input experiment",
        py::arg("sjp"),
        py::arg("thresholds") = DEFAULT_THRESHOLDS,
        py::arg("add_denovo") = DEFAULT_BUILD_DENOVO_JUNCTIONS)
    .def("pass_known_inplace",
        &majiq::GroupJunctionsGenerator::UpdateKnownInplace,
        "Update known junctions with their build status (ignore denovos)",
        py::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def_property_readonly("num_experiments",
        &majiq::GroupJunctionsGenerator::num_experiments,
        "Number of experiments that have been added to this group")
    .def_property_readonly("num_known",
        &majiq::GroupJunctionsGenerator::num_annotated,
        "Number of junctions known when constructed")
    .def_property_readonly("num_denovo",
        &majiq::GroupJunctionsGenerator::num_denovo,
        "Number of denovo junctions passing experiment-filters for 1+ inputs")
    .def("__len__", &majiq::GroupJunctionsGenerator::size);
}

void init_pyGroupIntronsGen(pyGroupIntronsGen_t& pyGroupIntronsGen) {
  using majiq::GroupIntronsGenerator;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyGroupIntronsGen
    .def(
        py::init<const std::shared_ptr<majiq::GeneIntrons>&>(),
        "Accumulate SJIntronsBins for build group of introns",
        py::arg("gene_introns"))
    .def_property_readonly("num_experiments",
        &GroupIntronsGenerator::num_experiments,
        "Number of experiments in current group")
    .def_property_readonly("_introns", &GroupIntronsGenerator::introns)
    .def_property_readonly("num_passed",
        [](py::object& self_obj) -> py::array_t<size_t> {
        GroupIntronsGenerator& self = self_obj.cast<GroupIntronsGenerator&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            self.num_passed(), 0, self_obj);
        },
        "Number of experiments or which each intron has passed")
    .def("df",
        [](py::object& self) -> py::object {
        auto base = self.attr("_introns").attr("df")();
        auto get_xr = [&self](py::str key) {
          return py::module_::import("xarray").attr("DataArray")(
              self.attr(key), py::arg("dims") = "intron_idx");
        };
        return base.attr("assign_coords")(
            py::arg("num_passed") = get_xr("num_passed"));
        },
        "View of gene introns and how many times they have passed")
    .def("__len__", &GroupIntronsGenerator::size)
    .def("add_experiment", &GroupIntronsGenerator::AddExperiment,
        "Add SJIntronsBins to build group",
        py::arg("sj"),
        py::arg("thresholds") = DEFAULT_THRESHOLDS)
    .def("update_introns", &GroupIntronsGenerator::UpdateInplace,
        "Pass introns that pass min-experiments in place, reset for next group",
        py::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS);
}

void init_pyPassedJunctionsGen(pyPassedJunctionsGen_t& pyPassedJunctionsGen) {
  pyPassedJunctionsGen
    .def(
      py::init<const std::shared_ptr<majiq::GeneJunctions>&>(),
      R"pbdoc(
      Accumulator of GroupJunctionsGenerator for different build groups
      )pbdoc",
      py::arg("junctions"))
    .def("add_group", &majiq::PassedJunctionsGenerator::AddGroup,
        "Combine passed junctions from build group of experiments",
        py::arg("group"),
        py::arg("min_experiments") = DEFAULT_BUILD_MINEXPERIMENTS)
    .def("get_passed", &majiq::PassedJunctionsGenerator::PassedJunctions,
        "Get static, array-based representation of passed junctions")
    .def_property_readonly("num_known",
        &majiq::PassedJunctionsGenerator::num_annotated,
        "Number of junctions known when constructed")
    .def_property_readonly("num_denovo",
        &majiq::PassedJunctionsGenerator::num_denovo,
        "Number of denovo junctions passing filters, previously not known")
    .def("__len__", &majiq::PassedJunctionsGenerator::size);
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
        "Create splicegraph from input GFF3 file",
        py::arg("gff3_path"),
        py::arg("process_ir") = DEFAULT_BUILD_PROCESS_IR,
        py::arg_v("gff3_types", majiq::gff3::default_gff3_types,
            "new_majiq._default_gff3_types()"))
    // XXX debug
    .def_static("infer_exons", &SpliceGraph::InferExons,
        "Infer exons from base annotated exons and junctions",
        py::arg("base_exons"), py::arg("junctions"))
    .def("make_group_junctions", &SpliceGraph::MakeGroupGenerator,
        "Create GroupJunctionsGenerator for the splicegraph junctions/exons")
    .def("make_build_junctions", &SpliceGraph::MakePassedGenerator,
        "Create PassedJunctionsGenerator for splicegraph junctions")
    .def("build_junctions", &SpliceGraph::BuildJunctionExons,
        "New splicegraph with updated junctions/exons using build junctions",
        py::arg("build_junctions"))
    .def("_pass_all",
        [](py::object& sg) {
        sg.attr("_junctions").attr("_pass_all")();
        sg.attr("_introns").attr("_pass_all")(); },
        "Pass all junctions and introns in the splicegraph")
    .def("_simplify_all",
        [](py::object& sg) {
        sg.attr("_junctions").attr("_simplify_all")();
        sg.attr("_introns").attr("_simplify_all")(); },
        "Simplify all junctions and introns in the splicegraph")
    .def("_unsimplify_all",
        [](py::object& sg) {
        sg.attr("_junctions").attr("_unsimplify_all")();
        sg.attr("_introns").attr("_unsimplify_all")(); },
        "Unsimplify all junctions and introns in the splicegraph")
    .def("close_to_annotated_exon",
        [](SpliceGraph& sg, size_t gene_idx, position_t x, bool to_following) {
        majiq::KnownGene g = (*sg.genes())[gene_idx];
        const majiq::Exons& exons = *sg.exons();
        return to_following
          ? majiq::detail::CloseToFollowingAnnotatedExon(exons, g, x)
          : majiq::detail::CloseToPrecedingAnnotatedExon(exons, g, x);
        },
        "True if position close to following/preceding annotated exon in gene",
        py::arg("gene_idx"),
        py::arg("x"),
        py::arg("to_following") = true)
    // access underlying data
    .def_property_readonly("_exons", &SpliceGraph::exons,
        "Access the splicegraph's exons")
    .def_property_readonly("_introns", &SpliceGraph::introns,
        "Access the splicegraph's introns")
    .def_property_readonly("_junctions", &SpliceGraph::junctions,
        "Access the splicegraph's junctions")
    .def_property_readonly("_genes", &SpliceGraph::genes,
        "Access the splicegraph's genes")
    .def_property_readonly("_contigs", &SpliceGraph::contigs,
        "Access the splicegraph's contigs")
    .def_property_readonly("_exon_connections", &SpliceGraph::exon_connections,
        "Access the splicegraph's exon connections")
    // get sj introns on which coverage may be read
    .def("sj_introns", [](SpliceGraph& sg, bool stranded) {
        return majiq::SJIntrons::FromGeneExonsAndIntrons(
            *sg.exons(), *sg.introns(), stranded);
        },
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
        });
}

void init_pyExperimentThresholds(
    pyExperimentThresholds_t& pyExperimentThresholds) {
  using majiq::ExperimentThresholds;
  using majiq::junction_ct_t;
  using majiq::junction_pos_t;
  using majiq::real_t;
  pyExperimentThresholds
    .def(py::init<junction_ct_t, junction_ct_t, junction_pos_t, real_t, real_t, real_t>(),
        "per-experiment thresholds for junctions",
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
        "Minimum number of reads for an annotated junction to pass")
    .def_property_readonly("mindenovo",
        [](const ExperimentThresholds& x) { return x.mindenovo_; },
        "Minimum number of reads for a denovo junction to pass")
    .def_property_readonly("minpos",
        [](const ExperimentThresholds& x) { return x.minpos_; },
        "Minimum number of nonzero positions for a junction to pass")
    .def_property_readonly("max_pctbins",
        [](const ExperimentThresholds& x) { return x.max_pctbins_; },
        "Maximum percentage of bins to require coverage in for intron to pass")
    .def_property_readonly("junction_acceptance_probability",
        [](const ExperimentThresholds& x) {
        return x.junction_acceptance_probability_; },
        "Set intron thresholds to match junction readrate with probability")
    .def_property_readonly("intron_acceptance_probability",
        [](const ExperimentThresholds& x) {
        return x.intron_acceptance_probability_; },
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
        })
    .def("intron_thresholds_generator",
        &ExperimentThresholds::intron_thresholds_generator,
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
        "Get intron thresholds for specified intron lengths",
        py::arg("intron_lengths"));
}

void init_SpliceGraphAll(py::module_& m) {
  rng = majiq::rng_t{};  // initialize PRNG
  m.def("set_seed",
      [](size_t x) { rng.seed(x); },
      "set random seed for new_majiq random number generation",
      py::arg("seed"));

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
  auto pyEvents = pyEvents_t(
      m, "Events", "Events from reference exon with junctions/introns");
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
    = pyExperimentThresholds_t(m, "ExperimentThresholds");
  auto pyIntronThresholdsGenerator
    = pyIntronThresholdsGenerator_t(m, "IntronThresholdsGenerator");

  init_pyExperimentThresholds(pyExperimentThresholds);
  init_pyIntronThresholdsGenerator(pyIntronThresholdsGenerator);
  init_Contigs(pyContigs);
  init_Genes(pyGenes);
  init_Exons(pyExons);
  init_GeneJunctions(pyGeneJunctions);
  init_GeneIntrons(pyGeneIntrons);
  init_SJJunctions(pySJJunctions);
  init_SJJunctionsBins(pySJJunctionsBins);
  init_SJIntrons(pySJIntrons);
  init_SJIntronsBins(pySJIntronsBins);
  init_pyGroupJunctionsGen(pyGroupJunctionsGen);
  init_pyPassedJunctionsGen(pyPassedJunctionsGen);
  init_pyGroupIntronsGen(pyGroupIntronsGen);
  init_PyEvents(pyEvents);
  init_PyEventsCoverage(pyEventsCoverage);
  init_PySpliceGraphReads(pySpliceGraphReads);
  init_pyExonConnections(pyExonConnections);
  init_PySimplifierGroup(pySimplifierGroup);
  init_SpliceGraph(pySpliceGraph);
}
