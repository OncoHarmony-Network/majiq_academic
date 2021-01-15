/**
 * pyContigs.cpp
 *
 * Set up python bindings to contigs
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <memory>
#include <cstddef>
#include <stdexcept>

#include "pybind_utils.hpp"
#include "internals/Contigs.hpp"
#include "internals/Genes.hpp"
#include "internals/Exons.hpp"
#include "internals/GeneJunctions.hpp"
#include "internals/Introns.hpp"
#include "internals/SpliceGraph.hpp"
#include "internals/GFF3.hpp"


constexpr char CONTIGS_NC_GROUP[] = "contigs";
constexpr char GENES_NC_GROUP[] = "genes";
constexpr char EXONS_NC_GROUP[] = "exons";
constexpr char JUNCTIONS_NC_GROUP[] = "junctions";
constexpr char INTRONS_NC_GROUP[] = "introns";

namespace py = pybind11;


std::shared_ptr<majiq::Contigs> ContigsFromNetcdf(py::str netcdf_path) {
  auto xr_contigs = majiq_pybind::OpenXarrayDataset(
      netcdf_path, py::str(CONTIGS_NC_GROUP));
  auto result = std::make_shared<majiq::Contigs>();
  py::list seqids = xr_contigs
    .attr("__getitem__")("seqid").attr("values").attr("tolist")();
  for (auto seqid : seqids) {
    result->add(seqid.cast<majiq::seqid_t>());
  }
  return result;
}

std::shared_ptr<majiq::Genes> GenesFromNetcdf(
    std::shared_ptr<majiq::Contigs> contigs, py::str netcdf_path) {
  using majiq::position_t;
  using majiq::geneid_t;
  using majiq::genename_t;
  auto xr_genes = majiq_pybind::OpenXarrayDataset(
      netcdf_path, py::str(GENES_NC_GROUP));
  auto get_array = [&xr_genes](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_genes.attr("__getitem__")(key));
  };
  py::array_t<size_t> _contig_idx = get_array("contig_idx");
  auto contig_idx = _contig_idx.unchecked<1>();
  py::array_t<position_t> _start = get_array("start");
  auto start = _start.unchecked<1>();
  py::array_t<position_t> _end = get_array("end");
  auto end = _end.unchecked<1>();
  py::array_t<std::array<char, 1>> _strand = get_array("strand");
  auto strand = _strand.unchecked<1>();
  py::list geneid
    = xr_genes.attr("__getitem__")("gene_id").attr("values").attr("tolist")();
  py::list genename
    = xr_genes.attr("__getitem__")("gene_name").attr("values").attr("tolist")();
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
  return std::make_shared<majiq::Genes>(std::move(gene_vec));
}

std::shared_ptr<majiq::Exons> ExonsFromNetcdf(
    std::shared_ptr<majiq::Genes> genes, py::str netcdf_path) {
  using majiq::position_t;
  auto xr_exons = majiq_pybind::OpenXarrayDataset(
      netcdf_path, py::str(EXONS_NC_GROUP));
  auto get_array = [&xr_exons](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_exons.attr("__getitem__")(key));
  };
  // extract usable values from dataset
  py::array_t<size_t> arr_gene_idx = get_array("gene_idx");
  py::array_t<position_t> arr_start = get_array("start");
  py::array_t<position_t> arr_end = get_array("end");
  py::array_t<position_t> arr_ann_start = get_array("annotated_start");
  py::array_t<position_t> arr_ann_end = get_array("annotated_end");
  // unchecked accesses to numpy array
  auto gene_idx = arr_gene_idx.unchecked<1>();
  auto start = arr_start.unchecked<1>();
  auto end = arr_end.unchecked<1>();
  auto ann_start = arr_ann_start.unchecked<1>();
  auto ann_end = arr_ann_end.unchecked<1>();
  // create vector of genes matching input arrays
  std::vector<majiq::Exon> exon_vec{};
  exon_vec.reserve(gene_idx.shape(0));
  for (py::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
    exon_vec.push_back(majiq::Exon{
        majiq::KnownGene{gene_idx(i), genes},
        majiq::ClosedInterval{start(i), end(i)},
        majiq::ClosedInterval{ann_start(i), ann_end(i)}});
  }
  return std::make_shared<majiq::Exons>(std::move(exon_vec));
}

template <class Connections, const char group_str[]>
std::shared_ptr<Connections> ConnectionsFromNetcdf(
    std::shared_ptr<majiq::Genes> genes, py::str netcdf_path) {
  using IntervalT = typename Connections::IntervalT;
  using RegionT = typename Connections::value_type;
  using majiq::position_t;
  auto xr_connections = majiq_pybind::OpenXarrayDataset(
      netcdf_path, py::str(group_str));
  auto get_array = [&xr_connections](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_connections.attr("__getitem__")(key));
  };
  // extract usable values from dataset
  py::array_t<size_t> arr_gene_idx = get_array("gene_idx");
  py::array_t<position_t> arr_start = get_array("start");
  py::array_t<position_t> arr_end = get_array("end");
  py::array_t<bool> arr_denovo = get_array("denovo");
  py::array_t<bool> arr_passed_build = get_array("passed_build");
  py::array_t<bool> arr_simplified = get_array("simplified");
  // unchecked accesses to numpy array
  auto gene_idx = arr_gene_idx.unchecked<1>();
  auto start = arr_start.unchecked<1>();
  auto end = arr_end.unchecked<1>();
  auto denovo = arr_denovo.unchecked<1>();
  auto passed_build = arr_passed_build.unchecked<1>();
  auto simplified = arr_simplified.unchecked<1>();
  // create vector of genes matching input arrays
  using majiq::KnownGene;

  std::vector<RegionT> connection_vec{};
  connection_vec.reserve(gene_idx.shape(0));
  for (py::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
    connection_vec.push_back(RegionT{
        majiq::KnownGene{gene_idx(i), genes}, IntervalT{start(i), end(i)},
        denovo(i), passed_build(i), simplified(i)});
  }
  return std::make_shared<Connections>(std::move(connection_vec));
}


void init_Contigs(
    py::class_<majiq::Contigs, std::shared_ptr<majiq::Contigs>>& pyContigs) {
  using majiq::seqid_t;
  using majiq::Contigs;
  using namespace py::literals;
  pyContigs
    .def_static("from_netcdf",
        [](py::str x) { return ContigsFromNetcdf(x); },
        "Load contigs from netcdf file", py::arg("netcdf_path"))
    .def("to_netcdf",
        [](py::object& self, py::str out, py::str mode) {
        using namespace py::literals;
        self.attr("df")().attr("to_netcdf")(
            out, "mode"_a = mode, "group"_a = py::str(CONTIGS_NC_GROUP));
        return;
        },
        "Save contigs to netcdf file",
        py::arg("netcdf_path"), py::arg("mode"))
    .def_property_readonly("seqid", &Contigs::seqids,
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("df",
        [](py::object& contigs) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(contigs, "contig_idx", {"seqid"});
        },
        "View on contig information as xarray Dataset")
    .def("__repr__", [](const Contigs& self) -> std::string {
        std::ostringstream oss;
        oss << self;
        return oss.str();
        })
    .def("__len__", &Contigs::size)
    .def("__contains__",
        [](const Contigs& s, seqid_t x) -> bool { return s.contains(x); });
}

void init_Genes(
    py::class_<majiq::Genes, std::shared_ptr<majiq::Genes>>& pyGenes) {
  using majiq_pybind::ArrayFromVectorAndOffset;
  using majiq::Contigs;
  using majiq::Genes;
  using majiq::position_t;
  using majiq::geneid_t;
  pyGenes
    .def_static("from_netcdf",
        [](py::str x, std::shared_ptr<Contigs> contigs) {
        return GenesFromNetcdf(contigs, x);
        },
        "Load genes from netcdf file",
        py::arg("netcdf_path"), py::arg("contigs"))
    .def("to_netcdf",
        [](py::object& self, py::str out, py::str mode) {
        using namespace py::literals;
        self.attr("df")().attr("to_netcdf")(
            out, "mode"_a = mode, "group"_a = py::str(GENES_NC_GROUP));
        return;
        },
        "Save contigs to netcdf file",
        py::arg("netcdf_path"), py::arg("mode"))
    .def_property_readonly("strand",
        [](py::object& genes_obj) -> py::array_t<std::array<char, 1>> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetof(majiq::Gene, strand);
        return ArrayFromVectorAndOffset<std::array<char, 1>, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[char] of characters indicating strand of each gene")
    .def_property_readonly("contig_idx",
        [](py::object& genes_obj) -> py::array_t<size_t> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetof(majiq::Gene, contig.contig_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[int] of indexes indicating contig gene belongs to")
    .def_property_readonly("start",
        [](py::object& genes_obj) -> py::array_t<position_t> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetof(majiq::Gene, interval.start);
        return ArrayFromVectorAndOffset<position_t, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[int] of gene starts matching gene_idx")
    .def_property_readonly("end",
        [](py::object& genes_obj) -> py::array_t<position_t> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetof(majiq::Gene, interval.end);
        return ArrayFromVectorAndOffset<position_t, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[int] of gene ends matching gene_idx")
    .def_property_readonly("gene_id", &majiq::Genes::geneids,
        "Sequence[str] of gene ids in order matching gene_idx")
    .def_property_readonly("gene_name", &majiq::Genes::genenames,
        "Sequence[str] of gene names in order matching gene_idx")
    .def("df",
        [](py::object& genes) -> py::object {
        using majiq_pybind::XarrayDatasetFromObject;
        return XarrayDatasetFromObject(genes, "gene_idx",
            {"contig_idx", "start", "end", "strand", "gene_id", "gene_name"});
        },
        "View on gene information as xarray Dataset")
    .def_property_readonly("is_sorted", &majiq::Genes::is_sorted,
        "True if the genes sorted in contig/coordinate order")
    .def("__repr__", [](const majiq::Genes& self) -> std::string {
        std::ostringstream oss;
        oss << "Genes<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &majiq::Genes::size)
    .def("__contains__",
        [](const Genes& self, geneid_t x) { return self.contains(x); })
    .def("__getitem__",
        [](const Genes& self, geneid_t x) { return self.get_gene_idx(x); },
        "gene_idx for specified gene_id", py::arg("gene_id"));
}

void init_Exons(
    py::class_<majiq::Exons, std::shared_ptr<majiq::Exons>>& pyExons) {
  using majiq::Exons;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyExons
    .def_static("from_netcdf",
        [](py::str x, std::shared_ptr<majiq::Genes> genes) {
        return ExonsFromNetcdf(genes, x);
        },
        "Load exons from netcdf file",
        py::arg("netcdf_path"), py::arg("genes"))
    .def("to_netcdf",
        [](py::object& self, py::str out, py::str mode) {
        using namespace py::literals;
        self.attr("df")().attr("to_netcdf")(
            out, "mode"_a = mode, "group"_a = py::str(EXONS_NC_GROUP));
        return;
        },
        "Save exons to netcdf file",
        py::arg("netcdf_path"), py::arg("mode"))
    .def_property_readonly("gene_idx",
        [](py::object& exons_obj) -> py::array_t<size_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, gene.gene_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of indexes indicating gene exon belongs to")
    .def_property_readonly("start",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of exon starts")
    .def_property_readonly("end",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetof(majiq::Exon, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of exon ends")
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
    .def("df",
        [](py::object& exons) -> py::object {
        return majiq_pybind::XarrayDatasetFromObject(exons, "exon_idx",
            {"gene_idx", "start", "end", "annotated_start", "annotated_end"});
        },
        "View on exon information as xarray Dataset")
    .def("__repr__", [](const Exons& self) -> std::string {
        std::ostringstream oss;
        oss << "Exons<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &Exons::size);
}

void init_GeneJunctions(
    py::class_<majiq::GeneJunctions, std::shared_ptr<majiq::GeneJunctions>>& pyGeneJunctions) {
  using majiq::GeneJunctions;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyGeneJunctions
    .def_static("from_netcdf",
        [](py::str x, std::shared_ptr<majiq::Genes> genes) {
        return ConnectionsFromNetcdf<GeneJunctions, JUNCTIONS_NC_GROUP>(
            genes, x);
        },
        "Load junctions from netcdf file",
        py::arg("netcdf_path"), py::arg("genes"))
    .def("to_netcdf",
        [](py::object& self, py::str out, py::str mode) {
        using namespace py::literals;
        self.attr("df")().attr("to_netcdf")(
            out, "mode"_a = mode, "group"_a = py::str(JUNCTIONS_NC_GROUP));
        return;
        },
        "Save junctions to netcdf file",
        py::arg("netcdf_path"), py::arg("mode"))
    .def_property_readonly("gene_idx",
        [](py::object& junctions_obj) -> py::array_t<size_t> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetof(majiq::GeneJunction, gene.gene_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[int] of indexes indicating gene junction belongs to")
    .def_property_readonly("start",
        [](py::object& junctions_obj) -> py::array_t<position_t> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetof(majiq::GeneJunction, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[int] of junction starts")
    .def_property_readonly("end",
        [](py::object& junctions_obj) -> py::array_t<position_t> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetof(majiq::GeneJunction, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[int] of junction ends")
    .def_property_readonly("denovo",
        [](py::object& junctions_obj) -> py::array_t<bool> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetof(majiq::GeneJunction, data.denovo);
        return ArrayFromVectorAndOffset<bool, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[bool] indicating if connection was not found in annotations")
    .def_property_readonly("passed_build",
        [](py::object& junctions_obj) -> py::array_t<bool> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetof(majiq::GeneJunction, data.passed_build);
        return ArrayFromVectorAndOffset<bool, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[bool] indicating if passed build criteria to be in LSV")
    .def_property_readonly("simplified",
        [](py::object& junctions_obj) -> py::array_t<bool> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetof(majiq::GeneJunction, data.simplified);
        return ArrayFromVectorAndOffset<bool, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[bool] indicating if the connection is simplified")
    .def("df",
        [](py::object& junctions) -> py::object {
        return majiq_pybind::XarrayDatasetFromObject(junctions, "junction_idx",
            {"gene_idx", "start", "end",
            "denovo", "passed_build", "simplified"});
        },
        "View on junction information as xarray Dataset")
    .def("__repr__", [](const GeneJunctions& self) -> std::string {
        std::ostringstream oss;
        oss << "GeneJunctions<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &GeneJunctions::size);
}

void init_Introns(
    py::class_<majiq::Introns, std::shared_ptr<majiq::Introns>>& pyIntrons) {
  using majiq::Introns;
  using majiq::position_t;
  using majiq_pybind::ArrayFromVectorAndOffset;
  pyIntrons
    .def_static("from_netcdf",
        [](py::str x, std::shared_ptr<majiq::Genes> genes) {
        return ConnectionsFromNetcdf<Introns, INTRONS_NC_GROUP>(genes, x);
        },
        "Load introns from netcdf file",
        py::arg("netcdf_path"), py::arg("genes"))
    .def("to_netcdf",
        [](py::object& self, py::str out, py::str mode) {
        using namespace py::literals;
        self.attr("df")().attr("to_netcdf")(
            out, "mode"_a = mode, "group"_a = py::str(INTRONS_NC_GROUP));
        return;
        },
        "Save introns to netcdf file",
        py::arg("netcdf_path"), py::arg("mode"))
    .def_property_readonly("gene_idx",
        [](py::object& introns_obj) -> py::array_t<size_t> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetof(majiq::Intron, gene.gene_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[int] of indexes indicating gene intron belongs to")
    .def_property_readonly("start",
        [](py::object& introns_obj) -> py::array_t<position_t> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetof(majiq::Intron, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[int] of intron starts")
    .def_property_readonly("end",
        [](py::object& introns_obj) -> py::array_t<position_t> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetof(majiq::Intron, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[int] of intron ends")
    .def_property_readonly("denovo",
        [](py::object& introns_obj) -> py::array_t<bool> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetof(majiq::Intron, data.denovo);
        return ArrayFromVectorAndOffset<bool, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if connection was not found in annotations")
    .def_property_readonly("passed_build",
        [](py::object& introns_obj) -> py::array_t<bool> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetof(majiq::Intron, data.passed_build);
        return ArrayFromVectorAndOffset<bool, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if passed build criteria to be in LSV")
    .def_property_readonly("simplified",
        [](py::object& introns_obj) -> py::array_t<bool> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetof(majiq::Intron, data.simplified);
        return ArrayFromVectorAndOffset<bool, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if the connection is simplified")
    .def("df",
        [](py::object& introns) -> py::object {
        return majiq_pybind::XarrayDatasetFromObject(introns, "intron_idx",
            {"gene_idx", "start", "end",
            "denovo", "passed_build", "simplified"});
        },
        "View on intron information as xarray Dataset")
    .def("__repr__", [](const Introns& self) -> std::string {
        std::ostringstream oss;
        oss << "Introns<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &Introns::size);
}

void init_SpliceGraph(py::class_<majiq::SpliceGraph>& pySpliceGraph) {
  using majiq::SpliceGraph;
  pySpliceGraph
    // constructors from netcdf, gff3
    .def_static("from_netcdf",
        [](py::str netcdf_path) {
        py::function load_dataset
          = py::module_::import("xarray").attr("load_dataset");
        auto contigs = ContigsFromNetcdf(netcdf_path);
        auto genes = GenesFromNetcdf(contigs, netcdf_path);
        auto exons = ExonsFromNetcdf(genes, netcdf_path);
        auto junctions
          = ConnectionsFromNetcdf<majiq::GeneJunctions, JUNCTIONS_NC_GROUP>(
              genes, netcdf_path);
        auto introns
          = ConnectionsFromNetcdf<majiq::Introns, INTRONS_NC_GROUP>(
              genes, netcdf_path);
        return SpliceGraph{contigs, genes, exons, junctions, introns};
        },
        "Load splicegraph from netCDF file",
        py::arg("netcdf_path"))
    .def_static("from_gff3",
        [](std::string gff3_path, bool process_ir) {
          using majiq::gff3::SpliceGraphBuilder;
          SpliceGraphBuilder builder{};
          return builder.from_gff3(gff3_path, process_ir);
        },
        "Create splicegraph from input GFF3 file",
        py::arg("gff3_path"), py::arg("process_ir"))
    // save to file
    .def("to_netcdf",
        [](py::object& sg, py::str output_path) {
        // don't write to existing file
        py::object Path = py::module_::import("pathlib").attr("Path");
        if (Path(output_path).attr("exists")().cast<bool>()) {
          std::ostringstream oss;
          oss << "Cannot save result to already existing file "
              << output_path.cast<std::string>();
          throw std::invalid_argument(oss.str());
        }
        sg.attr("_exons").attr("to_netcdf")(output_path, "w");
        sg.attr("_introns").attr("to_netcdf")(output_path, "a");
        sg.attr("_junctions").attr("to_netcdf")(output_path, "a");
        sg.attr("_genes").attr("to_netcdf")(output_path, "a");
        sg.attr("_contigs").attr("to_netcdf")(output_path, "a");
        return;
        },
        R"pbdoc(
        Serialize splicegraph to netCDF file

        Parameters
        ----------
        output_path: str
            Path for resulting file. Raises error if file already exists.
        )pbdoc",
        py::arg("output_path"))
    // access underlying data
    .def_property_readonly("_exons", &SpliceGraph::exons,
        "Access the splicegraph's exons")
    .def_property_readonly("_introns", &SpliceGraph::introns,
        "Access the splicegraph's introns")
    .def_property_readonly("_junctions", &SpliceGraph::junctions,
        "Access the splicegraph's junctions")
    .def_property_readonly("_genes", &SpliceGraph::genes,
        "Access the splicegraph's genes")
    .def_property_readonly("_overgenes", &SpliceGraph::overgenes,
        "Access the spicegraph's overgenes")
    .def_property_readonly("_contigs", &SpliceGraph::contigs,
        "Access the splicegraph's contigs")
    // access underlying data as xarray datasets
    .def_property_readonly("exons",
        [](py::object& sg) { return sg.attr("_exons").attr("df")(); },
        "xr.Dataset view of splicegraph's exons")
    .def_property_readonly("introns",
        [](py::object& sg) { return sg.attr("_introns").attr("df")(); },
        "xr.Dataset view of splicegraph's introns")
    .def_property_readonly("junctions",
        [](py::object& sg) { return sg.attr("_junctions").attr("df")(); },
        "xr.Dataset view of splicegraph's junctions")
    .def_property_readonly("genes",
        [](py::object& sg) { return sg.attr("_genes").attr("df")(); },
        "xr.Dataset view of splicegraph's genes")
    .def_property_readonly("overgenes",
        [](py::object& sg) { return sg.attr("_overgenes").attr("df")(); },
        "xr.Dataset view of splicegraph's overgenes")
    .def_property_readonly("contigs",
        [](py::object& sg) { return sg.attr("_contigs").attr("df")(); },
        "xr.Dataset view of splicegraph's contigs")
    // computation
    .def("potential_introns", &SpliceGraph::potential_introns,
        "potential introns between splicegraph exons")
    // string representation of splicegraph
    .def("__repr__", [](const SpliceGraph& sg) -> std::string {
        std::ostringstream oss;
        oss << sg;
        return oss.str();
        });
}

void init_SpliceGraphAll(py::module_& m) {
  using majiq::Contigs;
  using majiq::Genes;
  using majiq::Exons;
  using majiq::Introns;
  using majiq::GeneJunctions;
  using majiq::SpliceGraph;
  auto pyContigs = py::class_<Contigs, std::shared_ptr<Contigs>>(m, "Contigs",
      "Splicegraph contigs");
  auto pyGenes = py::class_<Genes, std::shared_ptr<Genes>>(m, "Genes",
      "Splicegraph genes");
  auto pyExons = py::class_<Exons, std::shared_ptr<Exons>>(m, "Exons",
      "Splicegraph exons");
  auto pyIntrons = py::class_<Introns, std::shared_ptr<Introns>>(m, "Introns",
      "Splicegraph introns");
  auto pyGeneJunctions = py::class_<
      GeneJunctions, std::shared_ptr<GeneJunctions>
    >(m, "GeneJunctions", "Splicegraph junctions");
  auto pySpliceGraph = py::class_<SpliceGraph>(m, "SpliceGraph",
      "Splicegraph managing exons, junctions, and introns within genes");

  init_Contigs(pyContigs);
  init_Genes(pyGenes);
  init_Exons(pyExons);
  init_GeneJunctions(pyGeneJunctions);
  init_Introns(pyIntrons);
  init_SpliceGraph(pySpliceGraph);
}
