/**
 * new_majiq.cpp
 *
 * Python bindings to new MAJIQ code
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include <string>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <array>
#include <initializer_list>
#include <cstddef>

#include "MajiqTypes.hpp"
#include "SpliceGraph.hpp"
#include "GFF3.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;  // bring in _a literal

/*
 * Create read-only array view into vector with offset (i.e. for struct member)
 */
template <class OutputT, class InputT>
py::array_t<OutputT> ArrayFromVectorAndOffset(
    const std::vector<InputT>& src,
    size_t offset,
    py::object& handle) {
  // pointer to first element of src after adding specified memory offset
  const OutputT* first = reinterpret_cast<const OutputT*>(
      reinterpret_cast<const char*>(src.data()) + offset);
  // construct array
  py::array_t<OutputT> result = py::array_t(
      // shape
      {src.size()},
      // strides
      {sizeof(InputT)},
      // pointer to first element
      first,
      // object to handle array
      handle);
  // set array to readonly -- discourage Python from editing it
  py::detail::array_proxy(result.ptr())->flags
    &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  // return resulting array
  return result;
}
/*
 * Create read-only array view into vector with offset (i.e. for struct member)
 */
template <class OutputT>
py::array_t<OutputT> ArrayFromOffsetsVector(
    const std::vector<OutputT>& src,
    bool is_start,
    py::object& handle) {
  // pointer to first element of src after adding offset for start or end
  const OutputT* first = src.data() + (is_start ? 0 : 1);
  // construct array
  py::array_t<OutputT> result = py::array_t(
      // shape
      {src.size() - 1},
      // strides
      {sizeof(OutputT)},
      // pointer to first element
      first,
      // object to handle array
      handle);
  // set array to readonly -- discourage Python from editing it
  py::detail::array_proxy(result.ptr())->flags
    &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  // return resulting array
  return result;
}

/**
 * Extract attrs from features to build xarray Dataset
 */
py::object XarrayDatasetFromObject(
    py::object features,
    py::str idx_name,
    std::initializer_list<py::str> attrs) {
  py::function np_arange = py::module_::import("numpy").attr("arange");
  py::module_ xr = py::module_::import("xarray");
  py::function xr_Dataset = xr.attr("Dataset");
  py::function xr_DataArray = xr.attr("DataArray");
  // coordinates for dataset in dictionary
  py::dict coordinates;
  // define the index dimension
  coordinates[idx_name] = xr_DataArray(
      np_arange(features.attr("__len__")()), "dims"_a = idx_name);
  // extract remaining attributes
  for (auto x : attrs) {
    coordinates[x] = xr_DataArray(features.attr(x), "dims"_a = idx_name);
  }
  return xr_Dataset("coords"_a = coordinates);
}

using majiq::position_t;
using majiq::seqid_t;
using majiq::geneid_t;
using majiq::genename_t;

using majiq::SpliceGraph;
using majiq::Exons;
using majiq::Introns;
using majiq::GeneJunctions;
using majiq::Genes;
using majiq::Contigs;
using majiq::OverGenes;

std::shared_ptr<Contigs> ContigsFromDataset(py::object xr_contigs) {
  auto result = std::make_shared<Contigs>();
  py::list seqids = xr_contigs
    .attr("__getitem__")("seqid")
    .attr("values").attr("tolist")();
  for (auto seqid : seqids) {
    result->add(seqid.cast<std::string>());
  }
  return result;
}
std::shared_ptr<Genes> GenesFromDataset(
    std::shared_ptr<Contigs> contigs, py::object xr_genes) {
  auto get_array = [&xr_genes](py::str key) {
    py::function np_array = py::module_::import("numpy").attr("array");
    return np_array(xr_genes.attr("__getitem__")(key));
  };
  // extract usable values from dataset
  py::array_t<size_t> arr_contig_idx = get_array("contig_idx");
  py::array_t<position_t> arr_start = get_array("start");
  py::array_t<position_t> arr_end = get_array("end");
  py::array_t<std::array<char, 1>> arr_strand = get_array("strand");
  py::list geneid
    = xr_genes.attr("__getitem__")("gene_id").attr("values").attr("tolist")();
  py::list genename
    = xr_genes.attr("__getitem__")("gene_name").attr("values").attr("tolist")();
  // unchecked accesses to numpy array
  auto contig_idx = arr_contig_idx.unchecked<1>();
  auto start = arr_start.unchecked<1>();
  auto end = arr_end.unchecked<1>();
  auto strand = arr_strand.unchecked<1>();
  // create vector of genes matching input arrays
  using majiq::Gene;
  using majiq::KnownContig;
  using majiq::ClosedInterval;
  using majiq::GeneStrandness;

  std::vector<Gene> gene_vec{};
  gene_vec.reserve(geneid.size());
  for (size_t i = 0; i < geneid.size(); ++i) {
    gene_vec.push_back(Gene{
        KnownContig{contig_idx(i), contigs}, ClosedInterval{start(i), end(i)},
        static_cast<GeneStrandness>(strand(i)[0]),
        geneid[i].cast<geneid_t>(), genename[i].cast<genename_t>()});
  }
  return std::make_shared<Genes>(std::move(gene_vec));
}
std::shared_ptr<Exons> ExonsFromDataset(
    std::shared_ptr<Genes> genes, py::object xr_exons) {
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
  using majiq::Exon;
  using majiq::KnownGene;
  using majiq::ClosedInterval;

  std::vector<Exon> exon_vec{};
  exon_vec.reserve(gene_idx.shape(0));
  for (py::ssize_t i = 0; i < gene_idx.shape(0); ++i) {
    exon_vec.push_back(Exon{
        KnownGene{gene_idx(i), genes}, ClosedInterval{start(i), end(i)},
        ClosedInterval{ann_start(i), ann_end(i)}});
  }
  return std::make_shared<Exons>(std::move(exon_vec));
}
template <class Connections>
std::shared_ptr<Connections> ConnectionsFromDataset(
    std::shared_ptr<Genes> genes, py::object xr_connections) {
  using IntervalT = typename Connections::IntervalT;
  using RegionT = typename Connections::value_type;
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
        KnownGene{gene_idx(i), genes}, IntervalT{start(i), end(i)},
        denovo(i), passed_build(i), simplified(i)});
  }
  return std::make_shared<Connections>(std::move(connection_vec));
}

PYBIND11_MODULE(new_majiq, m) {
  // documentation of module
  m.doc() = R"pbdoc(
  MAJIQ classes for splicegraph

  .. current_module::new_majiq
  .. autosummary::
    :toctree: _generate

    ClosedInterval
  )pbdoc";


  // vector<string> will be bound to Python list
  py::bind_vector<std::vector<std::string>>(m, "VectorString");

  // allow nonconst operations
  auto pySpliceGraph = py::class_<SpliceGraph>(m, "SpliceGraph",
      "Splicegraph managing exons, junctions, and introns within genes");

  // should only be used to provide views into splicegraph
  // and for de/serialization
  // NOTE: pybind11 casts away const-ness, so have to be careful
  // <https://pybind11.readthedocs.io/en/stable/limitations.html?highlight=const#design-choices>
  auto pyExons = py::class_<Exons, std::shared_ptr<Exons>>(m, "Exons",
      "Splicegraph exons");
  auto pyIntrons = py::class_<Introns, std::shared_ptr<Introns>>(m, "Introns",
      "Splicegraph introns");
  auto pyGeneJunctions = py::class_<
      GeneJunctions, std::shared_ptr<GeneJunctions>
    >(m, "GeneJunctions", "Splicegraph junctions");
  auto pyGenes = py::class_<Genes, std::shared_ptr<Genes>>(m, "Genes",
      "Splicegraph genes");
  auto pyOverGenes = py::class_<OverGenes, std::shared_ptr<OverGenes>>(
      m, "OverGenes", "Splicegraph groups of overlapping genes");
  auto pyContigs = py::class_<Contigs, std::shared_ptr<Contigs>>(m, "Contigs",
      "Splicegraph contigs");

  pyExons
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
        return XarrayDatasetFromObject(exons, "exon_idx",
            {"gene_idx", "start", "end", "annotated_start", "annotated_end"});
        },
        "View on exon information as xarray Dataset")
    .def("__repr__", [](const Exons& self) -> std::string {
        std::ostringstream oss;
        oss << "Exons<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &Exons::size);

  pyIntrons
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
        return XarrayDatasetFromObject(introns, "intron_idx",
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

  pyGeneJunctions
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
        return XarrayDatasetFromObject(junctions, "junction_idx",
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

  pyGenes
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
    .def_property_readonly("gene_id", &Genes::geneids,
        "Sequence[str] of gene ids in order matching gene_idx")
    .def_property_readonly("gene_name", &Genes::genenames,
        "Sequence[str] of gene names in order matching gene_idx")
    .def("df",
        [](py::object& genes) -> py::object {
        return XarrayDatasetFromObject(genes, "gene_idx",
            {"contig_idx", "start", "end", "strand", "gene_id", "gene_name"});
        },
        "View on gene information as xarray Dataset")
    .def("__repr__", [](const Genes& self) -> std::string {
        std::ostringstream oss;
        oss << "Genes<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &Genes::size)
    .def("__contains__",
        [](const Genes& self, geneid_t x) { return self.contains(x); })
    .def("__getitem__",
        [](const Genes& self, geneid_t x) { return self.get_gene_idx(x); },
        "gene_idx for specified gene_id", py::arg("gene_id"));

  pyOverGenes
    .def_property_readonly("contig_idx",
        [](py::object& overgenes_obj) -> py::array_t<size_t> {
        OverGenes& overgenes = overgenes_obj.cast<OverGenes&>();
        const size_t offset = offsetof(majiq::OverGene, contig.contig_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::OverGene>(
            overgenes.overgenes(), offset, overgenes_obj);
        },
        "array[int] of indexes indicating contig overgene belongs to")
    .def_property_readonly("start",
        [](py::object& overgenes_obj) -> py::array_t<position_t> {
        OverGenes& overgenes = overgenes_obj.cast<OverGenes&>();
        const size_t offset = offsetof(majiq::OverGene, coordinates.start);
        return ArrayFromVectorAndOffset<position_t, majiq::OverGene>(
            overgenes.overgenes(), offset, overgenes_obj);
        },
        "array[int] of overgene starts")
    .def_property_readonly("end",
        [](py::object& overgenes_obj) -> py::array_t<position_t> {
        OverGenes& overgenes = overgenes_obj.cast<OverGenes&>();
        const size_t offset = offsetof(majiq::OverGene, coordinates.end);
        return ArrayFromVectorAndOffset<position_t, majiq::OverGene>(
            overgenes.overgenes(), offset, overgenes_obj);
        },
        "array[int] of overgene ends")
    .def_property_readonly("gene_idx_begin",
        [](py::object& overgenes_obj) -> py::array_t<size_t> {
        OverGenes& overgenes = overgenes_obj.cast<OverGenes&>();
        return ArrayFromOffsetsVector(
            overgenes.overgene_offsets(), true, overgenes_obj);
        },
        "array[int] of gene_idx for first gene in the overgene")
    .def_property_readonly("gene_idx_end",
        [](py::object& overgenes_obj) -> py::array_t<size_t> {
        OverGenes& overgenes = overgenes_obj.cast<OverGenes&>();
        return ArrayFromOffsetsVector(
            overgenes.overgene_offsets(), false, overgenes_obj);
        },
        "array[int] of gene_idx following that of last gene in the overgene")
    .def("df",
        [](py::object& overgenes) -> py::object {
        return XarrayDatasetFromObject(overgenes, "overgene_idx",
            {"contig_idx", "start", "end", "gene_idx_begin", "gene_idx_end"});
        },
        "View on overgene information as xarray Dataset")
    .def("__repr__", [](const OverGenes& self) -> std::string {
        std::ostringstream oss;
        oss << "OverGenes<" << self.size() << " total>";
        return oss.str();
        })
    .def("__len__", &OverGenes::size);


  pyContigs
    .def_property_readonly("seqid", &Contigs::seqids,
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("df",
        [](py::object& contigs) -> py::object {
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

  pySpliceGraph
    // constructor from GFF3
    .def(py::init(
          [](std::string gff3_path, bool process_ir) {
            using majiq::gff3::SpliceGraphBuilder;
            SpliceGraphBuilder builder{};
            return builder.from_gff3(gff3_path, process_ir);
          }),
        "Create splicegraph from input GFF3 file",
        py::arg("gff3_path"), py::arg("process_ir"))
    // TODO(jaicher) Make more useful constructors
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
        auto save = [&sg, &output_path](py::str key, py::str mode) {
          sg.attr(key).attr("to_netcdf")(
              output_path, "mode"_a = mode, "group"_a = key);
        };
        save("exons", "w");
        save("introns", "a");
        save("junctions", "a");
        save("genes", "a");
        save("contigs", "a");
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
    .def(py::init([](py::str netcdf_path) {
          py::function load_dataset
            = py::module_::import("xarray").attr("load_dataset");
          auto contigs = ContigsFromDataset(
              load_dataset(netcdf_path, "group"_a = "contigs"));
          auto genes = GenesFromDataset(
              contigs, load_dataset(netcdf_path, "group"_a = "genes"));
          auto exons = ExonsFromDataset(
              genes, load_dataset(netcdf_path, "group"_a = "exons"));
          auto junctions = ConnectionsFromDataset<GeneJunctions>(
              genes, load_dataset(netcdf_path, "group"_a = "junctions"));
          auto introns = ConnectionsFromDataset<Introns>(
              genes, load_dataset(netcdf_path, "group"_a = "introns"));
          return SpliceGraph{
              std::move(contigs), std::move(genes), std::move(exons),
              std::move(junctions), std::move(introns)};
          }),
        "Load splicegraph from netCDF file", py::arg("netcdf_path"))
    // string representation of splicegraph
    .def("__repr__", [](const SpliceGraph& sg) -> std::string {
        std::ostringstream oss;
        oss << sg;
        return oss.str();
        });


#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}  // PYBIND11_MODULE(new_majiq, m)
