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

#include "MajiqTypes.hpp"
#include "SpliceGraph.hpp"
#include "GFF3.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;  // bring in _a literal


// <https://stackoverflow.com/questions/13180842/how-to-calculate-offset-of-a-class-member-at-compile-time>
template <typename T, typename U>
constexpr size_t offsetOf(U T::*member) {
  return reinterpret_cast<char*>(&(static_cast<T*>(nullptr)->*member))
    - static_cast<char*>(nullptr);
}
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


PYBIND11_MODULE(new_majiq, m) {
  // documentation of module
  m.doc() = R"pbdoc(
  MAJIQ classes for splicegraph

  .. current_module::new_majiq
  .. autosummary::
    :toctree: _generate

    ClosedInterval
  )pbdoc";

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
  auto pyContigs = py::class_<Contigs, std::shared_ptr<Contigs>>(m, "Contigs",
      "Splicegraph contigs");

  pyExons
    .def_property_readonly("gene_idx",
        [](py::object& exons_obj) -> py::array_t<size_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetOf(&majiq::Exon::gene)
            + offsetOf(&majiq::KnownGene::gene_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of indexes indicating gene exon belongs to")
    .def_property_readonly("start",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetOf(&majiq::Exon::coordinates)
            + offsetOf(&majiq::ClosedInterval::start);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of exon starts")
    .def_property_readonly("end",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetOf(&majiq::Exon::coordinates)
            + offsetOf(&majiq::ClosedInterval::end);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of exon ends")
    .def_property_readonly("annotated_start",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetOf(&majiq::Exon::annotated_coordinates)
            + offsetOf(&majiq::ClosedInterval::start);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of annotated exon starts")
    .def_property_readonly("annotated_end",
        [](py::object& exons_obj) -> py::array_t<position_t> {
        Exons& exons = exons_obj.cast<Exons&>();
        const size_t offset = offsetOf(&majiq::Exon::annotated_coordinates)
            + offsetOf(&majiq::ClosedInterval::end);
        return ArrayFromVectorAndOffset<position_t, majiq::Exon>(
            exons.data(), offset, exons_obj);
        },
        "array[int] of annotated exon ends")
    .def("df",
        [](py::object& exons) -> py::object {
        py::module_ np = py::module_::import("numpy");
        py::module_ xr = py::module_::import("xarray");
        py::function xr_DataArray = xr.attr("DataArray");
        py::function xr_Dataset = xr.attr("Dataset");
        py::function np_arange = np.attr("arange");
        // get coordinates for exon
        auto index = py::tuple(2);
        index[0] = py::str("exon_idx");
        index[1] = np_arange(exons.attr("__len__")());
        // list of coordinates for DataArray
        auto index_list = py::list();
        index_list.append(index);
        // so the values we want out are:
        auto coordinates = py::dict(
            "gene_idx"_a = xr_DataArray(exons.attr("gene_idx"), index_list),
            "start"_a = xr_DataArray(exons.attr("start"), index_list),
            "end"_a = xr_DataArray(exons.attr("end"), index_list),
            "annotated_start"_a = xr_DataArray(exons.attr("annotated_start"), index_list),
            "annotated_end"_a = xr_DataArray(exons.attr("annotated_end"), index_list));
        return xr_Dataset("coords"_a = coordinates);
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
        const size_t offset = offsetOf(&majiq::Intron::gene)
            + offsetOf(&majiq::KnownGene::gene_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[int] of indexes indicating gene intron belongs to")
    .def_property_readonly("start",
        [](py::object& introns_obj) -> py::array_t<position_t> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetOf(&majiq::Intron::coordinates)
            + offsetOf(&majiq::ClosedInterval::start);
        return ArrayFromVectorAndOffset<position_t, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[int] of intron starts")
    .def_property_readonly("end",
        [](py::object& introns_obj) -> py::array_t<position_t> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetOf(&majiq::Intron::coordinates)
            + offsetOf(&majiq::ClosedInterval::end);
        return ArrayFromVectorAndOffset<position_t, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[int] of intron ends")
    .def_property_readonly("denovo",
        [](py::object& introns_obj) -> py::array_t<bool> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetOf(&majiq::Intron::denovo);
        return ArrayFromVectorAndOffset<bool, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if connection was not found in annotations")
    .def_property_readonly("passed_build",
        [](py::object& introns_obj) -> py::array_t<bool> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetOf(&majiq::Intron::passed_build);
        return ArrayFromVectorAndOffset<bool, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if passed build criteria to be in LSV")
    .def_property_readonly("simplified",
        [](py::object& introns_obj) -> py::array_t<bool> {
        Introns& introns = introns_obj.cast<Introns&>();
        const size_t offset = offsetOf(&majiq::Intron::simplified);
        return ArrayFromVectorAndOffset<bool, majiq::Intron>(
            introns.data(), offset, introns_obj);
        },
        "array[bool] indicating if the connection is simplified")
    .def("df",
        [](py::object& introns) -> py::object {
        py::module_ np = py::module_::import("numpy");
        py::module_ xr = py::module_::import("xarray");
        py::function xr_DataArray = xr.attr("DataArray");
        py::function xr_Dataset = xr.attr("Dataset");
        py::function np_arange = np.attr("arange");
        // get coordinates for intron
        auto index = py::tuple(2);
        index[0] = py::str("intron_idx");
        index[1] = np_arange(introns.attr("__len__")());
        // list of coordinates for DataArray
        auto index_list = py::list();
        index_list.append(index);
        // so the values we want out are:
        auto coordinates = py::dict(
            "gene_idx"_a = xr_DataArray(introns.attr("gene_idx"), index_list),
            "start"_a = xr_DataArray(introns.attr("start"), index_list),
            "end"_a = xr_DataArray(introns.attr("end"), index_list),
            "denovo"_a = xr_DataArray(introns.attr("denovo"), index_list),
            "passed_build"_a = xr_DataArray(introns.attr("passed_build"), index_list),
            "simplified"_a = xr_DataArray(introns.attr("simplified"), index_list));
        return xr_Dataset("coords"_a = coordinates);
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
        const size_t offset = offsetOf(&majiq::GeneJunction::gene)
            + offsetOf(&majiq::KnownGene::gene_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[int] of indexes indicating gene junction belongs to")
    .def_property_readonly("start",
        [](py::object& junctions_obj) -> py::array_t<position_t> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetOf(&majiq::GeneJunction::coordinates)
            + offsetOf(&majiq::OpenInterval::start);
        return ArrayFromVectorAndOffset<position_t, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[int] of junction starts")
    .def_property_readonly("end",
        [](py::object& junctions_obj) -> py::array_t<position_t> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetOf(&majiq::GeneJunction::coordinates)
            + offsetOf(&majiq::OpenInterval::end);
        return ArrayFromVectorAndOffset<position_t, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[int] of junction ends")
    .def_property_readonly("denovo",
        [](py::object& junctions_obj) -> py::array_t<bool> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetOf(&majiq::GeneJunction::denovo);
        return ArrayFromVectorAndOffset<bool, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[bool] indicating if connection was not found in annotations")
    .def_property_readonly("passed_build",
        [](py::object& junctions_obj) -> py::array_t<bool> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetOf(&majiq::GeneJunction::passed_build);
        return ArrayFromVectorAndOffset<bool, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[bool] indicating if passed build criteria to be in LSV")
    .def_property_readonly("simplified",
        [](py::object& junctions_obj) -> py::array_t<bool> {
        GeneJunctions& junctions = junctions_obj.cast<GeneJunctions&>();
        const size_t offset = offsetOf(&majiq::GeneJunction::simplified);
        return ArrayFromVectorAndOffset<bool, majiq::GeneJunction>(
            junctions.data(), offset, junctions_obj);
        },
        "array[bool] indicating if the connection is simplified")
    .def("df",
        [](py::object& junctions) -> py::object {
        py::module_ np = py::module_::import("numpy");
        py::module_ xr = py::module_::import("xarray");
        py::function xr_DataArray = xr.attr("DataArray");
        py::function xr_Dataset = xr.attr("Dataset");
        py::function np_arange = np.attr("arange");
        // get coordinates for junction
        auto index = py::tuple(2);
        index[0] = py::str("junction_idx");
        index[1] = np_arange(junctions.attr("__len__")());
        // list of coordinates for DataArray
        auto index_list = py::list();
        index_list.append(index);
        // so the values we want out are:
        auto coordinates = py::dict(
            "gene_idx"_a = xr_DataArray(junctions.attr("gene_idx"), index_list),
            "start"_a = xr_DataArray(junctions.attr("start"), index_list),
            "end"_a = xr_DataArray(junctions.attr("end"), index_list),
            "denovo"_a = xr_DataArray(junctions.attr("denovo"), index_list),
            "passed_build"_a = xr_DataArray(junctions.attr("passed_build"), index_list),
            "simplified"_a = xr_DataArray(junctions.attr("simplified"), index_list));
        return xr_Dataset("coords"_a = coordinates);
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
        const size_t offset = offsetOf(&majiq::Gene::strand);
        return ArrayFromVectorAndOffset<std::array<char, 1>, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[char] of characters indicating strand of each gene")
    .def_property_readonly("contig_idx",
        [](py::object& genes_obj) -> py::array_t<size_t> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetOf(&majiq::Gene::contig)
            + offsetOf(&majiq::KnownContig::contig_idx);
        return ArrayFromVectorAndOffset<size_t, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[int] of indexes indicating contig gene belongs to")
    .def_property_readonly("start",
        [](py::object& genes_obj) -> py::array_t<position_t> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetOf(&majiq::Gene::interval)
            + offsetOf(&majiq::ClosedInterval::start);
        return ArrayFromVectorAndOffset<position_t, majiq::Gene>(
            genes.data(), offset, genes_obj);
        },
        "array[int] of gene starts matching gene_idx")
    .def_property_readonly("end",
        [](py::object& genes_obj) -> py::array_t<position_t> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetOf(&majiq::Gene::interval)
            + offsetOf(&majiq::ClosedInterval::end);
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
        py::module_ np = py::module_::import("numpy");
        py::module_ xr = py::module_::import("xarray");
        py::function xr_DataArray = xr.attr("DataArray");
        py::function xr_Dataset = xr.attr("Dataset");
        py::function np_arange = np.attr("arange");
        // get coordinates for gene
        auto index = py::tuple(2);
        index[0] = py::str("gene_idx");
        index[1] = np_arange(genes.attr("__len__")());
        // list of coordinates for DataArray
        auto index_list = py::list();
        index_list.append(index);
        // so the values we want out are:
        auto coordinates = py::dict(
            "contig_idx"_a = xr_DataArray(genes.attr("contig_idx"), index_list),
            "start"_a = xr_DataArray(genes.attr("start"), index_list),
            "end"_a = xr_DataArray(genes.attr("end"), index_list),
            "strand"_a = xr_DataArray(genes.attr("strand"), index_list),
            "gene_id"_a = xr_DataArray(genes.attr("gene_id"), index_list),
            "gene_name"_a = xr_DataArray(genes.attr("gene_name"), index_list));
        return xr_Dataset("coords"_a = coordinates);
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


  pyContigs
    .def_property_readonly("seqid", &Contigs::seqids,
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("df",
        [](py::object& contigs) -> py::object {
        py::module_ np = py::module_::import("numpy");
        py::module_ xr = py::module_::import("xarray");
        py::function xr_DataArray = xr.attr("DataArray");
        py::function xr_Dataset = xr.attr("Dataset");
        py::function np_arange = np.attr("arange");
        // get coordinates for contig
        auto index = py::tuple(2);
        index[0] = py::str("contig_idx");
        index[1] = np_arange(contigs.attr("__len__")());
        // list of coordinates for DataArray
        auto index_list = py::list();
        index_list.append(index);
        // so the values we want out are:
        auto coordinates = py::dict(
            "seqid"_a = xr_DataArray(contigs.attr("seqid"), index_list));
        return xr_Dataset("coords"_a = coordinates);
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
    .def_property_readonly("contigs",
        [](py::object& sg) { return sg.attr("_contigs").attr("df")(); },
        "xr.Dataset view of splicegraph's contigs")
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
