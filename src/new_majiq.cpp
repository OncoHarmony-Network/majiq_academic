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

#include "MajiqTypes.hpp"
#include "SpliceGraph.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;


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

  pyGenes
    .def_property_readonly("strand",
        [](py::object& genes_obj) -> py::array_t<char> {
        Genes& genes = genes_obj.cast<Genes&>();
        const size_t offset = offsetOf(&majiq::Gene::strand);
        return ArrayFromVectorAndOffset<char, majiq::Gene>(
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
    .def("__repr__", [](const Contigs& self) -> std::string {
        std::ostringstream oss;
        oss << self;
        return oss.str();
        })
    .def("__len__", &Contigs::size)
    .def("__contains__",
        [](const Contigs& s, seqid_t x) -> bool { return s.contains(x); });

  pySpliceGraph
    // empty constructor
    .def(py::init<>(), "Create empty SpliceGraph")
    // TODO(jaicher) Make more useful constructors
    // string representation of splicegraph
    .def("__repr__", [](const SpliceGraph& sg) -> std::string {
        std::ostringstream oss;
        oss << sg;
        return oss.str();
        })
    // access underlying data
    .def_property_readonly("exons", &SpliceGraph::exons,
        "Access the splicegraph's exons")
    .def_property_readonly("introns", &SpliceGraph::introns,
        "Access the splicegraph's introns")
    .def_property_readonly("junctions", &SpliceGraph::junctions,
        "Access the splicegraph's junctions")
    .def_property_readonly("genes", &SpliceGraph::genes,
        "Access the splicegraph's genes")
    .def_property_readonly("contigs", &SpliceGraph::contigs,
        "Access the splicegraph's contigs")
    // add elements one at a time
    // (we probably won't be using these outside of early debugging)
    .def("add_gene", &SpliceGraph::AddGene,
        R"pbdoc(
        Register the specified gene as valid in the splicegraph

        Parameters
        ----------
        gene_id: str
            Unique ID of specified gene (e.g. ENSG0000...)
        seqid: str
            Identifier of contig/chromosome gene is on
        start, end: int
            Coordinates from annotation for gene
        strand_forward: bool
            True if forward strand ('+'), False if reverse strand ('-')
        gene_name: str
            Alternative symbol for gene (e.g. HGNC ID)

        Notes
        -----
        Calls with previously added gene_id will be ignored.
        )pbdoc",
        py::arg("gene_id"), py::arg("seqid"), py::arg("start"), py::arg("end"),
        py::arg("strand_forward"), py::arg("gene_name"))
    .def("add_exon", &SpliceGraph::AddExon,
        "Add exon to splicegraph. Unregistered gene_id will raise error.",
        py::arg("gene_id"), py::arg("start"), py::arg("end"))
    .def("add_junction", &SpliceGraph::AddJunction,
        "Add junction to splicegraph. Unregistered gene_id will raise error.",
        py::arg("gene_id"), py::arg("start"), py::arg("end"))
    .def("add_intron", &SpliceGraph::AddIntron,
        "Add intron to splicegraph. Unregistered gene_id will raise error.",
        py::arg("gene_id"), py::arg("start"), py::arg("end"));


#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}  // PYBIND11_MODULE(new_majiq, m)
