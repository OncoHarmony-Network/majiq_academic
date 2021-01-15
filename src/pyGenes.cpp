/**
 * pyGenes.cpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <sstream>

#include "pybind_utils.hpp"
#include "internals/Genes.hpp"

namespace py = pybind11;
using majiq::Genes;

void init_Genes(py::class_<Genes, std::shared_ptr<Genes>>& pyGenes) {
  using majiq_pybind::GENES_NC_GROUP;
  using majiq_pybind::ArrayFromVectorAndOffset;
  using majiq_pybind::XarrayDatasetFromObject;
  using majiq::Contigs;
  using majiq::position_t;
  using majiq::geneid_t;
  pyGenes
    .def_static("from_netcdf",
        [](py::str x, std::shared_ptr<Contigs> contigs) {
        return majiq_pybind::GenesFromNetcdf(contigs, x);
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
    .def_property_readonly("is_sorted", &majiq::Genes::is_sorted,
        "True if the genes sorted in contig/coordinate order")
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

}
