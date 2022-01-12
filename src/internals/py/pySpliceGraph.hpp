/**
 * pySpliceGraph.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSPLICEGRAPH_HPP
#define MAJIQ_PYBIND_PYSPLICEGRAPH_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <iostream>
#include <memory>
#include <string>

#include "../GFF3.hpp"
#include "../SpliceGraph.hpp"

#include "constants.hpp"

namespace majiq {
namespace bindings {

using pySpliceGraph_t = pybind11::class_<SpliceGraph>;

inline void init_GFF3Types(pybind11::module_& m) {
  auto pyGFFFeatureTypes = pybind11::enum_<majiq::gff3::FeatureType>(
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
  pybind11::bind_map<majiq::gff3::featuretype_map_t>(m, "GFF3Types");
}

inline void init_SpliceGraph(pySpliceGraph_t& pySpliceGraph) {
  pySpliceGraph
    // expose constructor from individual components
    .def(
        pybind11::init<
          const std::shared_ptr<majiq::Contigs>&,
          const std::shared_ptr<majiq::Genes>&,
          const std::shared_ptr<majiq::Exons>&,
          const std::shared_ptr<majiq::GeneJunctions>&,
          const std::shared_ptr<majiq::GeneIntrons>&
        >(),
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        Initialize splicegraph from components

        Initialize splicegraph from components. Typically will want to use the
        factory methods `from_gff3` to create all components
        )pbdoc",
        pybind11::arg("contigs"), pybind11::arg("genes"), pybind11::arg("exons"),
        pybind11::arg("junctions"), pybind11::arg("introns"))
    // constructors from gff3
    .def_static("from_gff3",
        [](std::string gff3_path, majiq::gff3::featuretype_map_t gff3_types,
            bool process_ir,
            pybind11::object log_function) {
          using majiq::gff3::GFF3ExonHierarchy;
          using majiq::gff3::GFF3TranscriptModels;
          using majiq::gff3::ToTranscriptModels;
          // load gff3 exon hierarchy, convert to MAJIQ gene/transcript/exons
          auto gff3_models
            = ToTranscriptModels(GFF3ExonHierarchy{gff3_path, gff3_types});
          {
            pybind11::gil_scoped_acquire gil;  // acquire GIL to use Python
            if (!log_function.is_none()) {
              for (const auto& [tx_type, tx_ct]
                  : gff3_models.skipped_transcript_type_ct_) {
                log_function(tx_type, "parent", tx_ct);
              }
              for (const auto& [gene_type, gene_ct]
                  : gff3_models.skipped_gene_type_ct_) {
                log_function(gene_type, "top-level ancestor", gene_ct);
              }
            }  // log skipped exons if log_function is not none
          }
          // convert to splicegraph
          return gff3_models.models_.ToSpliceGraph(process_ir);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        Create splicegraph from input GFF3 file

        Parameters
        ----------
        gff3_path: str
            Path to GFF3 file to process (supports gzipped files)
        gff3_types: GFF3Types
            Map from GFF3 type to feature type it should be interpreted in
            building transcript models of each gene
        process_ir: bool
            Should annotated retained introns be assessed
        log_function: Optional[Callable[[str, str, int]]]
            If not None, it will be called for each unaccepted parent
            (transcript) or top-level ancestor (gene) type that was neither
            explicitly accepted or ignored when building transcript models.
            First argument indicates the unaccepted GFF3 type,
            second argument indicates if it was unaccepted as a parent or
            top-level ancestor,
            third argument indicates how many unique cases the GFF3 type was
            missing in this way.
        )pbdoc"
        "Create splicegraph from input GFF3 file",
        pybind11::arg("gff3_path"),
        pybind11::arg("gff3_types"),
        pybind11::arg("process_ir") = DEFAULT_BUILD_PROCESS_IR,
        pybind11::arg("log_function") = pybind11::none())
    .def_static("infer_exons", &SpliceGraph::InferExons,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Infer exons from base annotated exons and junctions",
        pybind11::arg("base_exons"), pybind11::arg("junctions"))
    .def("make_group_junctions", &SpliceGraph::MakeGroupGenerator,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Create GroupJunctionsGenerator for the splicegraph junctions/exons")
    .def("make_build_junctions", &SpliceGraph::MakePassedGenerator,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
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
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "True if position close to following/preceding annotated exon in gene",
        pybind11::arg("gene_idx"),
        pybind11::arg("x"),
        pybind11::arg("to_following") = true)
    // access underlying data
    .def_property_readonly("_exons", &SpliceGraph::exons,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Access the splicegraph's exons")
    .def_property_readonly("_introns", &SpliceGraph::introns,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Access the splicegraph's introns")
    .def_property_readonly("_junctions", &SpliceGraph::junctions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Access the splicegraph's junctions")
    .def_property_readonly("_genes", &SpliceGraph::genes,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Access the splicegraph's genes")
    .def_property_readonly("_contigs", &SpliceGraph::contigs,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Access the splicegraph's contigs")
    .def_property_readonly("_exon_connections", &SpliceGraph::exon_connections,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Access the splicegraph's exon connections")
    // get sj introns on which coverage may be read
    .def("sj_introns", [](SpliceGraph& sg, bool stranded) {
        return majiq::SJIntrons::FromGeneExonsAndIntrons(
            *sg.exons(), *sg.introns(), stranded);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get contig introns (by strand or not) for splicegraph",
        pybind11::arg("stranded"))
    // that's really for debugging because it depends on experiment. Instead,
    // just create them while loading BAM
    .def("sj_introns_from_bam",
        [](SpliceGraph& sg, const char* infile, majiq::junction_pos_t num_bins,
          majiq::ExperimentStrandness exp_strandness, int nthreads) {
        return majiq::SJIntronsBins::FromBam(infile, num_bins, *sg.exons(),
            *sg.introns(), exp_strandness, nthreads);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
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
        pybind11::arg("bam_path"),
        pybind11::arg("num_bins"),
        pybind11::arg("experiment_strandness") = DEFAULT_BAM_STRANDNESS,
        pybind11::arg("nthreads") = DEFAULT_BAM_NTHREADS);
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSPLICEGRAPH_HPP
