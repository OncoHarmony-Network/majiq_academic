/**
 * module.cpp
 *
 * Python bindings for module with MAJIQ internals
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <string>
#include <vector>

#include "../MajiqTypes.hpp"

#include "globalRNG.hpp"
#include "pyConnections.hpp"
#include "pyContigs.hpp"
#include "pyEvents.hpp"
#include "pyEventsAlign.hpp"
#include "pyEventsCoverage.hpp"
#include "pyExonConnections.hpp"
#include "pyExons.hpp"
#include "pyExperimentThresholds.hpp"
#include "pyGenes.hpp"
#include "pyGeneJunctionsAccumulator.hpp"
#include "pyGeneModules.hpp"
#include "pyPassedIntrons.hpp"
#include "pyPassedJunctions.hpp"
#include "pySimplifierGroup.hpp"
#include "pySpliceGraph.hpp"
#include "pySpliceGraphValues.hpp"
#include "pySJBinsReads.hpp"
#include "pySJIntrons.hpp"
#include "pySJJunctions.hpp"


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


PYBIND11_MODULE(internals, m) {
  m.attr("__name__") = "new_majiq.internals";
#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
  m.doc() = R"pbdoc(
  MAJIQ internals C++ API

  .. current_module::new_majiq.internals2
  .. autosummary::
     :toctree: _generate

  )pbdoc";

  // vector<string> will be bound to Python list
  pybind11::bind_vector<std::vector<std::string>>(m, "VectorString");

  // define types bound to m
  auto pyGeneStrandness = pybind11::enum_<majiq::GeneStrandness>(
      m, "GeneStrandness")
    .value("forward", majiq::GeneStrandness::FORWARD)
    .value("reverse", majiq::GeneStrandness::REVERSE)
    .value("ambiguous", majiq::GeneStrandness::AMBIGUOUS);
  auto pyExperimentStrandness = pybind11::enum_<majiq::ExperimentStrandness>(
      m, "ExperimentStrandness")
    .value("FORWARD", majiq::ExperimentStrandness::FORWARD)
    .value("REVERSE", majiq::ExperimentStrandness::REVERSE)
    .value("NONE", majiq::ExperimentStrandness::NONE);
  auto pyExperimentThresholds = majiq::bindings::pyExperimentThresholds_t(
      m, "ExperimentThresholds",
      "Thresholds on intron/junction coverage for inclusion in splicegraph");
  auto pyIntronThresholdsGenerator
    = majiq::bindings::pyIntronThresholdsGenerator_t(
      m, "IntronThresholdsGenerator",
      "Functor to return intron thresholds given intron length");
  auto pyContigs = majiq::bindings::pyContigs_t(
      m, "Contigs",
      "Independent contigs on a genome");
  auto pyGenes = majiq::bindings::pyGenes_t(
      m, "Genes",
      "Independent genes defined on contigs");
  auto pyExons = majiq::bindings::pyExons_t(
      m, "Exons",
      "Exons defined on genes");
  auto pyGeneIntrons = majiq::bindings::pyGeneIntrons_t(
      m, "GeneIntrons",
      "Retained introns defined on genes between exons");
  auto pyGeneJunctions = majiq::bindings::pyGeneJunctions_t(
      m, "GeneJunctions",
      "Junctions defined on genes connecting exons");
  auto pyGeneJunctionsAccumulator
    = majiq::bindings::pyGeneJunctionsAccumulator_t(
      m, "GeneJunctionsAccumulator", "Accumulator of GeneJunctions");
  auto pySJIntrons = majiq::bindings::pySJIntrons_t(
      m, "SJIntrons",
      "Introns defined on contigs for which coverage can be defined");
  auto pySJIntronsBins = majiq::bindings::pySJIntronsBins_t(
      m, "SJIntronsBins",
      "Summarized and per-bin counts for introns from an experiment");
  auto pySJJunctions = majiq::bindings::pySJJunctions_t(
      m, "SJJunctions",
      "Junctions defined on contigs for which coverage can be defined");
  auto pySJJunctionsBins = majiq::bindings::pySJJunctionsBins_t(
      m, "SJJunctionsBins",
      "Summarized and per-bin counts for junctions from an experiment");
  auto pyGroupJunctionsGenerator = majiq::bindings::pyGroupJunctionsGenerator_t(
      m, "GroupJunctionsGenerator",
      "Accumulator of SJJunctions in the same build group");
  auto pyPassedJunctionsGenerator
    = majiq::bindings::pyPassedJunctionsGenerator_t(
      m, "PassedJunctionsGenerator",
      "Accumulator of GroupJunctionsGenerator from multiple build groups");
  auto pyGroupIntronsGenerator = majiq::bindings::pyGroupIntronsGenerator_t(
      m, "GroupIntronsGenerator",
      "Accumulator of SJIntronsBins in the same build group");
  auto pyEvents = majiq::bindings::pyEvents_t(
      m, "Events",
      "Events from reference exon with junctions/introns");
  auto pyEventsAlign = majiq::bindings::pyEventsAlign_t(
      m, "EventsAlign",
      "Indexes to shared events between two Events containers");
  auto pyEventsCoverage = majiq::bindings::pyEventsCoverage_t(
      m, "EventsCoverage",
      "Coverage over events for a single experiment");
  auto pySpliceGraphReads = majiq::bindings::pySpliceGraphReads_t(
      m, "SpliceGraphReads",
      "Raw readrates for each intron and junction");
  auto pyExonConnections = majiq::bindings::pyExonConnections_t(
      m, "ExonConnections",
      "Connections from exons to junctions, introns");
  auto pySimplifierGroup = majiq::bindings::pySimplifierGroup_t(
      m, "SimplifierGroup",
      "Accumulator of SpliceGraphReads to unsimplify introns and junctions");
  auto pySpliceGraphMask = majiq::bindings::pySpliceGraphMask_t(
      m, "SpliceGraphMask",
      "Boolean mask over splicegraph junctions and introns");
  auto pyGeneModules = majiq::bindings::pyGeneModules_t(
      m, "GeneModules",
      "Splicing modules for each gene");
  auto pySpliceGraph = majiq::bindings::pySpliceGraph_t(
      m, "SpliceGraph",
      "SpliceGraph managing exons, junctions, and introns within genes");

  // make updates to m
  majiq::bindings::init_RNGfunctions(m);
  majiq::bindings::init_GFF3Types(m);
  majiq::bindings::init_ExperimentThresholds(pyExperimentThresholds);
  majiq::bindings::init_IntronThresholdsGenerator(
      pyIntronThresholdsGenerator);
  majiq::bindings::init_Contigs(pyContigs);
  majiq::bindings::init_Genes(pyGenes);
  majiq::bindings::init_Exons(pyExons);
  majiq::bindings::init_GeneJunctions(pyGeneJunctions);
  majiq::bindings::init_GeneJunctionsAccumulator(pyGeneJunctionsAccumulator);
  majiq::bindings::init_GeneIntrons(pyGeneIntrons);
  majiq::bindings::init_SJJunctions(pySJJunctions);
  majiq::bindings::init_SJJunctionsBins(pySJJunctionsBins);
  majiq::bindings::init_SJIntrons(pySJIntrons);
  majiq::bindings::init_SJIntronsBins(pySJIntronsBins);
  majiq::bindings::init_GroupJunctionsGenerator(pyGroupJunctionsGenerator);
  majiq::bindings::init_PassedJunctionsGenerator(pyPassedJunctionsGenerator);
  majiq::bindings::init_GroupIntronsGenerator(pyGroupIntronsGenerator);
  majiq::bindings::init_Events(pyEvents);
  majiq::bindings::init_EventsAlign(pyEventsAlign);
  majiq::bindings::init_EventsCoverage(pyEventsCoverage);
  majiq::bindings::init_SpliceGraphReads(pySpliceGraphReads);
  majiq::bindings::init_ExonConnections(pyExonConnections);
  majiq::bindings::init_SimplifierGroup(pySimplifierGroup);
  majiq::bindings::init_SpliceGraphMask(pySpliceGraphMask);
  majiq::bindings::init_GeneModules(pyGeneModules);
  majiq::bindings::init_SpliceGraph(pySpliceGraph);
}  // PYBIND11_MODULE
