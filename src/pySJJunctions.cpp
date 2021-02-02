/**
 * pySJJunctions.cpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <array>
#include <vector>
#include <memory>
#include <cstddef>

#include "pybind_utils.hpp"
#include "internals/SJJunctions.hpp"
#include "internals/SJJunctionsPositions.hpp"
#include "internals/IOBamJunctions.hpp"


namespace py = pybind11;
using PySJJunctionsPositionsT
  = py::class_<majiq::SJJunctionsPositions,
      std::shared_ptr<majiq::SJJunctionsPositions>>;


void enable_IOBamJunctions(PySJJunctionsPositionsT& pySJJunctionsPositions) {
  using majiq::ExperimentStrandness;
  pySJJunctionsPositions
    .def_static("from_bam",
        [](const char* infile, ExperimentStrandness strand, int nthreads) {
        return majiq::JunctionsFromBam<BAM_MIN_OVERHANG>(
            infile, nthreads, strand);
        },
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
        py::arg("experiment_strandness") = ExperimentStrandness::NONE,
        py::arg("nthreads") = 1);
}
