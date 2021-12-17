/**
 * pySJBinsReads.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSJBINSREADS_HPP
#define MAJIQ_PYBIND_PYSJBINSREADS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <utility>
#include <vector>

#include "../SJBinsReads.hpp"

#include "constants.hpp"
#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pySJIntronsBins_t = pyClassShared_t<SJIntronsBins>;
using pySJJunctionsBins_t = pyClassShared_t<SJJunctionsBins>;

template <
  typename SJBinsT,
  typename RegionsT = base_t<decltype(*(std::declval<SJBinsT>().regions()))>,
  typename BinReads = base_t<decltype(std::declval<SJBinsT>().reads()[0])>,
  typename CountT = decltype(std::declval<BinReads>().bin_reads)>
inline void define_sjbins_properties(pyClassShared_t<SJBinsT>& pySJBins) {
  pySJBins
    .def_property_readonly("_regions", &SJBinsT::regions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying regions bin reads are defined over")
    .def_property_readonly("total_bins", &SJBinsT::total_bins,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "the total number of bins possible (positions for junctions, too)")
    .def_property_readonly("bin_reads",
        [](pybind11::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        const size_t offset = offsetof(BinReads, bin_reads);
        return ArrayFromVectorAndOffset<CountT, BinReads>(
            sj.reads(), offset, sj_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of reads for each bin")
    .def_property_readonly("bin_idx",
        [](pybind11::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        const size_t offset = offsetof(BinReads, bin_idx);
        return ArrayFromVectorAndOffset<majiq::junction_pos_t, BinReads>(
            sj.reads(), offset, sj_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Bin index for the each bin reads")
    .def_property_readonly("_offsets",
        [](pybind11::object& sj_obj) {
        SJBinsT& sj = sj_obj.cast<SJBinsT&>();
        return ArrayFromVectorAndOffset<size_t, size_t>(
            sj.offsets(), 0, sj_obj);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Raw offsets for regions into bin reads")
    .def("numstacks",
        [](const SJBinsT& self,
          pybind11::array_t<size_t> idx, pybind11::array_t<majiq::real_t> pvalue) {
        auto f = [&self](size_t i, majiq::real_t p) {
          if (i >= self.num_regions()) {
            throw std::invalid_argument("idx has out of range values");
          }
          return self.numstacks(i, p); };
        return pybind11::vectorize(f)(idx, pvalue);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Get number of stacks for the specified regions given threshold",
        pybind11::arg("region_idx"),
        pybind11::arg("pvalue_threshold") = DEFAULT_BUILD_STACK_PVALUE)
    .def("numbins",
        [](const SJBinsT& self,
          pybind11::array_t<size_t> idx,
          pybind11::array_t<CountT> minreads) {
        auto f = [&self](size_t i, CountT r) {
          if (i >= self.num_regions()) {
            throw std::invalid_argument("idx has out of range values");
          }
          return self.numbins_minreads(i, r); };
        return pybind11::vectorize(f)(idx, minreads);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of bins for regions with more than specified number of reads",
        pybind11::arg("region_idx"), pybind11::arg("minreads"))
    .def("numreads",
        [](const SJBinsT& self,
          pybind11::array_t<size_t> idx,
          pybind11::array_t<majiq::junction_pos_t> num_stacks) {
        auto f = [&self](size_t i, majiq::junction_pos_t n) {
          if (i >= self.num_regions()) {
            throw std::invalid_argument("idx has out of range values");
          }
          return self.numreads(i, n); };
        return pybind11::vectorize(f)(idx, num_stacks);
        },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Number of reads for regions given known number of stacks",
        pybind11::arg("region_idx"), pybind11::arg("num_stacks"))
    .def("__len__", &SJBinsT::size,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Total number of bin reads")
    .def(pybind11::init([](
            std::shared_ptr<RegionsT> regions,
            pybind11::array_t<CountT> _bin_reads,
            pybind11::array_t<majiq::junction_pos_t> _bin_idx,
            pybind11::array_t<size_t> _offsets,
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
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Initialize bins over specified regions with per-bin coverage",
        pybind11::arg("sj"),
        pybind11::arg("bin_reads"), pybind11::arg("bin_idx"), pybind11::arg("_offsets"),
        pybind11::arg("total_bins"));
}

inline void init_SJIntronsBins(pySJIntronsBins_t& pySJIntronsBins) {
  define_sjbins_properties(pySJIntronsBins);
  pySJIntronsBins
    .def_static("from_bam", &SJIntronsBins::FromBam,
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
        exons: Exons
            Gene exons defining potential introns for coverage
        gene_introns: GeneIntrons
            Gene introns indicating annotated introns for coverage
        experiment_strandness: ExperimentStrandness
            Strandness of RNA-seq library
        nthreads: int
            Number of threads to use when reading in BAM file
        )pbdoc",
        pybind11::arg("bam_path"),
        pybind11::arg("num_bins"),
        pybind11::arg("exons"),
        pybind11::arg("gene_introns"),
        pybind11::arg("experiment_strandness") = DEFAULT_BAM_STRANDNESS,
        pybind11::arg("nthreads") = DEFAULT_BAM_NTHREADS);
}

inline void init_SJJunctionsBins(pySJJunctionsBins_t& pySJJunctionsBins) {
  define_sjbins_properties(pySJJunctionsBins);
  pySJJunctionsBins
    .def("project_reads", &SJJunctionsBins::ProjectReads,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Project reads to a different set of SJJunctions")
    .def_static("from_bam", &SJJunctionsBins::FromBam,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
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
        pybind11::arg("bam_path"),
        pybind11::arg("experiment_strandness") = DEFAULT_BAM_STRANDNESS,
        pybind11::arg("nthreads") = DEFAULT_BAM_NTHREADS);
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSJBINSREADS_HPP
