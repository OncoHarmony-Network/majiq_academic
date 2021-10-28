/**
 * PassedJunctions.cpp
 *
 * Implementation of assignment of SJJunctions to genes
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "PassedJunctions.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <numeric>


namespace majiq {


/**
 * Assign denovo junction junction to potential matched genes in range [first,
 * last)
 *
 * Track genes in og that junction could be assigned to. So, we prioritize
 * junctions that are close (or overlapping) an exon in the gene, preferring
 * genes where the junction is close on both ends over just one end vs neither
 */
void AssignDenovoJunction(const SJJunction& junction, const KnownGene& first,
    const KnownGene& last, const Exons& exons,
    std::map<GeneJunction, size_t>& dst) {
  // when junction is close to exons on both junction ends as appropriate
  std::vector<KnownGene> close_x2;
  // when only one end of junction close to a gene exon
  std::vector<KnownGene> close_x1;
  // when within gene boundaries but not close to exons at all
  std::vector<KnownGene> far;

  // check all of the overgenes
  for (KnownGene gene = first; gene != last; ++gene) {
    // if strand doesn't match, ignore
    if (!(junction.strand == GeneStrandness::AMBIGUOUS
          || junction.strand == gene.strand())) {
      continue;
    }
    // if junction cannot be placed inside of gene, ignore
    {
      const auto exon_first = exons.begin_parent(gene);
      const auto exon_last = exons.end_parent(gene);
      if (exon_first == exon_last
          || !IntervalSubsets(junction.coordinates,
            ClosedInterval{exon_first->coordinates.first_pos(),
              std::prev(exon_last)->coordinates.last_pos()})) {
        continue;
      }
    }
    // count how many of start/end are close to exon
    int n_close = 0;  // how many of start/end close to an exon?
    if (detail::CloseToPrecedingAnnotatedExon(
          exons, gene, junction.coordinates.start)) {
      ++n_close;
    }
    if (detail::CloseToFollowingAnnotatedExon(
          exons, gene, junction.coordinates.end)) {
      ++n_close;
    }
    // priority of this gene/where we track it
    auto& gene_priority_vec
      = n_close == 0 ? far : (n_close == 1 ? close_x1 : close_x2);
    gene_priority_vec.push_back(gene);
  }  // done looking at all potential genes
  // so get the highest priority vector of genes that is nonempty
  const auto& matched_genes
    = close_x2.empty() ? (close_x1.empty() ? far : close_x1) : close_x2;
  // add count to dst for this junction in those genes, if any
  std::for_each(matched_genes.begin(), matched_genes.end(),
      [&junction, &dst](const KnownGene& gene) {
      auto key = GeneJunction{gene, junction.coordinates, true, true, false};
      ++dst[key];
      });
  return;
}

void GroupJunctionsGenerator::AddExperiment(const SJJunctionsBins& sjp,
    const ExperimentThresholds& thresholds, bool process_denovo) {
  std::shared_lock group_lock{group_mutex_};
  {  // increment number of experiments
    std::lock_guard num_experiments_lock{num_experiments_mutex_};
    ++num_experiments_;
  }
  const SJJunctions& sj = *(sjp.regions());
  if (sj.empty()) { return; }  // no junctions to add

  const std::shared_ptr<Genes>& genes = exons_->parents();

  // for each contig we have genes for
  for (size_t dst_contig_idx = 0;
      dst_contig_idx < genes->contigs()->size();
      ++dst_contig_idx) {
    // try to get matching contigs in sj (since often does not share contigs)
    const auto opt_sj_contig_idx = sj.parents()->safe_idx(
        genes->contigs()->get(dst_contig_idx));
    if (!opt_sj_contig_idx.has_value()) { continue; }
    const size_t& sj_contig_idx = *opt_sj_contig_idx;

    // get gene junctions (or indexes to them) in contig sorted order
    // so we iterate over sj, gene junctions (known) and genes (for denovo)
    auto sj_it = sj.begin_parent(sj_contig_idx);
    const auto sj_it_end = sj.end_parent(sj_contig_idx);
    if (sj_it == sj_it_end) { continue; }  // no sj junctions for contig
    KnownGene g_it_start = genes->begin_contig(dst_contig_idx);
    const KnownGene g_it_end = genes->end_contig(dst_contig_idx);
    if (g_it_start == g_it_end) { continue; }  // no genes for contig

    // (indexes to) known genejunctions for this contig sorted by coordinates
    std::vector<size_t> known_idx(
        // number of known junctions for the contig
        known_->parent_idx_offsets_[g_it_end.idx_]
        - known_->parent_idx_offsets_[g_it_start.idx_]);
    std::iota(known_idx.begin(), known_idx.end(),
        known_->parent_idx_offsets_[g_it_start.idx_]);
    std::sort(known_idx.begin(), known_idx.end(),
        [&x = *known_](size_t i, size_t j) {
        // ignore gene, strand information, just coordinates within a contig
        return x[i].coordinates < x[j].coordinates; });
    auto k_it_start = known_idx.begin();

    // get denovos num passed for specified contig
    auto& denovos_num_passed = contig_denovos_num_passed_[dst_contig_idx];

    // lock on contig being updated
    std::lock_guard contig_lock{contig_mutex_[dst_contig_idx]};

    for (; sj_it != sj_it_end; ++sj_it) {
      const SJJunction& junction = *sj_it;
      JunctionPassedStatus status = sjp.passed(sj_it - sj.begin(), thresholds);
      if (status == JunctionPassedStatus::NOT_PASSED) { continue; }
      // advance to first known junction that could overlap with sj_it
      k_it_start = std::find_if(k_it_start, known_idx.end(),
          [&x = *known_, &junction](size_t i) {
          return !(x[i].coordinates < junction.coordinates); });
      bool found = false;  // was junction found in known?
      // try to find matching known junction
      for (auto k_it = k_it_start;
          // still have known junctions
          k_it != known_idx.end()
          // that have matching coordinates to junction
          && (*known_)[*k_it].coordinates == junction.coordinates;
          ++k_it) {
        const GeneJunction& kj = (*known_)[*k_it];  // gene junction for k_it
        // if strand matches and junction passed status ~ denovo status
        if ((junction.strand == GeneStrandness::AMBIGUOUS
             || junction.strand == kj.gene.strand())
            && (!kj.denovo() || status == JunctionPassedStatus::DENOVO_PASSED)
            ) {
          ++known_num_passed_[*k_it];
          found = true;
        }
      }  // loop over potentially matching known junctions
      if (!process_denovo) {
        // if no denovos, end of known junctions in contig is stopping condition
        if (k_it_start == known_idx.end()) { break; }
        // otherwise, we are processing junctions not previously known to us
      } else if (!found && status == JunctionPassedStatus::DENOVO_PASSED) {
        // advance to first gene that could overlap with junction
        g_it_start = std::find_if(g_it_start, g_it_end,
            [&junction](const KnownGene& x) {
            return !IntervalPrecedes(x.coordinates(), junction.coordinates); });
        if (g_it_start == g_it_end) { break; }  // no more genes for contig
        // gene can be from there until first gene that starts after junction
        auto g_it_after = std::find_if(g_it_start, g_it_end,
            [&junction](const KnownGene& x) {
            return IntervalPrecedes(junction.coordinates, x.coordinates()); });
        if (g_it_start != g_it_after) {
          AssignDenovoJunction(
              junction, g_it_start, g_it_after, *exons_, denovos_num_passed);
        }
      }
    }  // end loop over sj junctions on contig
  }  // end loop over contig
}

}  // namespace majiq
