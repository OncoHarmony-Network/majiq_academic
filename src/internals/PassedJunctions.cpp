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

#include "MajiqConstants.hpp"


namespace majiq {

inline bool CloseToPrecedingAnnotatedExon(
    const Exons& exons, const KnownGene& gene, position_t x) {
  // get first exon past x
  auto it = exons.overlap_upper_bound(gene, x);
  --it;  // first exon behind that (overlaps or behind)
  // keep going backwards until full exon, too far away, or different gene
  for (; it != exons.begin(); --it) {
    if (it->gene != gene
        || (it->is_full_exon() && !(it->is_denovo()))
        || (it->annotated_coordinates().last_pos()
          < x - MAX_DENOVO_DIFFERENCE)) {
      break;
    }
  }
  return (it->gene == gene
      && (it->is_full_exon() && !(it->is_denovo()))
      && x >= it->annotated_coordinates().start
      && x <= it->annotated_coordinates().end + MAX_DENOVO_DIFFERENCE);
}
inline bool CloseToFollowingAnnotatedExon(
    const Exons& exons, const KnownGene& gene, position_t x) {
  // get first exon overlapping or past x
  auto it = exons.overlap_lower_bound(gene, x);
  // keep going forwards until full exon, too far away, or different gene
  for (; it != exons.end(); ++it) {
    if (it->gene != gene
        || (it->is_full_exon() && !(it->is_denovo()))
        || (it->annotated_coordinates().first_pos()
          > x + MAX_DENOVO_DIFFERENCE)) {
      break;
    }
  }
  return (it != exons.end()
      && it->gene == gene
      && (it->is_full_exon() && !(it->is_denovo()))
      && x >= it->annotated_coordinates().start - MAX_DENOVO_DIFFERENCE
      && x <= it->annotated_coordinates().end);
}

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
    if (!IntervalSubsets(junction.coordinates, gene.coordinates())) {
      continue;
    }
    // count how many of start/end are close to exon
    int n_close = 0;  // how many of start/end close to an exon?
    if (CloseToPrecedingAnnotatedExon(
          exons, gene, junction.coordinates.start)) {
      ++n_close;
    }
    if (CloseToFollowingAnnotatedExon(
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

void GroupJunctionsGenerator::AddExperiment(
    const SJJunctions& sj, junction_ct_t minreads, junction_ct_t mindenovo,
    junction_pos_t minpos, bool process_denovo) {
  ++num_experiments_;  // increment number of experiments
  if (sj.empty()) { return; }  // no junctions to add, terminate early

  mindenovo = std::max(minreads, mindenovo);  // denovo junctions are stricter

  const std::shared_ptr<Genes>& genes = exons_->parents();

  // for each contig we have genes for
  for (size_t dst_contig_idx = 0; dst_contig_idx < genes->contigs()->size();
      ++dst_contig_idx) {
    const auto opt_sj_contig_idx
      = sj.parents()->safe_idx(genes->contigs()->get(dst_contig_idx));
    // if sj doesn't share dst contig, just skip it
    if (!opt_sj_contig_idx.has_value()) { continue; }
    const size_t& sj_contig_idx = *opt_sj_contig_idx;
    // first overgene for the contig
    auto og = OverGene{genes->begin_contig(dst_contig_idx)};
    const auto og_end = OverGene{genes->end_contig(dst_contig_idx)};

    // process each sj junction in order on this contig.
    for (auto sj_it = sj.begin_parent(sj_contig_idx);
        sj_it != sj.end_parent(sj_contig_idx);
        ++sj_it) {
      const SJJunction& junction = *sj_it;
      JunctionPassedStatus status
        = junction.passed(minreads, mindenovo, minpos);
      if (status == JunctionPassedStatus::NOT_PASSED) { continue; }
      // update og until could overlap with junction
      for (; og != og_end
          && og.coordinates.end < junction.coordinates.start;
          ++og) { }
      // if we have run out of overgenes in contig, futile to look at any more
      // junctions on the contigs
      if (og == og_end) { break; }
      // junction can only be added if it's part of the overgene
      if (!IntervalSubsets(junction.coordinates, og.coordinates)) { continue; }
      // handle annotated junctions for each gene in overgene
      bool known_found = false;
      std::for_each(og.first(), og.last(),
          [&known_found, this, &junction, &status](const KnownGene& gene) {
          auto key = GeneJunction{gene, junction.coordinates};
          auto match = known_->find(
              GeneJunction{gene, junction.coordinates});
          if (match != known_->end()
              // it's conceivable that we know denovo junctions, so if denovo,
              // we have to check that it passed denovo filters
              && (!match->denovo()
                || status == JunctionPassedStatus::DENOVO_PASSED)) {
            ++known_num_passed_[match - known_->begin()];
            known_found = true;
          }
          });
      // add if we are processing denovo junctions, is denovo, and passed denovo
      if (process_denovo
          && status == JunctionPassedStatus::DENOVO_PASSED
          && !known_found) {
        AssignDenovoJunction(
            junction, og.first(), og.last(), *exons_, denovos_num_passed_);
      }
    }  // loop over sj junctions for the contig
  }  // end loop over contigs for which we have genes
}

}  // namespace majiq
