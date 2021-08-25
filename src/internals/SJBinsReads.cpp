/**
 * SJBinsReads.cpp
 *
 * Implementation of non-inline functions in SJBinsReads
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "SJBinsReads.hpp"

#include <stdexcept>
#include <optional>
#include <memory>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>

#include "SJJunctions.hpp"
#include "bam/Alignments.hpp"
#include "bam/CigarRegions.hpp"


namespace majiq {

SJJunctionsBins SJJunctionsBins::FromBam(
    const char* infile, ExperimentStrandness exp_strandness, int nthreads) {
  auto contigs = Contigs::create();
  junction_pos_t max_read_length = 0;
  size_t num_junctions = 0;
  size_t num_junction_positions = 0;
  // counts per contig, junction, strand, position
  std::vector<  // per contig
    std::map<  // per junction/strand
      std::pair<OpenInterval, GeneStrandness>,
      std::map<junction_pos_t, junction_ct_t>
    >
  > counts;

  // open up bam and process
  {
    // open up bam for reading with specified number of threads
    bam::AlignmentFile in(infile, nthreads);

    // track contigs -- just copy over contigs from alignment file
    std::for_each(in.target_name(), in.target_name() + in.n_targets(),
        [&contigs](char* seqid) { contigs->add(Contig{std::string{seqid}}); });
    // expand junctions for all contigs
    counts.resize(in.n_targets());
    // iterate over alignments
    bam::AlignmentRecord aln{};
    int r = 0;
    // while we successfully read in a new alignment
    while ((r = in.NextRead(aln)) >= 0) {
      // update max read length
      max_read_length = std::max(
          max_read_length, static_cast<junction_pos_t>(aln.read_length()));
      // only process if unique/mapped
      if (!aln.unique_mapped()) { continue; }
      // parse junctions from uniquely mapped alignments
      bam::CigarRegions aln_regions = aln.cigar_regions();
      for (auto it = aln_regions.begin(); it != aln_regions.end(); ++it) {
        const auto& region = *it;
        if (region.position_
            > aln_regions.alignment_length_ - USE_MIN_OVERHANG) {
          break;
        } else if (region.position_ >= USE_MIN_OVERHANG
            && (region.type_
              == bam::CigarRegions::RegionType::OFF_GENOME_JUNCTION)) {
          // counts for this junction
          auto& junction_counts
            = counts[aln.tid()][std::make_pair(
                region.coordinates_.AsOpen(), aln.strandness(exp_strandness))];
          // increment num_junctions if we haven't seen it before
          if (junction_counts.empty()) { ++num_junctions; }
          // increment position, and num_junction_positions if new
          auto use_pos = region.position_
            + aln_regions.clipping_lengths_.first - USE_MIN_OVERHANG;
          if (1 == ++junction_counts[use_pos]) {
            ++num_junction_positions;
          }
        }
      }  // added each junction
    }  // while successfully loading next read
    if (r < -1) {
      std::ostringstream oss;
      oss << "Error parsing read from BAM " << infile
        << " (error code: " << r << ")";
      // throw std::runtime_error(oss.str());
      std::cerr << oss.str() << std::endl;
    }  // otherwise r = 0 and at end of file
  }  // done processing BAM file

  // number of valid positions (note that type may be unsigned)
  junction_pos_t eff_len
    = (max_read_length + 1 > 2 * USE_MIN_OVERHANG)
    ? max_read_length + 1 - 2 * USE_MIN_OVERHANG : 0;
  // initialize vectors for outputs
  std::vector<SJJunction> junctions(num_junctions);
  std::vector<BinReads<junction_ct_t>> reads(num_junction_positions);
  std::vector<size_t> reads_offsets(num_junctions + 1);
  size_t jidx = 0;
  size_t jpidx = 0;
  for (size_t contig_idx = 0; contig_idx < counts.size(); ++contig_idx) {
    if (counts[contig_idx].empty()) { continue; }
    KnownContig contig = (*contigs)[contig_idx];
    for (auto& [pos_strand, position_reads] : counts[contig_idx]) {
      const junction_pos_t num_positions = position_reads.size();
      // copy over position reads and sort in increasing order of read counts
      std::transform(position_reads.begin(), position_reads.end(),
          reads.begin() + jpidx,
          [](const std::pair<junction_pos_t, junction_ct_t>& x) {
          return BinReads<junction_ct_t>{x.first, x.second};
          });
      std::sort(reads.begin() + jpidx, reads.begin() + jpidx + num_positions);
      // update junctions, offsets
      junctions[jidx] = SJJunction{contig, pos_strand.first, pos_strand.second};
      reads_offsets[jidx + 1] = reads_offsets[jidx] + num_positions;
      // increment jidx, jpidx
      ++jidx;
      jpidx += num_positions;
    }
  }
  // finally, create desired result
  return SJJunctionsBins{
      std::make_shared<SJJunctions>(contigs, std::move(junctions)),
      std::move(reads), std::move(reads_offsets), eff_len};
}

SJIntronsBins SJIntronsBins::FromBam(const char* infile,
    const junction_pos_t num_bins, const Exons& exons,
    const GeneIntrons& gene_introns, ExperimentStrandness exp_strandness,
    int nthreads) {
  // get introns we will be quantifying against
  const auto introns_ptr = std::make_shared<SJIntrons>(
      SJIntrons::FromGeneExonsAndIntrons(exons, gene_introns,
        exp_strandness != ExperimentStrandness::NONE));
  const SJIntrons& introns = *introns_ptr;
  // get contigs for these introns
  auto contigs = introns.parents();
  // count number of reads per intron bin
  std::vector<junction_ct_t> counts2d(introns.size() * num_bins);
  auto counts_mut = [&counts2d, &num_bins](
      size_t intron_idx, junction_pos_t bin_idx) -> junction_ct_t& {
    return counts2d[intron_idx * num_bins + bin_idx];
  };
  // count number of nonzero bins
  size_t ct_nonzero = 0;

  // open up bam and process
  {
    bam::AlignmentFile in(infile, nthreads);
    // map target ids in BAM to known contigs
    std::vector<std::optional<KnownContig>> target_contigs(in.n_targets());
    for (size_t tid = 0; tid < target_contigs.size(); ++tid) {
      auto opt_contig_idx = contigs->safe_idx(in.target_name()[tid]);
      if (opt_contig_idx.has_value()) {
        target_contigs[tid] = (*contigs)[*opt_contig_idx];
      }
    }
    // iterate over alignments
    bam::AlignmentRecord aln{};
    int r = 0;
    while ((r = in.NextRead(aln)) >= 0) {
      // only process if unique/mapped
      if (!aln.unique_mapped()) { continue; }
      // only process if on target in contigs
      const auto& opt_contig = target_contigs[aln.tid()];
      if (!opt_contig.has_value()) { continue; }
      const KnownContig& contig = *opt_contig;
      // get strandness of alignment
      GeneStrandness aln_strand = aln.strandness(exp_strandness);
      // get regions from alignment
      bam::CigarRegions aln_regions = aln.cigar_regions();
      auto aln_it = aln_regions.begin();
      // get first intron that could overlap with this first region
      // (start - 1) because must handle zero-length introns where end < start
      auto intron_it
        = introns.overlap_lower_bound(contig, aln_it->coordinates_.start - 1);
      // If we have intron and last cigar region that intersects it, how to add?
      auto AddIntron = [&counts_mut, &ct_nonzero, &aln_regions, &num_bins](
          size_t intron_idx,
          const bam::CigarRegions::Region& cigar_region,
          const SJIntron& intron) {
        // assumes checked other conditions for validity besides overhang
        // (i.e. no junctions, correct strand, etc.)
        const junction_pos_t raw_alignment_position
          = cigar_region.position_
          + (cigar_region.type_ == bam::CigarRegions::RegionType::ON_GENOME
              ? intron.coordinates.end - cigar_region.coordinates_.start
              // i.e. non-junction deletions crossing intron end get last
              // position in intron before deletion
              : 0);
        // 1. alignment overlaps intron end by 1 + raw_alignment_position
        // bases
        // 2. there are intron_length + alignment length - 1 possible
        // positions, and overlap intron start by that - raw position
        const auto intron_length = intron.coordinates.length();
        if (
            // sufficient overlap on intron end
            1 + raw_alignment_position >= USE_MIN_OVERHANG
            // suffcient overlap on intron start
            && ((
                intron_length + aln_regions.alignment_length_ - 1
                - raw_alignment_position)
              >= USE_MIN_OVERHANG)) {
          // get raw read position, adjusting for overhang
          const junction_pos_t raw_read_position
            = raw_alignment_position + aln_regions.clipping_lengths_.first
            - (USE_MIN_OVERHANG - 1);
          // map to appropriate intron bin
          using detail::IntronCoarseBins;
          const junction_pos_t resulting_bin = (
              IntronCoarseBins(
                IntronCoarseBins::num_raw_positions(intron_length, num_bins),
                num_bins)
              .resulting_bin(raw_read_position));
          // update counts for this intron/bin
          if (1 == (++counts_mut(intron_idx, resulting_bin))) {
            // if first count for this intron/bin, update ct_nonzero
            ++ct_nonzero;
          }
        }
      };
      // loop over cigar regions
      // need last intersecting region when last region has end before intron's
      bam::CigarRegions::Region last_intersecting_region;
      while (aln_it != aln_regions.end()) {
        if (intron_it == introns.end() || intron_it->contig != contig) {
          // no more introns possible for read
          break;
        } else if (intron_it->strand != aln_strand) {
          // this intron is wrong strand
          ++intron_it;
          continue;
        }
        const bool intron_precedes_region
          = IntervalPrecedes(intron_it->coordinates, aln_it->coordinates_);
        if (intron_precedes_region) {
          // next intron might intersect or come after region
          ++intron_it;
          continue;
        }
        const bool intron_intersects_region
          = IntervalIntersects(aln_it->coordinates_, intron_it->coordinates);
        if (!intron_intersects_region) {
          // since we also know that intron does not precede region, we know
          // that region precedes current intron, so we must advance region
          ++aln_it;
          continue;
        } else if (aln_it->type_
            == bam::CigarRegions::RegionType::OFF_GENOME_JUNCTION) {
          // cannot count introns where junction is present inside it
          ++intron_it;
          continue;
        } else if (aln_it->coordinates_.end < intron_it->coordinates.end) {
          // we might have a better region with which to determine intron
          // position (or a junction that disqualifies the intron)
          last_intersecting_region = *aln_it;  // keep in case is best region
          ++aln_it;
          continue;
        } else {
          // this is definitely last cigar region intersecting the intron
          // try to add it
          AddIntron(intron_it - introns.begin(), *aln_it, *intron_it);
          ++intron_it;  // go to next intron
          continue;
        }
      }  // loop over cigar regions (can exit from break or ++aln_it)
      if (intron_it != introns.end() && intron_it->contig == contig
          && IntervalIntersects(
            last_intersecting_region.coordinates_, intron_it->coordinates)) {
        // intron intersects but continues past end of alignment, add it here
        AddIntron(
            intron_it - introns.begin(), last_intersecting_region, *intron_it);
      }
    }  // end reading alignments
    if (r < -1) {
      std::ostringstream oss;
      oss << "Error parsing read from BAM " << infile
        << " (error code: " << r << ")";
      // throw std::runtime_error(oss.str());
      std::cerr << oss.str() << std::endl;
    }  // otherwise r = 0 and at end of file
  }  // done processing BAM file
  // get offsets vector
  std::vector<size_t> offsets(1 + introns.size());
  std::vector<BinReads<intron_ct_t>> reads(ct_nonzero);
  size_t read_idx = 0;
  for (size_t intron_idx = 0; intron_idx < introns.size(); ++intron_idx) {
    const size_t start_idx = read_idx;
    // get scaling to scale raw counts to average number of positions per bin
    using detail::IntronCoarseBins;
    const auto total_positions = IntronCoarseBins::num_raw_positions(
        introns[intron_idx].coordinates.length(), num_bins);
    const real_t avg_num_positions
      = static_cast<real_t>(total_positions) / num_bins;
    const auto binning = IntronCoarseBins(total_positions, num_bins);
    for (junction_pos_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
      // get count for intron/bin
      const junction_ct_t& bin_ct = counts_mut(intron_idx, bin_idx);
      if (bin_ct > 0) {
        // nonzero counts get added to reads
        reads[read_idx++] = BinReads<intron_ct_t>{
            bin_idx,
            // scale counts to average number of positions
            bin_ct * avg_num_positions / binning.bin_num_positions(bin_idx)};
      }
    }
    // get reads in sorted order per intron
    std::sort(reads.begin() + start_idx, reads.begin() + read_idx);
    // offset for this intron
    offsets[1 + intron_idx] = read_idx;
  }
  return SJIntronsBins(
      introns_ptr, std::move(reads), std::move(offsets), num_bins);
}

}  // namespace majiq