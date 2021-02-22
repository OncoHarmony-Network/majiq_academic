/**
 * SJJunctionsPositions.cpp
 *
 * Implementation of non-inline functions in SJJunctionsPositions
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#include "SJJunctionsPositions.hpp"

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

bool SJJunctionsPositions::is_valid(
    const std::shared_ptr<SJJunctions>& junctions,
    const std::vector<PositionReads>& reads,
    const std::vector<size_t>& offsets,
    junction_pos_t num_positions) {
  if (
      // must have junctions
      junctions == nullptr
      // offsets[i], offsets[i + 1] correspond to junctions[i] --> match sizes
      || junctions->size() != offsets.size() - 1
      // last offset must correspond to size of reads
      || offsets.back() != reads.size()) {
    return false;
  }
  // for each junction
  for (size_t jidx = 0; jidx < junctions->size(); ++jidx) {
    // offsets into junction positions
    auto jp_start = reads.begin() + offsets[jidx];
    auto jp_end = reads.begin() + offsets[jidx + 1];
    if (
        // offsets must be non-decreasing
        jp_end < jp_start
        // difference in offsets correspond to numpos
        || jp_end - jp_start != (*junctions)[jidx].numpos()
        // can't have numpos greater than total possible
        || (*junctions)[jidx].numpos() > num_positions
        // all junction positions must be in [0, num_positions)
        || !std::all_of(jp_start, jp_end,
          [num_positions](const PositionReads& x) {
          return x.pos < num_positions;
          })
        // sum of position reads equals junction numreads
        || std::accumulate(jp_start, jp_end, junction_ct_t{},
          [](junction_ct_t s, const PositionReads& x) { return s + x.reads; })
          != (*junctions)[jidx].numreads()
        // PositionReads in sorted order (by number of reads)
        || !std::is_sorted(jp_start, jp_end)) {
      return false;
    }
  }
  return true;
}

SJJunctionsPositions SJJunctionsPositions::FromBam(
    const char* infile, ExperimentStrandness exp_strandness, int nthreads) {
  auto contigs = std::make_shared<Contigs>();
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
  std::vector<PositionReads> reads(num_junction_positions);
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
          return PositionReads{x.first, x.second};
          });
      std::sort(reads.begin() + jpidx, reads.begin() + jpidx + num_positions);
      // summarize total reads
      const junction_ct_t num_reads = std::accumulate(
          reads.begin() + jpidx, reads.begin() + jpidx + num_positions,
          junction_ct_t{},
          [](junction_ct_t acc, const PositionReads& x) {
          return acc + x.reads;
          });
      // update junctions, offsets
      junctions[jidx] = SJJunction{contig, pos_strand.first, pos_strand.second,
                                   ExperimentCounts{num_reads, num_positions}};
      reads_offsets[jidx + 1] = reads_offsets[jidx] + num_positions;
      // increment jidx, jpidx
      ++jidx;
      jpidx += num_positions;
    }
  }
  // finally, create desired result
  return SJJunctionsPositions{
      std::make_shared<SJJunctions>(std::move(junctions)),
      reads, reads_offsets, eff_len};
}

void SJIntronsBins::check_valid() const {
  if (offsets_.empty()) {
    throw std::runtime_error("SJIntronsBins cannot have empty offsets");
  } else if (introns_.size() != offsets_.size() - 1) {
    throw std::runtime_error(
        "SJIntronsBins offsets do not correspond to number of introns");
  } else if (offsets_.back() != reads_.size()) {
    throw std::runtime_error(
        "SJIntronsBins offsets do not correspond to reads per bin");
  } else {
    for (size_t jidx = 0; jidx < introns_.size(); ++jidx) {
      if (offsets_[jidx] > offsets_[jidx + 1]) {
        throw std::runtime_error(
            "SJIntronsBins has nonmonotonically increasing offsets");
      } else if (offsets_[jidx + 1] - offsets_[jidx] > num_bins_) {
        throw std::runtime_error(
            "SJIntronsBins has reads for more bins than possible");
      } else if (
          !std::all_of(
            reads_.begin() + offsets_[jidx],
            reads_.begin() + offsets_[jidx + 1],
            [num_bins = num_bins_](const PositionReads& x) {
            return x.pos < num_bins; })) {
        throw std::runtime_error(
            "SJIntronsBins has intron with invalid read position");
      }
    }  // end loop over introns
  }
}

SJIntronsBins SJIntronsBins::FromBam(const char* infile,
    const junction_pos_t num_bins, const Exons& exons,
    const Introns& gene_introns, ExperimentStrandness exp_strandness,
    int nthreads) {
  // get introns we will be quantifying against
  ContigIntrons introns = ContigIntrons::FromGeneExonsAndIntrons(
      exons, gene_introns, exp_strandness != ExperimentStrandness::NONE);
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
          const ContigIntron& intron) {
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
          const junction_pos_t min_pos_per_bin
            = (num_bins + intron_length) / num_bins;
          const junction_pos_t bins_with_extra
            = (num_bins + intron_length) % num_bins;
          const junction_pos_t first_pos_without_extra
            = bins_with_extra * (1 + min_pos_per_bin);
          const junction_pos_t resulting_bin
            = raw_read_position < first_pos_without_extra
            ? raw_read_position / (1 + min_pos_per_bin)
            : (bins_with_extra
                + ((raw_read_position - first_pos_without_extra)
                  / min_pos_per_bin));
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
  std::vector<PositionReads> reads(ct_nonzero);
  size_t read_idx = 0;
  for (size_t intron_idx = 0; intron_idx < introns.size(); ++intron_idx) {
    const size_t start_idx = read_idx;
    for (junction_pos_t bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
      // get count for intron/bin
      const junction_ct_t& bin_ct = counts_mut(intron_idx, bin_idx);
      if (bin_ct > 0) {
        // nonzero counts get added to reads
        reads[read_idx++] = PositionReads{bin_idx, bin_ct};
      }
    }
    // get reads in sorted order per intron
    std::sort(reads.begin() + start_idx, reads.begin() + read_idx);
    // offset for this intron
    offsets[1 + intron_idx] = read_idx;
  }
  return SJIntronsBins(
      std::move(introns), std::move(reads), std::move(offsets), num_bins);
}

}  // namespace majiq
