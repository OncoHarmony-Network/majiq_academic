/**
 * SimplifierGroup.hpp
 *
 * Use SpliceGraphReads over a group of experiments to *unsimplify* introns and
 * junctions
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_SIMPLIFIERGROUP_HPP
#define MAJIQ_SIMPLIFIERGROUP_HPP

#include <memory>
#include <mutex>
#include <shared_mutex>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "SpliceGraphReads.hpp"
#include "ExonConnections.hpp"


namespace majiq {

// number of experiments for src/dst events with evidence for unsimplifying
struct SimplifierCount {
  size_t src_ct_unsimplify;
  size_t dst_ct_unsimplify;

  void reset() {
    src_ct_unsimplify = dst_ct_unsimplify = 0;
  }
  void increment(const EventType& type) {
    // increment appropriate value
    switch (type) {
      case EventType::SRC_EVENT:
        ++src_ct_unsimplify;
        break;
      case EventType::DST_EVENT:
        ++dst_ct_unsimplify;
        break;
    }
  }
  // we keep simplified if src, dst cts are both below min_experiments
  bool keep_simplified(size_t min_experiments) const {
    return (src_ct_unsimplify < min_experiments)
      && (dst_ct_unsimplify < min_experiments);
  }

  SimplifierCount() { reset(); }
  SimplifierCount(const SimplifierCount&) = default;
  SimplifierCount(SimplifierCount&&) = default;
  SimplifierCount& operator=(const SimplifierCount&) = default;
  SimplifierCount& operator=(SimplifierCount&&) = default;
};

class SimplifierGroup {
 private:
  const std::shared_ptr<ExonConnections> exon_connections_;

  size_t num_experiments_;
  std::vector<SimplifierCount> introns_passed_;
  std::vector<SimplifierCount> junctions_passed_;
  // mutexes
  std::mutex simplifier_mutex_;  // exclusive access to write to above variables

 public:
  explicit SimplifierGroup(
      const std::shared_ptr<ExonConnections>& exon_connections)
      : exon_connections_{exon_connections},
        num_experiments_{0},
        introns_passed_(
            exon_connections_ == nullptr
            ? 0 : exon_connections_->introns()->size()),
        junctions_passed_(
            exon_connections_ == nullptr
            ? 0 : exon_connections_->junctions()->size()),
        simplifier_mutex_{} {
    if (exon_connections_ == nullptr) {
      throw std::runtime_error("SimplifierGroup given null exon connections");
    }
  }

  const std::shared_ptr<ExonConnections>& exon_connections() const {
    return exon_connections_;
  }
  const size_t& num_experiments() const { return num_experiments_; }
  const std::vector<SimplifierCount>& introns_passed() const {
    return introns_passed_;
  }
  const std::vector<SimplifierCount>& junctions_passed() const {
    return junctions_passed_;
  }

  void UpdateInplace(real_t min_experiments_f) {
    std::lock_guard lock{simplifier_mutex_};
    size_t min_experiments = detail::min_experiments_from_float(
        num_experiments_, min_experiments_f);
    const GeneIntrons& introns = *(exon_connections_->introns());
    const GeneJunctions& junctions = *(exon_connections_->junctions());
    for (size_t i = 0; i < junctions.size(); ++i) {
      // if currently simplified but evidence says not to keep simplified...
      if (junctions[i].simplified()
          && !junctions_passed_[i].keep_simplified(min_experiments)) {
        junctions[i].simplified() = false;
      }
      // either way, reset the counts for this junction
      junctions_passed_[i].reset();
    }
    for (size_t i = 0; i < introns.size(); ++i) {
      // if currently simplified but evidence says not to keep simplified...
      if (introns[i].simplified()
          && !introns_passed_[i].keep_simplified(min_experiments)) {
        introns[i].simplified() = false;
      }
      // either way, reset the counts for this intron
      introns_passed_[i].reset();
    }
    // reset num_experiments_
    num_experiments_ = 0;
  }

  void AddExperiment(
      const SpliceGraphReads& sg_reads,
      real_t min_psi,
      real_t minreads_annotated_junction,
      real_t minreads_denovo_junction,
      real_t minreads_intron) {
    // match introns/junctions
    if (exon_connections_->introns() != sg_reads.introns()) {
      throw std::runtime_error(
          "SimplifierGroup given experiment reads with different introns");
    } else if (exon_connections_->junctions() != sg_reads.junctions()) {
      throw std::runtime_error(
          "SimplifierGroup given experiment reads with different junctions");
    }
    const GeneJunctions& junctions = *(exon_connections_->junctions());
    const Exons& exons = *(exon_connections_->exons());

    std::lock_guard lock{simplifier_mutex_};
    ++num_experiments_;  // update number of experiments

    // get evidence to unsimplify per potential event
    constexpr std::array<EventType, 2> TYPES
      = {EventType::SRC_EVENT, EventType::DST_EVENT};
    // for each reference exon/direction...
    for (size_t exon_idx = 0;
        exon_idx < exon_connections_->num_exons(); ++exon_idx) {
      for (const auto& type : TYPES) {
        Event event{exon_idx, type};
        const auto junction_idx_begin
          = exon_connections_->begin_junctions_for(event);
        const auto junction_idx_end
          = exon_connections_->end_junctions_for(event);
        const auto intron_idx_begin
          = exon_connections_->begin_introns_for(event);
        const auto intron_idx_end
          = exon_connections_->end_introns_for(event);
        // get total numreads to calculate min_psi
        // NOTE: this includes exitrons, simplified connections
        real_t total_reads
          = std::accumulate(junction_idx_begin, junction_idx_end, real_t{0},
              [&x = sg_reads.junctions_reads()](real_t s, size_t i) {
              return s + x[i]; })
          + std::accumulate(intron_idx_begin, intron_idx_end, real_t{0},
              [&x = sg_reads.introns_reads()](real_t s, size_t i) {
              return s + x[i]; });
        if (total_reads == 0) { continue; }  // no evidence for this event
        // if this is a half exon, total_reads does not consider intersecting
        // intron coverage. This coverage should be added to total_reads.
        if (!exons[exon_idx].is_full_exon()) {
          // get next full exon
          size_t full_exon_idx = 1 + exon_idx;
          bool full_exon_found = false;
          for (;
              full_exon_idx < exons.size()
              && exons[full_exon_idx].gene.idx_ == exons[exon_idx].gene.idx_;
              ++full_exon_idx) {
            if (exons[full_exon_idx].is_full_exon()) {
              full_exon_found = true;
              break;
            }
          }
          if (full_exon_found) {
            Event full_event{full_exon_idx,
                // if forward strand gene, then to get intron overlapping
                // exon_idx, we would need to look at this next exon as a
                // target event
                exons[full_exon_idx].strand_forward()
                ? EventType::DST_EVENT
                : EventType::SRC_EVENT};
            total_reads = std::accumulate(
                exon_connections_->begin_introns_for(full_event),
                exon_connections_->end_introns_for(full_event),
                total_reads,
                [&x = sg_reads.introns_reads()](real_t s, size_t i) {
                return s + x[i]; });
          }
        }
        // determine if experiment has evidence to unsimplify the connections
        for (auto it = junction_idx_begin; it != junction_idx_end; ++it) {
          // reads for this junction index
          const real_t& reads = sg_reads.junctions_reads()[*it];
          // minreads for this junction index
          const real_t& minreads = junctions[*it].denovo()
            ? minreads_denovo_junction : minreads_annotated_junction;
          // is there enough evidence to unsimplify?
          if (reads >= minreads && reads >= min_psi * total_reads) {
            junctions_passed_[*it].increment(type);
          }
        }  // increment evidence for junctions in event
        for (auto it = intron_idx_begin; it != intron_idx_end; ++it) {
          // reads for this intron index
          const real_t& reads = sg_reads.introns_reads()[*it];
          // is there enough evidence to unsimplify?
          if (reads >= minreads_intron && reads >= min_psi * total_reads) {
            introns_passed_[*it].increment(type);
          }
        }  // increment evidence for introns in event
      }  // loop over event types for a reference exon
    }  // loop over reference exons
  }
};

}  // namespace majiq


#endif  // MAJIQ_SIMPLIFIERGROUP_HPP
