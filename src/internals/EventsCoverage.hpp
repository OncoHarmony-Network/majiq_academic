/**
 * EventsCoverage.hpp
 *
 * Coverage for an experiment for selected Events
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_EVENTSCOVERAGE_HPP
#define MAJIQ_EVENTSCOVERAGE_HPP

#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Events.hpp"
#include "SJJunctionsPositions.hpp"


namespace majiq {

struct CoverageSummary {
  real_t numreads;
  real_t numbins;
  CoverageSummary(real_t _numreads, real_t _numbins)
      : numreads{_numreads}, numbins{_numbins} { }
  CoverageSummary() : CoverageSummary{0, 0} { }
  CoverageSummary(const CoverageSummary&) = default;
  CoverageSummary(CoverageSummary&&) = default;
  CoverageSummary& operator=(const CoverageSummary&) = default;
  CoverageSummary& operator=(CoverageSummary&&) = default;
};

class CoverageBootstraps {
 private:
  size_t num_bootstraps_;
  std::vector<real_t> coverage_;

 public:
  size_t num_bootstraps() const { return num_bootstraps_; }
  size_t num_connections() const { return coverage_.size() / num_bootstraps(); }

  const std::vector<real_t>& data() const { return coverage_; }
  const std::vector<real_t>& data() { return coverage_; }

  // random access into the matrix of bootstraps
  const real_t& operator()(size_t connection_idx, size_t bootstrap_idx) const {
    return coverage_[connection_idx * num_bootstraps_ + bootstrap_idx];
  }
  real_t& operator()(size_t connection_idx, size_t bootstrap_idx) {
    return coverage_[connection_idx * num_bootstraps_ + bootstrap_idx];
  }

  // slices for a single connection
  typename std::vector<real_t>::const_iterator
  begin_connection(size_t connection_idx) const {
    return coverage_.cbegin() + (connection_idx * num_bootstraps_);
  }
  typename std::vector<real_t>::const_iterator
  end_connection(size_t connection_idx) const {
    return begin_connection(1 + connection_idx);
  }
  typename std::vector<real_t>::iterator
  begin_connection(size_t connection_idx) {
    return coverage_.begin() + (connection_idx * num_bootstraps_);
  }
  typename std::vector<real_t>::iterator
  end_connection(size_t connection_idx) {
    return begin_connection(1 + connection_idx);
  }

  CoverageBootstraps(size_t num_connections, size_t num_bootstraps)
      : num_bootstraps_{num_bootstraps},
        coverage_(num_connections * num_bootstraps) { }
  CoverageBootstraps() : CoverageBootstraps{0, 0} { }
  CoverageBootstraps(const CoverageBootstraps&) = default;
  CoverageBootstraps(CoverageBootstraps&&) = default;
  CoverageBootstraps& operator=(const CoverageBootstraps&) = default;
  CoverageBootstraps& operator=(CoverageBootstraps&&) = default;
};

class EventsCoverage {
 private:
  std::shared_ptr<Events> events_;
  std::vector<CoverageSummary> summaries_;
  CoverageBootstraps bootstraps_;

 public:
  const std::shared_ptr<Events>& events() const { return events_; }
  const std::vector<CoverageSummary>& summaries() const {
    return summaries_;
  }
  const std::vector<CoverageSummary>& summaries() {
    return summaries_;
  }
  const CoverageBootstraps& bootstraps() const {
    return bootstraps_;
  }
  const CoverageBootstraps& bootstraps() {
    return bootstraps_;
  }
  size_t num_connections() const { return events_->num_connections(); }

  EventsCoverage(const std::shared_ptr<Events>& events,
      std::vector<CoverageSummary>&& summaries,
      CoverageBootstraps&& bootstraps)
      : events_{events},
        summaries_{std::move(summaries)},
        bootstraps_{std::move(bootstraps)} {
    if (events_ == nullptr) {
      throw std::runtime_error("Null events in EventsCoverage");
    } else if (summaries_.size() != events_->num_connections()) {
      throw std::runtime_error(
          "EventsCoverage summaries mismatched number of connections");
    } else if (bootstraps_.num_connections() != events_->num_connections()) {
      throw std::runtime_error(
          "EventsCoverage bootstraps mismatched number of connections");
    }
  }
  EventsCoverage(const EventsCoverage&) = default;
  EventsCoverage(EventsCoverage&&) = default;
  EventsCoverage& operator=(const EventsCoverage&) = default;
  EventsCoverage& operator=(EventsCoverage&&) = default;

  // how do we fill in coverage from raw junctions/introns?
 private:
  template <typename SJBinsT>
  static void merge_bins(
      const Events& events, const SJBinsT& sj,
      rng_t& generator, real_t pvalue_threshold,
      std::vector<CoverageSummary>& summaries, CoverageBootstraps& bootstraps) {
    static_assert(
        std::is_same_v<SJJunctionsPositions, SJBinsT>
        || std::is_same_v<SJIntronsBins, SJBinsT>,
        "SJBinsT must be SJJunctionsPositions or SJIntronsBins");
    constexpr bool IS_INTRON = std::is_same_v<SJIntronsBins, SJBinsT>;
    // get indexes of event connections to iterate over in contig order
    const auto connection_idx_begin = events.connection_idx_begin<IS_INTRON>();
    const auto connection_idx_end = events.connection_idx_end<IS_INTRON>();
    // how we convert that back to underlying junction/intron
    auto connection_at = [&events](size_t connection_idx) -> const auto& {
      return events.connection_at<IS_INTRON>(connection_idx); };
    // get regions from SJ we will be comparing to
    const auto sj_regions = *(sj.regions());
    // how we determine if we no longer have match in coordinates
    using RegionIntervalBeforeT = std::conditional_t<
      IS_INTRON,
      IntervalPrecedesT,
      std::less<decltype(sj_regions[0].coordinates)>>;

    // we will temporarily accumulate the total number of positions
    // contributing to each event connection. Afterwards, we will normalize to
    // total number of bins (junction rate)
    std::vector<junction_pos_t> weights(
        connection_idx_end - connection_idx_begin, 0);

    // loop over connection_idx. We determine the first contig and match that
    // in sj, keeping sj, connection_idx in sync/accumulating
    // weights/reads/bootstraps
    for (auto idx_it = connection_idx_begin; idx_it != connection_idx_end;) {
      const KnownContig& dst_contig = connection_at(*idx_it).contig();
      const auto opt_sj_contig_idx
        = sj_regions.parents()->safe_idx(dst_contig.get());
      if (opt_sj_contig_idx.has_value()) {
        const size_t& sj_contig_idx = *opt_sj_contig_idx;
        // loop over sj regions for this contig
        for (auto sj_it = sj_regions.begin_parent(sj_contig_idx);
            sj_it != sj_regions.end_parent(sj_contig_idx); ++sj_it) {
          // get first connection idx not behind sj_it
          idx_it = std::find_if(idx_it, connection_idx_end,
              [&connection_at, &sj_coordinates = sj_it->coordinates,
              &dst_contig](size_t connection_idx) {
              const auto& x = connection_at(connection_idx);
              return (x.contig() != dst_contig
                  || !RegionIntervalBeforeT{}(
                    x.coordinates, sj_coordinates)); });
          // if this is the end or a different contig, go to next contig
          if (idx_it == connection_idx_end
              || connection_at(*idx_it).contig() != dst_contig) { break; }
          // otherwise, are there connections that overlap?
          // if there are, we want to calculate the following coverage
          // information one and add to all matching connections
          struct _SJCoverageInfo {
            junction_pos_t weight_, num_stacks_;
            real_t numreads_, numbins_;
            std::vector<real_t> bootstraps_;
            _SJCoverageInfo(size_t num_bootstraps)
                : bootstraps_(num_bootstraps) { }
          };
          std::optional<_SJCoverageInfo> opt_coverage_info;
          for (auto overlap_it = idx_it;
              (overlap_it != connection_idx_end
               && connection_at(*overlap_it).contig() == dst_contig
               && !RegionIntervalBeforeT{}(
                 sj_it->coordinates, connection_at(*overlap_it).coordinates));
              ++overlap_it) {
            // potential match
            const auto& connection = connection_at(*overlap_it);
            // ignore non-matching strands
            if (sj_it->strand != GeneStrandness::AMBIGUOUS
                && sj_it->strand != connection.strand()) { continue; }
            // if first time, set values for opt_coverage_info
            if (!opt_coverage_info.has_value()) {
              opt_coverage_info = _SJCoverageInfo{bootstraps.num_bootstraps()};
              const auto sj_idx = sj_it - sj_regions.begin();
              opt_coverage_info->weight_ = sj.weight(sj_idx);
              opt_coverage_info->num_stacks_
                = sj.numstacks(sj_idx, pvalue_threshold);
              opt_coverage_info->numreads_
                = sj.numreads(sj_idx, opt_coverage_info->num_stacks_);
              opt_coverage_info->numbins_
                = sj.numbins_nonzero(sj_idx) - opt_coverage_info->num_stacks_;
              sj.bootstrap_coverage(sj_idx, opt_coverage_info->num_stacks_,
                  generator,
                  opt_coverage_info->bootstraps_.begin(),
                  opt_coverage_info->bootstraps_.end());
            }
            // add coverage_info into connection details
            const _SJCoverageInfo& coverage_info = *opt_coverage_info;
            weights[overlap_it - connection_idx_begin] += coverage_info.weight_;
            summaries[*overlap_it].numreads += coverage_info.numreads_;
            summaries[*overlap_it].numbins += coverage_info.numbins_;
            for (size_t j = 0; j < bootstraps.num_bootstraps(); ++j) {
              bootstraps(*overlap_it, j) += coverage_info.bootstraps_[j];
            }  // fill in all bootstraps
          }  // loop over overlapping event connections for current sj region
        }  // loop over sj regions for the contig
      }  // only iterate if both events and connections have contig
      // get to next idx_it that is on a different contig
      idx_it = std::find_if(idx_it, connection_idx_end,
          [&connection_at, &dst_contig](size_t connection_idx) {
          return connection_at(connection_idx).contig() != dst_contig; });
    }  // done looping over event connections

    // normalize by weights
    const real_t total_bins_f = static_cast<real_t>(sj.total_bins());
    for (auto it = connection_idx_begin; it != connection_idx_end; ++it) {
      const auto& weight = weights[it - connection_idx_begin];
      if (weight == 0) { continue; }
      const real_t scale = total_bins_f / weight;
      summaries[*it].numreads *= scale;
      summaries[*it].numbins *= scale;  // TODO(jaicher): should just be average
      std::transform(
          bootstraps.begin_connection(*it), bootstraps.end_connection(*it),
          bootstraps.begin_connection(*it),
          [&scale](const real_t& x) { return x * scale; });
    }  // done normalizing by weights
    return;
  }

 public:
  static EventsCoverage FromSJ(
      const std::shared_ptr<Events>& events,
      const SJJunctionsPositions& sj_junctions,
      const SJIntronsBins& sj_introns,
      size_t num_bootstraps,
      rng_t& generator,
      real_t pvalue_threshold) {
    if (events == nullptr) {
      throw std::runtime_error("EventsCoverage::FromSJ given null events");
    }
    std::vector<CoverageSummary> summaries(events->num_connections());
    CoverageBootstraps bootstraps{events->num_connections(), num_bootstraps};
    merge_bins(*events, sj_junctions,
        generator, pvalue_threshold, summaries, bootstraps);
    merge_bins(*events, sj_introns,
        generator, pvalue_threshold, summaries, bootstraps);
    return EventsCoverage{events, std::move(summaries), std::move(bootstraps)};
  }
};
}  // namespace majiq

#endif  // MAJIQ_EVENTSCOVERAGE_HPP
