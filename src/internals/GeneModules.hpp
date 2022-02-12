/**
 * GeneModules.hpp
 *
 * Modules for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_GENEMODULES_HPP
#define MAJIQ_GENEMODULES_HPP

#include <algorithm>
#include <map>
#include <memory>
#include <vector>

#include "Exons.hpp"
#include "ExonConnections.hpp"
#include "GeneIntrons.hpp"
#include "GeneJunctions.hpp"
#include "GeneRegion.hpp"
#include "Interval.hpp"
#include "Regions.hpp"
#include "SpliceGraphMask.hpp"

namespace majiq {

struct GeneModuleIndexes {
 public:
  size_t start_exon_idx;
  size_t end_exon_idx;
  size_t start_intron_idx;
  size_t num_introns;
  size_t start_junction_idx;
  size_t num_junctions;

  // constructors
  GeneModuleIndexes(size_t _start_exon_idx, size_t _end_exon_idx,
      size_t _start_intron_idx, size_t _num_introns,
      size_t _start_junction_idx, size_t _num_junctions)
    : start_exon_idx{_start_exon_idx},
      end_exon_idx{_end_exon_idx},
      start_intron_idx{_start_intron_idx},
      num_introns{_num_introns},
      start_junction_idx{_start_junction_idx},
      num_junctions{_num_junctions} {
    if (start_exon_idx > end_exon_idx) {
      throw std::invalid_argument("Module end exon cannot precede start");
    }
  }
  GeneModuleIndexes() : GeneModuleIndexes{0, 0, 0, 0, 0, 0} { }
  GeneModuleIndexes(const GeneModuleIndexes&) = default;
  GeneModuleIndexes(GeneModuleIndexes&&) = default;
  GeneModuleIndexes& operator=(const GeneModuleIndexes&) = default;
  GeneModuleIndexes& operator=(GeneModuleIndexes&&) = default;
};  // struct GeneModuleIndexes

struct GeneModule
    : public detail::GeneRegion<OpenInterval, GeneModuleIndexes> {
 public:
  using BaseT = detail::GeneRegion<OpenInterval, GeneModuleIndexes>;

  // constructors
  GeneModule(KnownGene gene, OpenInterval coordinates,
      GeneModuleIndexes indexes)
    : BaseT{gene, coordinates, indexes} { }
  GeneModule(KnownGene gene, OpenInterval coordinates)
    : GeneModule{gene, coordinates, GeneModuleIndexes{}} { }
  GeneModule() : GeneModule{KnownGene{}, OpenInterval{}} {}
  GeneModule(const GeneModule&) = default;
  GeneModule(GeneModule&&) = default;
  GeneModule& operator=(const GeneModule&) = default;
  GeneModule& operator=(GeneModule&&) = default;

  const size_t& start_exon_idx() const noexcept {
    return this->data.start_exon_idx;
  }
  const size_t& end_exon_idx() const noexcept {
    return this->data.end_exon_idx;
  }
  const size_t& start_intron_idx() const noexcept {
    return this->data.start_intron_idx;
  }
  const size_t& num_introns() const noexcept { return this->data.num_introns; }
  const size_t& start_junction_idx() const noexcept {
    return this->data.start_junction_idx;
  }
  const size_t& num_junctions() const noexcept {
    return this->data.num_junctions;
  }
};  // struct GeneModule

// false: no overlap between gene modules
class GeneModules : public detail::Regions<GeneModule, false> {
  using BaseT = detail::Regions<GeneModule, false>;

 private:
  // indexes into modules for each intron/junction
  const std::vector<size_t> introns_module_idx_;
  const std::vector<size_t> junctions_module_idx_;
  // exon connections on which modules live
  const std::shared_ptr<ExonConnections> exon_connections_;
  // mask over introns and junctions
  const std::shared_ptr<SpliceGraphMask> mask_;

  // extract vector of GeneModule for BaseT constructor
  struct IdentifiedModules {
    std::vector<GeneModule> modules;
    // indexes into modules for each intron/junction
    std::vector<size_t> introns_module_idx;
    std::vector<size_t> junctions_module_idx;
  };
  static IdentifiedModules IdentifyModules(
      const ExonConnections& exon_connections,
      const SpliceGraphMask& mask) {
    // check that they match
    if (exon_connections.introns() != mask.introns()) {
      throw std::invalid_argument(
          "IdentifyModules connections and mask do not share introns");
    }
    if (exon_connections.junctions() != mask.junctions()) {
      throw std::invalid_argument(
          "IdentifyModules connections and mask do not share junctions");
    }
    const auto& introns = *exon_connections.introns();
    const auto& junctions = *exon_connections.junctions();
    // populate vector of modules
    std::vector<GeneModule> result;
    // track which module each intron and junction belongs to
    std::vector<size_t> intron_module(introns.size());
    std::vector<size_t> junction_module(junctions.size());
    // track current gene, next valid intron, next valid junction.
    // we accumulate introns/junctions that belong to gene, identifying when
    // we have moved into the next gene or module

    bool module_started = false;
    // current gene, coordinates, data we are working with
    KnownGene gene{0, introns.parents()};
    OpenInterval coordinates{};
    GeneModuleIndexes indexes{};

    // get first intron_idx that isn't masked (or end of introns)
    size_t intron_idx = 0;
    for (; intron_idx < introns.size()
        && !static_cast<bool>(mask.introns_values()[intron_idx]);
        ++intron_idx) {
      intron_module[intron_idx] = 0;
    }
    // get first junction_idx that isn't masked (or end of junctions)
    size_t junction_idx = 0;
    for (; junction_idx < junctions.size()
        && !static_cast<bool>(mask.junctions_values()[junction_idx]);
        ++junction_idx) {
      junction_module[junction_idx] = 0;
    }

    /**
     * add module if one has already been started
     */
    auto AddModule
        = [&result, &module_started, &gene, &coordinates, &indexes,
          &intron_idx, &junction_idx, &intron_module, &junction_module]() {
      if (module_started) { result.emplace_back(gene, coordinates, indexes); }
      // update index to something reasonable.
      // most important quality is num_introns, num_junctions = 0.
      // using intron_idx, junction_idx keeps sequence of these monotonic.
      // this can help identify which module an intron or junction belongs to.
      indexes.start_exon_idx = indexes.end_exon_idx;
      indexes.start_intron_idx = intron_idx;
      indexes.num_introns = 0;
      indexes.start_junction_idx = junction_idx;
      indexes.num_junctions = 0;
      // set module_started flag to false
      module_started = false;
      // update current intron/junction module index to current module, if any
      if (intron_idx < intron_module.size()) {
        intron_module[intron_idx] = result.size();
      }
      if (junction_idx < junction_module.size()) {
        junction_module[junction_idx] = result.size();
      }
      return;
    };
    /**
     * Given intron, update coordinates (and maybe module) and intron_idx
     */
    auto TakeIntron
        = [&module_started, &AddModule, &coordinates, &introns, &indexes,
          &intron_idx, &introns_mask = mask.introns_values(),
          &result, &intron_module]() {
      // get intron coordinates as open interval (expand it by +/- 1)
      auto intron_open = introns[intron_idx].coordinates.AsOpen();
      // does the current intron occur after current module coordinates?
      if (!module_started || intron_open.start >= coordinates.end) {
        AddModule();
        coordinates = intron_open;
        indexes.start_intron_idx = intron_idx;
        indexes.start_exon_idx = introns[intron_idx].start_exon_idx();
        module_started = true;
      } else {
        // module is already started, and start of intron is in current module
        coordinates.end = std::max(coordinates.end, intron_open.end);
      }
      ++indexes.num_introns;  // increment count of valid (unmasked) introns
      indexes.end_exon_idx = std::max(
          indexes.end_exon_idx, introns[intron_idx].end_exon_idx());
      // update intron_idx to next valid intron
      for (intron_idx = 1 + intron_idx;
          intron_idx < introns.size()
          && !static_cast<bool>(introns_mask[intron_idx]);
          ++intron_idx) {
        intron_module[intron_idx] = result.size();
      }
      if (intron_idx < introns.size()) {
        intron_module[intron_idx] = result.size();
      }
      return;
    };
    /**
     * Given junction, update coordinates (and maybe module) and junction_idx
     */
    auto TakeJunction
        = [&module_started, &AddModule, &coordinates, &junctions, &indexes,
          &junction_idx, &junctions_mask = mask.junctions_values(),
          &result, &junction_module]() {
      // get junction coordinates
      const auto& junction_open = junctions[junction_idx].coordinates;
      // does the current junction occur after current module coordinates?
      if (!module_started || junction_open.start >= coordinates.end) {
        AddModule();
        coordinates = junction_open;
        indexes.start_junction_idx = junction_idx;
        indexes.start_exon_idx = junctions[junction_idx].start_exon_idx();
        module_started = true;
      } else {
        // module is already started, and start of junction is in current module
        coordinates.end = std::max(coordinates.end, junction_open.end);
      }
      ++indexes.num_junctions;  // increment count of valid (unmasked) junctions
      indexes.end_exon_idx = std::max(
          indexes.end_exon_idx, junctions[junction_idx].end_exon_idx());
      // update junction_idx to next valid junction
      for (junction_idx = 1 + junction_idx;
          junction_idx < junctions.size()
          && !static_cast<bool>(junctions_mask[junction_idx]);
          ++junction_idx) {
        junction_module[junction_idx] = result.size();
      }
      if (junction_idx < junctions.size()) {
        junction_module[junction_idx] = result.size();
      }
      return;
    };

    while (intron_idx < introns.size() || junction_idx < junctions.size()) {
      // 4 cases: both valid, only one valid, neither valid
      bool introns_valid = intron_idx < introns.size()
        && introns[intron_idx].gene == gene;
      bool junctions_valid = junction_idx < junctions.size()
        && junctions[junction_idx].gene == gene;
      if (introns_valid && junctions_valid) {
        // take first intron/junction (when treating introns as open intervals)
        if (junctions[junction_idx].coordinates
            < introns[intron_idx].coordinates.AsOpen()) {
          TakeJunction();
        } else {
          TakeIntron();
        }
      } else if (introns_valid) {
        TakeIntron();
      } else if (junctions_valid) {
        TakeJunction();
      } else {
        // add module, if appropriate
        AddModule();
        // update gene, if appropriate
        if (intron_idx >= introns.size()) {
          gene.idx_ = junctions[junction_idx].gene.idx_;
        } else if (junction_idx >= junctions.size()) {
          gene.idx_ = introns[intron_idx].gene.idx_;
        } else {
          gene.idx_ = std::min(junctions[junction_idx].gene.idx_,
              introns[intron_idx].gene.idx_);
        }  // update gene idx to same as next intron or junction
      }  // handle 4 cases to take next intron/junction or gene
    }  // iterate over all valid introns and junctions
    AddModule();  // add last module (if one was started)

    return IdentifiedModules{
        std::move(result),
        std::move(intron_module),
        std::move(junction_module)};
  }

 public:
  const std::vector<size_t>& introns_module_idx() const {
    return introns_module_idx_;
  }
  const std::vector<size_t>& junctions_module_idx() const {
    return junctions_module_idx_;
  }
  const std::shared_ptr<ExonConnections>& exon_connections() const {
    return exon_connections_;
  }
  const std::shared_ptr<SpliceGraphMask>& mask() const { return mask_; }

  bool event_is_masked(const Event& event) const {
    for (auto jit = exon_connections_->begin_junctions_for(event);
        jit != exon_connections_->end_junctions_for(event); ++jit) {
      if (mask_->junctions_values()[*jit]) {
        // there was a junction
        return false;
      }
    }
    // try the intron, if any
    for (auto iit = exon_connections_->begin_introns_for(event);
        iit != exon_connections_->end_introns_for(event); ++iit) {
      if (mask_->introns_values()[*iit]) {
        // there was an intron
        return false;
      }
    }
    // this event is entirely masked out
    return true;
  }

  size_t event_module_index(const Event& event) const {
    // infer module index of event by finding first junction in event that is
    // not masked out
    size_t last_module_idx{0};
    for (auto jit = exon_connections_->begin_junctions_for(event);
        jit != exon_connections_->end_junctions_for(event); ++jit) {
      last_module_idx = junctions_module_idx_[*jit];
      if (mask_->junctions_values()[*jit]) {
        // this junction wasn't masked out
        return last_module_idx;
      }
    }
    // try the intron, if any
    for (auto iit = exon_connections_->begin_introns_for(event);
        iit != exon_connections_->end_introns_for(event); ++iit) {
      last_module_idx = introns_module_idx_[*iit];
      if (mask_->introns_values()[*iit]) {
        // this intron wasn't masked out
        return last_module_idx;
      }
    }
    // this event is entirely masked out, so use the module the last visited
    // connection ostensibly maps to
    return last_module_idx;
  }

  // constructor
 private:
  GeneModules(
      IdentifiedModules&& identified,
      const std::shared_ptr<ExonConnections> exon_connections,
      const std::shared_ptr<SpliceGraphMask> mask)
    : BaseT{exon_connections->exons()->parents(),
      std::move(identified.modules)},
      introns_module_idx_{std::move(identified.introns_module_idx)},
      junctions_module_idx_{std::move(identified.junctions_module_idx)},
      exon_connections_{exon_connections},
      mask_{mask} {
    if (introns_module_idx_.size() != exon_connections_->introns()->size()) {
      throw std::invalid_argument(
          "Invalid size for index into modules for introns");
    }
    if (junctions_module_idx_.size() != exon_connections_->junctions()->size()) {
      throw std::invalid_argument(
          "Invalid size for index into modules for junctions");
    }
  }
 public:
  GeneModules(
      const std::shared_ptr<ExonConnections> exon_connections,
      const std::shared_ptr<SpliceGraphMask> mask)
    : GeneModules{
        std::move(IdentifyModules(*exon_connections, *mask)),
        exon_connections,
        mask} { }
};  // class GeneModules

}  // namespace majiq

#endif  // MAJIQ_GENEMODULES_HPP
