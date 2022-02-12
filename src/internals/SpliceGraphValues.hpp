/*
 * SpliceGraphValues.hpp
 *
 * Template base class for values of given type over gene introns and gene
 * junctions
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_SPLICEGRAPHVALUES_HPP
#define MAJIQ_SPLICEGRAPHVALUES_HPP

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "GeneIntrons.hpp"
#include "GeneJunctions.hpp"


namespace majiq {
namespace detail {

template <typename ValueT>
class SpliceGraphValues {
 public:
  using value_type = ValueT;

 private:
  const std::shared_ptr<GeneIntrons> introns_;
  const std::shared_ptr<GeneJunctions> junctions_;
  const std::vector<value_type> introns_values_;
  const std::vector<value_type> junctions_values_;

 public:
  SpliceGraphValues(
      const std::shared_ptr<GeneIntrons>& introns,
      const std::shared_ptr<GeneJunctions>& junctions,
      std::vector<value_type>&& introns_values,
      std::vector<value_type>&& junctions_values)
      : introns_{introns},
        junctions_{junctions},
        introns_values_{std::move(introns_values)},
        junctions_values_{std::move(junctions_values)} {
    if (introns_ == nullptr) {
      throw std::runtime_error("SpliceGraphValues given null introns");
    } else if (junctions_ == nullptr) {
      throw std::runtime_error("SpliceGraphValues given null junctions");
    } else if (introns_->size() != introns_values_.size()) {
      throw std::runtime_error(
          "SpliceGraphValues introns values do not match introns in size");
    } else if (junctions_->size() != junctions_values_.size()) {
      throw std::runtime_error(
          "SpliceGraphValues junctions values do not match junctions in size");
    }
  }
  SpliceGraphValues(const SpliceGraphValues&) = default;
  SpliceGraphValues(SpliceGraphValues&&) = default;
  SpliceGraphValues& operator=(const SpliceGraphValues&) = delete;
  SpliceGraphValues& operator=(SpliceGraphValues&&) = delete;

  const std::shared_ptr<GeneIntrons>& introns() const { return introns_; }
  const std::shared_ptr<GeneJunctions>& junctions() const { return junctions_; }
  const std::vector<value_type>& introns_values() const {
    return introns_values_;
  }
  const std::vector<value_type>& junctions_values() const {
    return junctions_values_;
  }
};  // class SpliceGraphValues

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_SPLICEGRAPHVALUES_HPP
