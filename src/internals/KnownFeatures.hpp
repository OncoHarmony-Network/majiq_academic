/**
 * KnownFeatures.hpp
 *
 * Parent template class for parent features (i.e. contigs, genes) for
 * referring to them by indexes instead of storing multiple copies (flypaper
 * pattern)
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_KNOWNFEATURES_HPP
#define MAJIQ_KNOWNFEATURES_HPP

#include <memory>
#include <utility>
#include <map>
#include <tuple>
#include <optional>
#include <functional>
#include <boost/functional/hash.hpp>


namespace majiq {
namespace detail {
template <typename ContainerT> class KnownFeatures;
template <typename KnownFeaturesT> struct KnownFeature {
  using FeatureT = typename KnownFeaturesT::FeatureT;
  size_t idx_;
  std::shared_ptr<KnownFeaturesT> ptr_;

  // get underlying feature
  const FeatureT& get() const { return ptr_->get(idx_); }
  KnownFeature(size_t idx, std::shared_ptr<KnownFeaturesT> ptr)
      : idx_{idx}, ptr_{ptr} { }
  KnownFeature() = default;
  KnownFeature(const KnownFeature&) = default;
  KnownFeature(KnownFeature&&) = default;
  KnownFeature& operator=(const KnownFeature&) = default;
  KnownFeature& operator=(KnownFeature&&) = default;

  friend inline bool operator<(const KnownFeature& x, const KnownFeature& y) noexcept {
    return std::tie(x.ptr_, x.idx_) < std::tie(y.ptr_, y.idx_);
  }
  friend inline bool operator==(const KnownFeature& x, const KnownFeature& y) noexcept {
    return x.idx_ == y.idx_ && x.ptr_ == y.ptr_;
  }
  friend inline bool operator!=(const KnownFeature& x, const KnownFeature& y) noexcept {
    return !(x == y);
  }
};

template <typename ContainerT>
class KnownFeatures {
 public:
  using FeatureT = typename ContainerT::value_type;
  using KeyT = decltype(std::declval<FeatureT>().unique_key());

 protected:
  std::map<KeyT, size_t> idx_map_;
  ContainerT features_;

 public:
  // retrieve underlying features
  const FeatureT& get(size_t idx) const { return features_[idx]; }

  // NOTE: need pointer to *derived* classes, so add appropriately to subclasses
  // // get "known" reference to these elements
  // using Known = KnownFeature<KnownFeatures>;
  // const Known operator[](size_t idx) { return Known{idx, this->shared_from_this()}; }

  // work with data that's there
  size_t size() const { return features_.size(); }
  size_t count(const KeyT& key) const { return idx_map_.count(key); }
  size_t count(const FeatureT& x) const { return count(x.unique_key()); }
  size_t get_idx(const KeyT& key) const { return idx_map_.at(key); }
  size_t get_idx(const FeatureT& x) const { return get_idx(x.unique_key()); }
  std::optional<size_t> safe_idx(const KeyT& key) {
    std::optional<size_t> result;
    auto match = idx_map_.find(key);
    if (match != idx_map_.end()) {
      result = match->second;
    }
    return result;
  }
  std::optional<size_t> safe_idx(const FeatureT& x) { return safe_idx(x.unique_key()); }

  explicit KnownFeatures(ContainerT&& features)
      : features_{std::move(features)} {
    for (size_t idx = 0; idx < features_.size(); ++idx) {
      idx_map_[features_[idx].unique_key()] = idx;
    }
  }
  KnownFeatures() = default;
  KnownFeatures(const KnownFeatures&) = default;
  KnownFeatures(KnownFeatures&&) = default;
  KnownFeatures& operator=(const KnownFeatures&) = delete;
  KnownFeatures& operator=(KnownFeatures&&) = delete;

  friend inline bool operator==(
      const KnownFeatures& x, const KnownFeatures& y) noexcept {
    return x.features_ == y.features_;
  }
};
template <typename T>
inline std::size_t hash_value(const KnownFeature<T>& x) noexcept {
  std::size_t result = std::hash<size_t>{}(x.idx_);
  boost::hash_combine(result, x.ptr_);
  return result;
}
}  // namespace detail
}  // namespace majiq

namespace std {
template <typename T> struct hash<majiq::detail::KnownFeature<T>> {
  std::size_t operator()(const majiq::detail::KnownFeature<T>& x) const noexcept {
    return majiq::detail::hash_value(x);
  }
};
}  // namespace std

#endif  // MAJIQ_KNOWNFEATURES_HPP
