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
#include <vector>
#include <boost/functional/hash.hpp>


namespace majiq {
namespace detail {
template <typename ContainerT> class KnownFeatures;

// NOTE: must put DerivedKnownT = name of class inheriting KnownFeature
template <typename KnownFeaturesT, typename DerivedKnownT>
struct KnownFeature {
  using FeatureT = decltype(std::declval<KnownFeaturesT>().get(0));
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

  friend inline bool operator<(
      const KnownFeature& x, const KnownFeature& y) noexcept {
    return x.idx_ < y.idx_;
  }
  friend inline bool operator==(
      const KnownFeature& x, const KnownFeature& y) noexcept {
    return x.idx_ == y.idx_;
  }

  // NOTE: hack to avoid reimplementing this for derived classes
  // make this a random_access_iterator of itself (over KnownFeaturesT
  using iterator_category = std::random_access_iterator_tag;
  using value_type = DerivedKnownT;  // itself!
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  reference operator*() noexcept {
    return *static_cast<DerivedKnownT*>(this);  // itself!
  }
  DerivedKnownT& operator++() noexcept {
    ++idx_;
    return *static_cast<DerivedKnownT*>(this);
  }
  DerivedKnownT operator++(int) noexcept {
    DerivedKnownT old = *static_cast<DerivedKnownT*>(this);
    operator++();
    return old;
  }
  DerivedKnownT& operator--() noexcept {
    --idx_;
    return *static_cast<DerivedKnownT*>(this);
  }
  DerivedKnownT operator--(int) noexcept {
    DerivedKnownT old = *static_cast<DerivedKnownT*>(this);
    operator--();
    return old;
  }
  DerivedKnownT& operator+=(difference_type n) noexcept {
    idx_ += n;
    return *static_cast<DerivedKnownT*>(this);
  }
  friend DerivedKnownT operator+(
      const DerivedKnownT& lhs, difference_type n) noexcept {
    return DerivedKnownT{lhs.idx_ + n, lhs.ptr_};
  }

  // derived
  friend inline bool operator>(
      const KnownFeature& x, const KnownFeature& y) noexcept {
    return y < x;
  }
  friend inline bool operator>=(
      const KnownFeature& x, const KnownFeature& y) noexcept {
    return !(x < y);
  }
  friend inline bool operator<=(
      const KnownFeature& x, const KnownFeature& y) noexcept {
    return !(y < x);
  }
  friend inline bool operator!=(
      const KnownFeature& x, const KnownFeature& y) noexcept {
    return !(x == y);
  }
  DerivedKnownT& operator-=(difference_type n) noexcept {
    return (*static_cast<DerivedKnownT*>(this) += -n);
  }
  friend DerivedKnownT operator-(
      const DerivedKnownT& lhs, difference_type n) noexcept {
    return lhs + (-n);
  }
  friend difference_type operator-(
      const DerivedKnownT& lhs, const DerivedKnownT& rhs) noexcept {
    return lhs.idx_ - rhs.idx_;
  }
  reference operator[](difference_type n) {
    return *(*static_cast<DerivedKnownT*>(this) + n);
  }
};

template <typename ContainerT>
class KnownFeatures {
 public:
  using FeatureT = std::remove_const_t<std::remove_reference_t<
    decltype(std::declval<ContainerT>()[0])>>;
  using KeyT = decltype(std::declval<FeatureT>().unique_key());

 protected:
  std::map<KeyT, size_t> idx_map_;
  ContainerT features_;

 public:
  // retrieve underlying features
  const FeatureT& get(size_t idx) const { return features_[idx]; }

  // work with data that's there
  size_t size() const { return features_.size(); }
  size_t count(const KeyT& key) const { return idx_map_.count(key); }
  size_t count(const FeatureT& x) const { return count(x.unique_key()); }
  size_t get_idx(const KeyT& key) const { return idx_map_.at(key); }
  size_t get_idx(const FeatureT& x) const { return get_idx(x.unique_key()); }
  std::optional<size_t> safe_idx(const KeyT& key) const {
    std::optional<size_t> result;
    auto match = idx_map_.find(key);
    if (match != idx_map_.end()) {
      result = match->second;
    }
    return result;
  }
  std::optional<size_t> safe_idx(const FeatureT& x) const {
    return safe_idx(x.unique_key());
  }

  // move constructor for ContainerT
  template <typename CT,
           std::enable_if_t<std::is_same_v<CT, ContainerT>, bool> = true>
  explicit KnownFeatures(CT&& features) : features_{std::move(features)} {
    for (size_t idx = 0; idx < features_.size(); ++idx) {
      idx_map_[features_[idx].unique_key()] = idx;
    }
  }
  KnownFeatures() = default;
  KnownFeatures(const KnownFeatures&) = default;
  KnownFeatures(KnownFeatures&&) = default;
  KnownFeatures& operator=(const KnownFeatures&) = default;
  KnownFeatures& operator=(KnownFeatures&&) = default;

  friend inline bool operator==(
      const KnownFeatures& x, const KnownFeatures& y) noexcept {
    return x.features_ == y.features_;
  }
};
template <typename T, typename U>
inline std::size_t hash_value(const KnownFeature<T, U>& x) noexcept {
  std::size_t result = std::hash<size_t>{}(x.idx_);
  boost::hash_combine(result, x.ptr_);
  return result;
}
}  // namespace detail
}  // namespace majiq

namespace std {
template <typename T, typename U>
struct hash<majiq::detail::KnownFeature<T, U>> {
  std::size_t operator()(
      const majiq::detail::KnownFeature<T, U>& x) const noexcept {
    return majiq::detail::hash_value(x);
  }
};
}  // namespace std

#endif  // MAJIQ_KNOWNFEATURES_HPP
