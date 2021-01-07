/**
 * ContiguousSet.hpp
 *
 * Set implementation that is lazily added to a sorted vector of unique using a
 * set (or unordered set) to hold elements.
 * Does not permit removal of elements.
 *
 * Copypright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_UTILS_CONTIGUOUS_SET_HPP
#define MAJIQ_UTILS_CONTIGUOUS_SET_HPP


#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <functional>


namespace majiq {
namespace utils {
template <class Key, class Compare = std::less<Key>>
struct DefaultSetUnion {
  using value_type = Key;
  using key_compare = Compare;
  template <class InputIt1, class InputIt2, class OutputIt>
  OutputIt operator()(
      InputIt1 first1, InputIt1 last1,
      InputIt2 first2, InputIt2 last2,
      OutputIt d_first) {
    static_assert(
        std::is_same<
          value_type,
          typename std::iterator_traits<InputIt1>::value_type>::value,
        "InputIt1 must iterate over type Key in DefautSetUnion");
    static_assert(
        std::is_same<
          value_type,
          typename std::iterator_traits<InputIt2>::value_type>::value,
        "InputIt2 must iterate over type Key in DefautSetUnion");
    return std::set_union(first1, last1, first2, last2, d_first, Compare{});
  }
};

template <
  class Key,
  class Compare = std::less<Key>,
  class SetUnion = DefaultSetUnion<Key, Compare>>
class contiguous_set {
 public:
  using value_type = Key;
  using key_compare = Compare;
  // check that set union matches
  static_assert(
      std::is_same<value_type, typename SetUnion::value_type>::value,
      "SetUnion must have matching Key in contiguous_set");
  static_assert(
      std::is_same<key_compare, typename SetUnion::key_compare>::value,
      "SetUnion must have matching Compare in contiguous set");

  // containers
  using vecT = std::vector<value_type>;
  using vecConstT = std::vector<const value_type>;
  // what we could use for the unsorted elements?
  using _ordered_setT = std::set<value_type, key_compare>;
  using _unordered_setT = std::unordered_set<value_type>;
  // what we are using (either one of them)
  using setT = _unordered_setT;

  // iterators
  using const_iterator = typename vecT::const_iterator;

 private:
  // static functions used by constructors, etc
  template <bool other_sorted = false, class otherT>
  static vecT _combine_sorted(vecT sorted, otherT other) {
    // otherT must be a container of value_type
    static_assert(
        std::is_same<value_type, typename otherT::value_type>::value,
        "otherT must have same type as contiguous_set");
    if (other.size() == 0) {
      return sorted;
    }
    vecT out;
    out.reserve(sorted.size() + other.size());
    if constexpr(other_sorted || std::is_same<otherT, _ordered_setT>::value) {
      // other is aready sorted
      SetUnion{}(
          sorted.begin(), sorted.end(), other.begin(), other.end(),
          std::back_inserter(out));
    } else {
      vecT other_vec{other.begin(), other.end()};
      std::sort(other_vec.begin(), other_vec.end(), Compare{});
      SetUnion{}(
          sorted.begin(), sorted.end(), other_vec.begin(), other_vec.end(),
          std::back_inserter(out));
    }
    return out;
  }

 protected:  // THE UNDERLYING CONTAINERS
  mutable vecT sorted_vec_;
  mutable setT lazy_set_;  // temporary elements waiting to be merged in

 public:
  size_t size() const noexcept {
    return sorted_vec_.size() + lazy_set_.size();
  }
  bool empty() const noexcept { return size() == 0; }

 protected:  // MAKE AND ACCESS CONTIGUOUS DATA
  // check/make sure it's sorted
  bool is_contiguous() const { return lazy_set_.size() == 0; }
  void make_contiguous() const {
    if (!is_contiguous()) {
      sorted_vec_ = _combine_sorted(sorted_vec_, lazy_set_);
      lazy_set_.clear();
    }
    return;
  }

 public:
  // access elements by index
  const value_type& operator[](size_t idx) const {
    make_contiguous();
    return sorted_vec_[idx];
  }
  // constant iterators over contiguous values
  const_iterator begin() const {
    make_contiguous();
    return sorted_vec_.cbegin();
  }
  const_iterator end() const {
    make_contiguous();
    return sorted_vec_.cend();
  }
  // view into vector for pybind11
  const vecT& data() {
    make_contiguous();
    return sorted_vec_;
  }
  const vecConstT& data() const {
    make_contiguous();
    return sorted_vec_;
  }

 protected:  // MEMBERSHIP IN SET
  // check containment in either set
  bool _sorted_contains(const value_type& key) const {
    auto eq_range = std::equal_range(
        sorted_vec_.begin(), sorted_vec_.end(), key, key_compare{});
    return eq_range.first != eq_range.second;
  }
  bool _lazy_contains(const value_type& key) const {
    return lazy_set_.count(key) > 0;
  }

 public:
  size_t count(const value_type& key) const {
    return (_sorted_contains(key) || _set_contains(key)) ? 1 : 0;
  }
  bool insert(const value_type& key) {
    if (sorted_vec_.size() == 0 || Compare{}(sorted_vec_.back(), key)) {
      // key is not in this, and it can be added to contiguous portion
      sorted_vec_.push_back(key);
      return true;  // we added it
    } else if (_sorted_contains(key)) {
      return false;  // didn't add because it was in contiguous portion
    } else {
      // insert into set, determine if it was already there (.second)
      return lazy_set_.insert(key).second;
    }
  }

  // comparison
  friend inline bool operator==(
      const contiguous_set& x, const contiguous_set& y) {
    x.make_contiguous();
    y.make_contiguous();
    return x.sorted_vec_ == y.sorted_vec_;
  }

  // constructors
  contiguous_set() = default;
  contiguous_set(const contiguous_set& x) = default;
  contiguous_set(contiguous_set&& x) = default;
  contiguous_set& operator=(const contiguous_set& x) = default;
  contiguous_set& operator=(contiguous_set&& x) = default;
};
}  // namespace utils
}  // namespace majiq

#endif  // MAJIQ_UTILS_CONTIGUOUS_SET_HPP
