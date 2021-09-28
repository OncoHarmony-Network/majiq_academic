/**
 * helpers.hpp
 *
 * MajiqGufuncs::detail helper functions
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQGUFUNCS_HELPERS_HPP
#define MAJIQGUFUNCS_HELPERS_HPP

#include <numpy/ndarraytypes.h>
#include <iterator>

namespace MajiqGufuncs {
namespace detail {

/**
 * iterator over 1d core dimensions in gufunc inner loops
 */
template <typename T>
class CoreIt {
 public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = T;
  using difference_type = npy_intp;
  using pointer = value_type*;
  using reference = value_type&;

 private:
  char* ptr_;
  npy_intp stride_;

  // private constructor because we can't allow zero stride
  CoreIt(char* ptr, npy_intp stride)
    : ptr_{ptr}, stride_{stride} { }

 public:
  // default constructors/assignment
  CoreIt() : CoreIt{nullptr, sizeof(value_type)} { }
  CoreIt(const CoreIt&) = default;
  CoreIt(CoreIt&&) = default;
  CoreIt& operator=(const CoreIt&) = default;
  CoreIt& operator=(CoreIt&&) = default;

  // directly apply stride to object
  void apply_stride(npy_intp s) noexcept {
    ptr_ += s;
    return;
  }

  // get iterator with different stride
  CoreIt with_stride(npy_intp s) noexcept {
    return CoreIt{ptr_, s};
  }

  pointer operator->() const noexcept {
    return reinterpret_cast<pointer>(ptr_);
  }
  reference operator*() const noexcept {
    return *operator->();
  }
  CoreIt& operator++() noexcept {
    ptr_ += stride_;
    return *this;
  }
  CoreIt operator++(int) noexcept {
    CoreIt temp = *this;
    ++(*this);
    return temp;
  }
  CoreIt& operator--() noexcept {
    ptr_ -= stride_;
    return *this;
  }
  CoreIt operator--(int) noexcept {
    CoreIt temp = *this;
    --(*this);
    return temp;
  }
  CoreIt& operator+=(difference_type n) noexcept {
    ptr_ += n * stride_;
    return *this;
  }
  CoreIt& operator-=(difference_type n) noexcept {
    ptr_ -= n * stride_;
    return *this;
  }
  friend inline CoreIt operator+(
      const CoreIt& x, difference_type n) noexcept {
    CoreIt temp = x;
    temp += n;
    return temp;
  }
  friend inline CoreIt operator+(
      difference_type n, const CoreIt& x) noexcept {
    return x + n;
  }
  friend inline CoreIt operator-(
      const CoreIt& x, difference_type n) noexcept {
    CoreIt temp = x;
    temp -= n;
    return temp;
  }
  friend inline CoreIt operator-(
      difference_type n, const CoreIt& x) noexcept {
    return x - n;
  }
  // x + (x - y) == y
  friend inline difference_type operator-(
      const CoreIt& x, const CoreIt& y) noexcept {
    return (x.ptr_ - y.ptr_) / x.stride_;
  }
  reference operator[](difference_type n) const noexcept {
    return *reinterpret_cast<pointer>(ptr_ + stride_ * n);
  }

  friend inline bool operator==(
      const CoreIt& x, const CoreIt& y) noexcept {
    return x.ptr_ == y.ptr_;
  }
  friend inline bool operator!=(
      const CoreIt& x, const CoreIt& y) noexcept {
    return !(x == y);
  }
  friend inline bool operator<(
      const CoreIt& x, const CoreIt& y) noexcept {
    return (x.stride_ > 0) ? (x.ptr_ < y.ptr_) : (y.ptr_ < x.ptr_);
  }
  friend inline bool operator>(
      const CoreIt& x, const CoreIt& y) noexcept {
    return y < x;
  }
  friend inline bool operator>=(
      const CoreIt& x, const CoreIt& y) noexcept {
    return !(x < y);
  }
  friend inline bool operator<=(
      const CoreIt& x, const CoreIt& y) noexcept {
    return !(x > y);
  }

  static CoreIt begin(char* ptr, npy_intp stride) {
    return CoreIt<T>{ptr, stride};
  }
  // XXX: this pattern does not work if n == 1 because numpy will set stride to
  // 0 for broadcasting
  static CoreIt end(const CoreIt& x, npy_intp n) {
    return x + n;
  }

  // how many values to fill starting at current position?
  void fill(npy_intp n, const T& value) {
    if (stride_ != 0) {
      std::fill(*this, *this + n, value);
    } else if (n != 0) {
      // stride is zero and we wanted some values to be set
      operator*() = value;
    }
  }
};

}  // namespace detail
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_HELPERS_HPP
