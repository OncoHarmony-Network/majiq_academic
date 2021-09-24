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

namespace MajiqGufuncs {
namespace detail {
/**
 * Get x[idx] given type and stride of array
 */
template <typename T>
inline T& get_value(char* x_ptr, npy_intp idx, npy_intp stride) {
  return *reinterpret_cast<T*>(x_ptr + stride * idx);
}
template <typename T>
inline T& get_value(char* x_ptr) {
  return *reinterpret_cast<T*>(x_ptr);
}


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
  npy_intp stride_;  // XXX: don't let this be 0.

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
    return get_value<value_type>(ptr_, n, stride_);
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
    if (stride == 0) {
      // make it possible to reach end iterator
      stride = sizeof(value_type);
    }
    return CoreIt<T>{ptr, stride};
  }
  static CoreIt end(const CoreIt& x, npy_intp n) {
    return x + n;
  }
};

}  // namespace detail
}  // namespace MajiqGufuncs

#endif  // MAJIQGUFUNCS_HELPERS_HPP
