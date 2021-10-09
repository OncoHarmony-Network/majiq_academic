/**
 * RNGPool.hpp
 *
 * Threadsafe pool of random number generators
 *
 * Copyright 2021 <University of Pennsylvania>
 *
 * Author: Joseph K Aicher
 */

#ifndef MAJIQINCLUDE_RNGPOOL_HPP
#define MAJIQINCLUDE_RNGPOOL_HPP

#include <condition_variable>
#include <mutex>
#include <memory>
#include <vector>
#include <utility>

namespace MajiqInclude {

// inspired by https://stackoverflow.com/a/27837534
template <typename Generator>
class RNGPool {
 private:
  // how many resources have we made?
  int64_t n_;
  // what was the initial seed?
  int64_t base_seed_;
  // shared pointer, of pointer to self. weak_ptr of this held when resource
  // acquired by another thread, so when thread finished, can know if pool
  // still exists (return back), otherwise delete
  std::shared_ptr<RNGPool<Generator>*> this_ptr_;
  std::vector<std::unique_ptr<Generator>> pool_;
  std::mutex pool_mutex_;
  std::condition_variable pool_cv_;

 public:
  bool empty() const { return pool_.empty(); }
  int64_t n_available() const { return pool_.size(); }
  int64_t n() const { return n_; }

  /**
   * Put the generator into the pool
   */
  void add(std::unique_ptr<Generator> g) {
    std::lock_guard<std::mutex> lock{pool_mutex_};
    pool_.push_back(std::move(g));
    pool_cv_.notify_one();
  }

  /**
   * have at least new_n generators (including those acquired)
   */
  void resize(int64_t new_n) {
    // create as many generators as requested
    for (; n_ < new_n; ++n_) {
      std::unique_ptr<Generator> g{new Generator{}};
      g->seed(base_seed_ + n_);
      add(std::move(g));
    }
    return;
  }
  /**
   * reseed available generators
   */
  void seed(int64_t seed) {
    std::lock_guard<std::mutex> lock{pool_mutex_};
    // reseed generators currently in the pool
    base_seed_ = seed;
    for (int64_t i = 0; i < n_available(); ++i) {
      pool_[i]->seed(base_seed_ + i);
    }
    return;
  }

  /**
   * Initialize pool with specified number of generators, initial seed
   */
  RNGPool(int64_t n, int64_t seed)
    : n_{0},
      base_seed_{seed},
      this_ptr_{new RNGPool<Generator>*{this}} {
    resize(n);
    return;
  }
  RNGPool() : RNGPool{1, 20211008} { }
  RNGPool(const RNGPool&) = delete;
  RNGPool(RNGPool&&) = default;
  RNGPool& operator=(const RNGPool&) = delete;
  RNGPool& operator=(RNGPool&&) = default;

  // To acquire/return objects back to pool, we give a unique pointer with
  // custom deleter that returns the acquired generator back to the pool (if
  // available)
 private:
  class ExternalDeleter {
   private:
    std::weak_ptr<RNGPool<Generator>*> pool_;  // where to return object

   public:
    explicit ExternalDeleter(std::weak_ptr<RNGPool<Generator>*> pool)
      : pool_{pool} { }
    void operator()(Generator* ptr) {
      // return it back if possible
      if (auto pool_ptr = pool_.lock()) {
        try {
          (*pool_ptr)->add(std::unique_ptr<Generator>{ptr});
          return;
        } catch(...) {
          // go to end of function, where we delete ptr
        }
      }
      std::default_delete<Generator>{}(ptr);
    }
  };

 public:
  using ExternalPtrT = std::unique_ptr<Generator, ExternalDeleter>;
  ExternalPtrT acquire() {
    std::unique_lock<std::mutex> lock{pool_mutex_};
    pool_cv_.wait(lock, [this]() -> bool { return !empty(); });
    ExternalPtrT result{pool_.back().release(),
      ExternalDeleter{std::weak_ptr<RNGPool<Generator>*>{this_ptr_}}};
    pool_.pop_back();
    lock.unlock();
    return std::move(result);
  }
};

}  // namespace MajiqInclude

#endif  // MAJIQINCLUDE_RNGPOOL_HPP
