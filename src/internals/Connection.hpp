/**
 * Connection.hpp
 *
 * Mixin class for junctions/introns, which can be denovo (or annotated),
 * passed (for LSVs), and simplified
 *
 * NOTE that we make passed_build and simplified mutable. This means that they
 * can be modified even in a const value (i.e. a const Connection can still
 * change these values). This helps us ensure that our junctions/introns/exons
 * maintain identity (const in Regions_) but able to be passed or simplified
 * inplace
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CONNECTION_HPP
#define MAJIQ_CONNECTION_HPP


namespace majiq {
namespace detail {
struct Connection {
 public:
  bool denovo;
  mutable bool passed_build;
  mutable bool simplified;
  mutable size_t start_exon_idx;
  mutable size_t end_exon_idx;

  // true if start/end exon idx match. (for use when connections made)
  const bool is_exitronic() const noexcept {
    return start_exon_idx == end_exon_idx;
  }

  // constructors
  Connection(bool _denovo, bool _passed_build, bool _simplified,
      size_t _start_exon_idx, size_t _end_exon_idx)
      : denovo{_denovo},
        passed_build{_passed_build},
        simplified{_simplified},
        start_exon_idx{_start_exon_idx},
        end_exon_idx{_end_exon_idx} { }
  Connection(const Connection& x, size_t _start_exon_idx, size_t _end_exon_idx)
      : Connection{x.denovo, x.passed_build, x.simplified,
          _start_exon_idx, _end_exon_idx} { }
  Connection(bool _denovo, bool _passed_build, bool _simplified)
      : Connection{_denovo, _passed_build, _simplified, 0, 0} { }
  Connection() : Connection{false, false, false} { }
  Connection(const Connection& x) = default;
  Connection(Connection&& x) = default;
  Connection& operator=(const Connection& x) = default;
  Connection& operator=(Connection&& x) = default;
};
}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_CONNECTION_HPP
