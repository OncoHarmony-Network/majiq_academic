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

  // constructors
  Connection(bool _denovo, bool _passed_build, bool _simplified)
      : denovo{_denovo}, passed_build{_passed_build}, simplified{_simplified} {
  }
  Connection() : Connection{false, false, false} { }
  Connection(const Connection& x) = default;
  Connection(Connection&& x) = default;
  Connection& operator=(const Connection& x) = default;
  Connection& operator=(Connection&& x) = default;
};
}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_CONNECTION_HPP
