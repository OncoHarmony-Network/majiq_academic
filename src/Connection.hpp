/**
 * Connection.hpp
 *
 * Mixin class for junctions/introns, which can be denovo (or annotated),
 * passed (for LSVs), and simplified
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
  bool passed_build;
  bool simplified;

  // constructors
  Connection(bool _denovo, bool _passed_build, bool _simplified)
      : denovo{_denovo}, passed_build{_passed_build}, simplified{_simplified} {
  }
  Connection() : Connection{false, false, false} { }
  Connection(const Connection& x) = default;
  Connection(Connection&& x) = default;
  Connection& operator=(const Connection& x) = default;
  Connection& operator=(Connection&& x) = default;
  Connection& operator&=(const Connection& rhs) noexcept;
};
Connection operator&(const Connection& x, const Connection& y) noexcept {
  return Connection{
    x.denovo && y.denovo,
    x.passed_build || y.passed_build,
    x.simplified && y.simplified};
}
Connection& Connection::operator&=(const Connection& rhs) noexcept {
  denovo = denovo && rhs.denovo;
  passed_build = passed_build || rhs.passed_build;
  simplified = simplified && rhs.simplified;
  return *this;
}
}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_CONNECTION_HPP
