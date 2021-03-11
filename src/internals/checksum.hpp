/**
 * checksum.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_CHECKSUM_HPP
#define MAJIQ_CHECKSUM_HPP

#include <boost/crc.hpp>
#include <cstdint>


namespace majiq {
namespace detail {

using checksum_gen_t = boost::crc_32_type;
using checksum_t = uint32_t;

}  // namespace detail
}  // namespace majiq

#endif  // MAJIQ_CHECKSUM_HPP
