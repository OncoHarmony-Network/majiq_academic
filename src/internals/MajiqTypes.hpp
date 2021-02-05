/**
 * MaijqTypes.hpp
 *
 * base types for MAJIQ
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_TYPES_HPP
#define MAJIQ_TYPES_HPP

#include <cstdint>
#include <string>

namespace majiq {

using real_t = float;  // for numerics
using position_t = int64_t;  // for coordinates, difference in coordinates
using junction_pos_t = uint32_t;  // position for junctions
using junction_ct_t = uint32_t;  // counts of reads with splits

using seqid_t = std::string;  // type for contig ids
using geneid_t = std::string;
using genename_t = std::string;

enum class GeneStrandness : unsigned char {
  FORWARD = '+',
  REVERSE = '-',
  AMBIGUOUS = '.',
};

enum class ExperimentStrandness : unsigned char {
  FORWARD = 'F',  // read 1 is in forward direction, Salmon ISF
  REVERSE = 'R',  // read 2 is in forward direction, Salmon ISR
  NONE = 'N',  // could be either way, Salmon IU
};

}  // namespace majiq

#endif  // MAJIQ_TYPES_HPP
