/**
 * Alignments.hpp
 *
 * Interface for reading individual alignments from BAM file
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_BAM_ALIGNMENTS_HPP
#define MAJIQ_BAM_ALIGNMENTS_HPP

#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <sstream>
#include <stdexcept>
#include "CigarIntrons.hpp"
#include "CigarJunctions.hpp"
#include "../MajiqTypes.hpp"

namespace majiq {
namespace bam {
class AlignmentFile;

class AlignmentRecord {
 private:
  bam1_t *b_;  // holds current alignment information

 public:
  inline int32_t read_length() const { return b_->core.l_qseq; }
  /**
   * 0-based leftmost coordinate
   */
  inline position_t pos() const { return b_->core.pos; }

  /**
   * target id as defined by header (index into contig information)
   */
  inline int32_t tid() const { return b_->core.tid; }

  /**
   * query name of read, null-terminated
   */
  inline const char* query_name() const { return bam_get_qname(b_); }

  /**
   * read is primary alignment, mapped to reference
   */
  inline bool unique_mapped() const {
    // cannot have unmapped or secondary alignment flag set
    return !(b_->core.flag & (BAM_FUNMAP | BAM_FSECONDARY));
  }

  /**
   * Is this the first read in the pair (or only read if single-stranded)?
   */
  inline bool is_read1() const {
    // read is not read2: two flags for read1/read2, so read2 only
    // if read1 is 0 and read2 is 1, which can be expressed by
    return (b_->core.flag & (BAM_FREAD1 | BAM_FREAD2)) != BAM_FREAD2;
  }

  /**
   * Is the aligned read marked as being reversed?
   */
  inline bool is_reversed() const {
    return static_cast<bool>(b_->core.flag & BAM_FREVERSE);
  }

  /**
   * information about junctions from cigar array
   */
  template <uint32_t min_overhang>
  CigarJunctions<min_overhang> cigar_junctions() const {
    return CigarJunctions<min_overhang>(
        pos(), read_length(), bam_get_cigar(b_), b_->core.n_cigar);
  }
  /**
   * information about introns from cigar array
   */
  template <uint32_t min_overhang>
  CigarIntrons<min_overhang> cigar_introns() const {
    return CigarIntrons<min_overhang>(
        pos(), read_length(), bam_get_cigar(b_), b_->core.n_cigar);
  }

  /**
   * given experimental protocol, was it +/- strand?
   */
  GeneStrandness strandness(ExperimentStrandness experiment) {
    // unstranded experiment has ambiguous strandness
    if (experiment == ExperimentStrandness::NONE) {
        return GeneStrandness::AMBIGUOUS;
    }
    // forward strand if read1 is in the same direction as the experiment
    const bool read1_forward = is_read1() != is_reversed();
    const bool experiment_forward = experiment == ExperimentStrandness::FORWARD;
    return (read1_forward == experiment_forward) ?
      GeneStrandness::FORWARD : GeneStrandness::REVERSE;
  }

  // initialize and destroy memory used by htslib alignment object
  AlignmentRecord() : b_(bam_init1()) { }
  ~AlignmentRecord() { bam_destroy1(b_); }
  friend class AlignmentFile;  // allow AlignmentFile direct access to pointer
};

class AlignmentFile {
 private:
  samFile *in_;
  sam_hdr_t *header_;
  htsThreadPool p_;

  /**
   * Free memory for in_, header_, p_ given error in constructor or destructor
   */
  void cleanup_hts();

 public:
  /**
   * array of targets
   */
  char** target_name() const { return header_->target_name; }
  int32_t n_targets() const { return header_->n_targets; }

  /**
   * Update AlignmentRecord by mutable reference with next alignment in file
   *
   * @return >=0 on success, -1 if end of stream, < -1 if error
   */
  int NextRead(AlignmentRecord& b) { return sam_read1(in_, header_, b.b_); }

  AlignmentFile(const char* input_path, int nthreads)
      : in_(nullptr),
        header_(nullptr),
        p_{nullptr, 0} {
    // load samfile
    if (nullptr == (in_ = sam_open(input_path, "r"))) {
      cleanup_hts();
      std::ostringstream oss;
      oss << "Failed to open alignment file " << input_path;
      throw std::runtime_error(oss.str());
    }
    // load header
    if (nullptr == (header_ = sam_hdr_read(in_))) {
      cleanup_hts();
      std::ostringstream oss;
      oss << "Failed to read header from " << input_path;
      throw std::runtime_error(oss.str());
    }
    // do we set up threading?
    if (nthreads > 1) {
      // we use 1 less thread for I/O
      if (nullptr == (p_.pool = hts_tpool_init(nthreads))) {
        cleanup_hts();
        std::ostringstream oss;
        oss << "Failed to create thread pool for reading " << input_path;
        throw std::runtime_error(oss.str());
      }
      // otherwise, p.pool has been set, so we need to use it
      // attach threadpool to bam file we are reading
      int set_opt_code = hts_set_opt(in_, HTS_OPT_THREAD_POOL, &p_);
      if (set_opt_code != 0) {
        cleanup_hts();
        std::ostringstream oss;
        oss << "Failed to attach thread pool when reading " << input_path;
        throw std::runtime_error(oss.str());
      }
    }
  }
};

// implementation cleaning up pointers in AlignmentFile
void AlignmentFile::cleanup_hts() {
  if (nullptr != in_) {
    int sam_close_ret = sam_close(in_);
    if (sam_close_ret < 0) {
      std::ostringstream oss;
      oss << "Error closing alignment file (code:" << sam_close_ret << ")";
      throw std::runtime_error(oss.str());
    }
  }
  if (nullptr != header_) {
    sam_hdr_destroy(header_);
  }
  if (nullptr != p_.pool) {
    hts_tpool_destroy(p_.pool);
  }
}

}  // namespace bam
}  // namespace majiq

#endif  // MAJIQ_BAM_ALIGNMENTS_HPP
