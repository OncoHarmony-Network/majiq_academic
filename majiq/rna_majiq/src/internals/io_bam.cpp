#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <random>
#include <string>
#include <cmath>
#include <algorithm>
#include "io_bam.hpp"
#include "majiq_utils.hpp"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "boost/math/distributions/poisson.hpp"


// initialize static random state
vector<std::mt19937> io_bam::IOBam::generators_ = vector<std::mt19937>(0);

using namespace std;
namespace io_bam {
    unsigned int eff_len_from_read_length(int read_length) {
        if (read_length <= 2 * MIN_BP_OVERLAP) {
            return 0;
        } else {
            return read_length + 1 - (2 * MIN_BP_OVERLAP);
        }
    }

    char IOBam::_get_strand(bam1_t * read){
        char strn = '.';

        if (strandness_ == FWD_STRANDED){
            strn = ((!is_read_reverse(read) && is_read1(read)) ||
                    (is_read_reverse(read) && is_read2(read)) ||
                    (!is_read_reverse(read) && (is_read1(read) == is_read2(read)))) ? '+' : '-' ;

        } else if (strandness_ == REV_STRANDED){

            strn = ((is_read_reverse(read) && is_read1(read)) ||
                    (!is_read_reverse(read) && is_read2(read)) ||
                    (is_read_reverse(read) && (is_read1(read) == is_read2(read)))) ? '+' : '-' ;
        }
        return (strn);
    }

    inline int _unique(bam1_t * read){
        return ((read->core.flag & 0x100) != 0x100) ;
    }

    inline int _unmapped(bam1_t * read){
        return ((read->core.flag & 0x4) == 0x4) ;
    }

    // returns true if associated cigar operation offsets position on read
    inline bool _cigar_advance_read(const char cigar_type) {
        return cigar_type & 1;
    }

    // returns true if associated cigar operation offsets position on reference
    inline bool _cigar_advance_reference(const char cigar_type) {
        return cigar_type & 2;
    }

    /**
     * adjust cigar_ptr, n_cigar to ignore clipping, return length of clipping on left and right
     *
     * @param n_cigar number of cigar operations passed by reference (adjusted by call)
     * @param cigar array of cigar operations passed by reference (adjusted by call)
     *
     */
    pair<int, int> _adjust_cigar_soft_clipping(int &n_cigar, uint32_t* &cigar) {
        // initialize lengths of soft clipping on left/right to return
        int left_length = 0;
        int right_length = 0;
        // remove clipping on the right
        for (int i = n_cigar - 1; i >= 0; --i) {  // iterate backwards
            const char cigar_op = bam_cigar_op(cigar[i]);
            if (cigar_op == BAM_CHARD_CLIP) {
                // ignore hard clipping cigar operations on right
                --n_cigar;
            } else if (cigar_op == BAM_CSOFT_CLIP) {
                // ignore soft clipping cigar operations, update right length
                --n_cigar;
                right_length += bam_cigar_oplen(cigar[i]);
            } else {
                break;
            }
        }
        // remove clipping on the left
        int lhs_clipping = 0;  // offset to apply to *cigar_ptr
        for (int i = 0; i < n_cigar; ++i) {
            const char cigar_op = bam_cigar_op(cigar[i]);
            if (cigar_op == BAM_CHARD_CLIP) {
                // ignore hard clipping cigar operations on left
                ++lhs_clipping;
            } else if (cigar_op == BAM_CSOFT_CLIP) {
                // ignore soft clipping cigar operations, update left length
                ++lhs_clipping;
                left_length += bam_cigar_oplen(cigar[i]);
            } else {
                break;
            }
        }
        if (lhs_clipping > 0) {
            n_cigar -= lhs_clipping;
            cigar = cigar + lhs_clipping;
        }
        return make_pair(left_length, right_length);
    }

    /**
     * Update effective length, appropriately resizing things
     *
     * @note this should not be called after introns are inserted because
     * positional information mapped to bins will lose meaning
     *
     * @note XXX this is not thread-safe! assumes that there are no parallel
     * threads accessing eff_len_or junc_vec. This is true when parsing BAM
     * files. Otherwise, we would need to add locks on eff_len_ and junc_vec.
     * Alternatively, we could use tbb::concurrent_vector to avoid the need for
     * locks if we were okay with adding a new library
     */
    void IOBam::update_eff_len(const unsigned int new_eff_len) {
        if (new_eff_len > eff_len_) {
            eff_len_ = new_eff_len;
            if (eff_len_ > buff_len_) {
                // we need to expand buffers
                // give safety factor but make sure no less than eff_len_
                buff_len_ = std::max(
                        eff_len_,
                        static_cast<unsigned int>(eff_len_ * BUFF_LEN_SAFETY_FACTOR)
                );
                // resize all the buffers in parallel
                #pragma omp parallel for num_threads(nthreads_)
                for (unsigned int i = 0; i < junc_vec.size(); ++i) {
                    junc_vec[i]->resize(buff_len_, 0);
                }
            }
        }
    }

    /**
     * Associate junction nreads_ptr with appropriate genes
     */


    void IOBam::find_junction_genes(
            string chrom, char strand, int start, int end,
            shared_ptr<vector<float>> nreads_ptr) {
        // find the set of overlapping genes that this junction could overlap with
        vector<overGene*>::iterator low = lower_bound(
                glist_[chrom].begin(),
                glist_[chrom].end(),
                start,
                _Region::func_comp);

        // what is the set of all possible genes?
        if (low == glist_[chrom].end()) {
            // junction is after all genes for the contig
            return;
        }
        // otherwise, the set of all possible genes is:
        const auto& possible_genes = (*low)->glist;

        // we can have 3 kinds of matches of our junction to the genes in this
        // list:
        // Stage 1: the gene already has the junction
        // Stage 2: the junction start and end are "close" to exons in the gene
        //  (within MAX_DENOVO_DIFFERENCE) of exon boundary in correct
        //  direction but not stage 1
        // Stage 3: the junction is intragenic (contained in gene boundaries)
        //  but not stage 1 or 2
        bool found_stage1 = false;  // did we already find stage 1 genes?
        vector<Gene*> stage2_genes;  // genes that match stage 2
        vector<Gene*> stage3_genes;  // genes that match stage 3
        // get key for junction to match against
        Junction junc(start, end, false, simpl_);
        const coord_key_t key = junc.get_key();

        // the junction could be in any of the genes from the first set
        for (const auto &gObj: possible_genes) {
            bool utr_gene_junc = false;
            // try to rule out being any stage...
            if (strand != '.' && strand != gObj->get_strand()) {
                // junction must be matched to gene that matches its strand
                continue;
            }

            if (start < gObj->get_utr_start() || end > gObj->get_utr_end()) {
                // if either junction end falls outside of the possible extended area around the gene, ignore it
                //cout << "SKIPPED JUNC UTR " << gObj->get_name() << ' ' << start << '-' << end << ' ' << gObj->get_utr_start()  << '-' << gObj->get_utr_end() << '\n';
                continue;
            }
            if (start < gObj->get_start() && end > gObj->get_end()) {
                // if both ends of the junction fall outside of the gene or UTR end, also ignore it
                //cout << "SKIPPED JUNC INTRA " << gObj->get_name() << ' ' << start << '-' << end << ' ' << gObj->get_start() << '-' << gObj->get_end() << '\n';
                continue;
            }
            if (start < gObj->get_start() || end > gObj->get_end()) {
                // if either junction end fall outside of the gene or UTR end, set as potential UTR junction
                utr_gene_junc = true;
            }

            // junction is matched to this gene in either stage 1, 2, or 3
            // Stage 1: is the junction already present for the gene?
            if (gObj->junc_map_.count(key) > 0) {
                found_stage1 = true;
                // stage 1 genes initialized as they are found
                gObj->initialize_junction(key, start, end, nreads_ptr, simpl_);
            } else if (found_stage1) {
                // there is no point in distinguishing between stage 2 and 3,
                // because we will only use stage 1 matches
                continue;
            } else {

                // this section searches the UTR junc
                if(!allow_full_intergene_ && utr_gene_junc){
                    bool start_in_exon = false;
                    bool end_in_exon = false;
                    for (const auto &ex: gObj->exon_map_) {
                        // compare junction to exon boundaries
                        const int exon_start = ex.second->get_start();
                        const int exon_end = ex.second->get_end();
                        // skip half exons
                        if (exon_start < 0 || exon_end < 0) {
                            continue;
                        }
                        // is the junction close the the full-exon boundaries?
                        if (!start_in_exon) {
                            if (start <= exon_start) {
                                // this nor subsequent exons will be close to
                                // start (not stage 2)
                                break;
                            } else if (start <= exon_end) {
                                start_in_exon = true;
                                // still need to see if end is near exon
                            }
                        }
                        if (end < exon_start) {
                            // this nor subsequent exons will be close to end (not
                            // stage 2)
                            break;
                        } else if (end < exon_end) {
                            end_in_exon = true;
                            break;
                        }
                    }

                    if(!start_in_exon && !end_in_exon){
                        // we ignore UTR junctions with neither end near an annotated exon
                        continue;
                    }
                }


                // unless we find a stage 1 gene, we need to know if this is
                // stage 2 or stage 3 gene. So check whether start and end are
                // close enough to any of the gene's exons
                bool start_near_exon = false;  // is junction start near exon?
                bool end_near_exon = false;  // is junction end near exon?
                for (const auto &ex: gObj->exon_map_) {
                    // compare junction to exon boundaries
                    const int exon_start = ex.second->get_start();
                    const int exon_end = ex.second->get_end();
                    // skip half exons
                    if (exon_start < 0 || exon_end < 0) {
                        continue;
                    }
                    // is the junction close the the full-exon boundaries?
                    if (!start_near_exon) {
                        if (start <= exon_start) {
                            // this nor subsequent exons will be close to
                            // start (not stage 2)
                            break;
                        } else if (start <= exon_end + MAX_DENOVO_DIFFERENCE) {
                            start_near_exon = true;
                            // still need to see if end is near exon
                        }
                    }
                    if (end < exon_start - MAX_DENOVO_DIFFERENCE) {
                        // this nor subsequent exons will be close to end (not
                        // stage 2)
                        break;
                    } else if (end < exon_end) {
                        end_near_exon = true;
                        break;
                    }
                }




                // was this stage 2 or stage 3?
                if (start_near_exon && end_near_exon) {
                    stage2_genes.push_back(gObj);
                } else {
                    stage3_genes.push_back(gObj);
                }
            }
        }
        // initialize genes for this junction if there were no stage 1 genes
        // (stage 1 genes are initialized online when found)
        if (!found_stage1) {
            // initialize stage2 genes if any found, otherwise any stage3 genes
            for (const auto &g: (stage2_genes.size() > 0 ? stage2_genes : stage3_genes)) {
                g->initialize_junction(key, start, end, nreads_ptr, simpl_);
            }
        }
        return;
    }


    /**
     * Add junction at position/strand to junction coverage array
     *
     * @param chrom sequence id where junction was found
     * @param strand information of if junction count was stranded/its direction
     * @param start, end boundaries of the junction (1-based coordinates)
     * @param offset position to add reads to
     * genomic-position-based offset.
     * @param sreads number of reads to add to junction at stranded position
     *
     */
    void IOBam::add_junction(string chrom, char strand, int start, int end, const unsigned int offset, int sreads) {
        if (offset >= eff_len_ || offset < 0) {
            // XXX -- should warn that we are losing reads (should handle offsets
            // earlier, this function should not need to check a second time
            return;
        }
        // get string-based key to map over junctions (index of 2d coverage vector)
        string key = chrom + ":" + strand + ":" + to_string(start) + "-" + to_string(end) ;

        // indicator if this is the first time we have seen this junction
        bool new_j = false ;
        // acquire lock for working with junc_map/junc_vec
        omp_set_lock(&map_lck_) ;
        {
            if (junc_map.count(key) == 0 ) {
                // we haven't seen this junction before
                junc_map[key] = junc_vec.size();  // index to add to junc_vec
                junc_vec.push_back(make_shared<vector<float>>(buff_len_, 0));
                // mark as new junction
                new_j = true ;
            }
        }
        omp_unset_lock(&map_lck_) ;
        shared_ptr<vector<float>> v_ptr = junc_vec[junc_map[key]];
        vector<float> &v = *v_ptr;

        if (new_j) {
            // associate junction with genes (prioritizing annotated, other stages for denovo)
            find_junction_genes(chrom, strand, start, end, v_ptr);
        }
        #pragma omp atomic
            v[offset] += sreads;  // add reads to junction array
        return ;
    }


    /**
     * Parse input read for split reads and add to position coverage array
     *
     * @param header SAM header with information sequence ids (chromosomes)
     * @param read SAM alignment being parsed
     *
     * @return 0 after adding junctions from primary alignments that have
     * sufficient overlap into the read
     *
     */
    int IOBam::parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) {
        if (!_unique(read) || _unmapped(read)) {
            // ignore unmapped reads and secondary alignments
            return 0;
        }

        // get genomic position of start of read (zero-based coordinate)
        string chrom(header->target_name[read->core.tid]);
        int read_pos = read->core.pos;

        // get length of read as initial length of alignment
        const int read_length = read->core.l_qseq;

        // get cigar operations, adjust them and alignment length to ignore clipping
        int n_cigar = read->core.n_cigar;
        uint32_t *cigar = bam_get_cigar(read) ;
        const pair<int, int> clipping_lengths = _adjust_cigar_soft_clipping(
                n_cigar, cigar
        );
        const int alignment_length = read_length - clipping_lengths.first - clipping_lengths.second;

        // track offsets on genomic position and alignment as we iterate over cigar
        int genomic_offset = 0;
        int alignment_offset = 0;

        // iterate over cigar operations, look for split reads
        for (int i = 0; i < n_cigar; ++i) {
            // extract operation and offset for the operation from cigar[i]
            const char cigar_op = bam_cigar_op(cigar[i]);
            const int cigar_oplen = bam_cigar_oplen(cigar[i]);
            // check if we have a junction we can add
            if (cigar_op == BAM_CREF_SKIP && alignment_offset >= MIN_BP_OVERLAP) {
                // this is a junction overlapping sufficiently with the alignment
                // (we don't need to check alignment_offset <= alignment_length - MIN_BP_OVERLAP
                // because we exit early when advancing alignment offset)
                // get genomic start and end of junction
                const int junction_start = read_pos + genomic_offset;
                const int junction_end = junction_start + cigar_oplen + 1;
                // get junction position relative to read offset (not alignment
                // offset). This requires adjustment for soft clipping on left
                const unsigned int junction_pos = alignment_offset
                    + clipping_lengths.first - MIN_BP_OVERLAP;
                // add the junction at the specified position
                try {
                    add_junction(chrom, _get_strand(read), junction_start,
                                 junction_end, junction_pos, 1);
                } catch (const std::logic_error& e) {
                    cout << "ERROR" << e.what() << '\n';
                }
            }
            // update alignment and genomic offsets
            const char cigar_type = bam_cigar_type(cigar[i]);
            if (_cigar_advance_read(cigar_type)) {
                alignment_offset += cigar_oplen;
                if (alignment_offset > alignment_length - MIN_BP_OVERLAP) {
                    break;  // future operations will not lead to valid junction
                }
            }
            if (_cigar_advance_reference(cigar_type)) {
                genomic_offset += cigar_oplen;
            }
        }
        // done processing this read
        return 0;
    }


    /**
     * Parse input read for introns and add to position coverage array
     *
     * @param header SAM header with information sequence ids (chromosomes)
     * @param read SAM alignment being parsed
     *
     * @return 0 after adding introns from primary alignments that have
     * sufficient overlap into the read
     *
     */
    int IOBam::parse_read_for_ir(bam_hdr_t *header, bam1_t *read) {
        if (!_unique(read) || _unmapped(read)) {
            // ignore unmapped reads and secondary alignments
            return 0;
        }

        // get genomic position of start of read (zero-based coordinate)
        string chrom(header->target_name[read->core.tid]) ;
        if (intronVec_.count(chrom) == 0) {
            // if there is nothing to count in this chromosome, don't bother
             return 0;
        }
        int read_pos = read->core.pos;

        // get iterator with all introns with end position after alignment start
        vector<Intron *>::iterator low = intronVec_[chrom].begin()
            // get lower-bound using cumulative maximum to meet conditions to
            // do correct O(log(n)) binary search using std::lower_bound
            + (lower_bound(intronVec_end_cummax_[chrom].begin(),
                           intronVec_end_cummax_[chrom].end(), read_pos)
               - intronVec_end_cummax_[chrom].begin());
        if (low == intronVec_[chrom].end()) {
            // this read can't contribute to any of the introns, so exit here
            return 0;
        }

        // fill records of genomic/read offsets to determine offset of introns
        vector<pair<int, int>> genomic_alignment_offsets;
        genomic_alignment_offsets.push_back({read_pos, 0});
        // fill records of junctions splitting read (ignore introns if intersect)
        vector<pair<int, int>> junc_record ;
        // first and last genomic positions an intron can include to be considered
        int first_pos = -10;  // -10 is placeholder value
        int last_pos = -10;

        // get length of read
        const int read_length = read->core.l_qseq;

        // get cigar operations, adjust them and alignment length to ignore clipping
        int n_cigar = read->core.n_cigar ;
        uint32_t *cigar = bam_get_cigar(read) ;
        const pair<int, int> clipping_lengths = _adjust_cigar_soft_clipping(
                n_cigar, cigar
        );
        const int alignment_length = read_length - clipping_lengths.first - clipping_lengths.second;

        // parse CIGAR operations for genomic_alignment_offsets, junc_record, {first,last}_pos
        int genomic_offset = 0;
        int alignment_offset = 0;
        for (int i = 0; i < n_cigar; ++i) {
            // extract operation and offset for the operation from cigar[i]
            const char cigar_op = bam_cigar_op(cigar[i]);
            const int cigar_oplen = bam_cigar_oplen(cigar[i]);
            if (cigar_op == BAM_CREF_SKIP) {
                // get start and end of junction
                const int junction_start = read_pos + genomic_offset;
                const int junction_end = junction_start + cigar_oplen + 1;
                try {
                    junc_record.push_back({junction_start, junction_end});
                } catch (const std::logic_error& e) {
                    cout << "ERROR" << e.what() << '\n';
                }
            }
            const char cigar_type = bam_cigar_type(cigar_op);
            const bool advance_read = _cigar_advance_read(cigar_type);
            const bool advance_reference = _cigar_advance_reference(cigar_type);
            if (advance_read) {
                // check if this operation gives us first/last position for valid intron
                if (
                        (alignment_offset < MIN_BP_OVERLAP)
                        && ((alignment_offset + cigar_oplen) >= MIN_BP_OVERLAP)
                ) {
                    // first valid value for intron_end
                    first_pos = read_pos + genomic_offset
                        + (advance_reference ? MIN_BP_OVERLAP - alignment_offset : 0);
                }
                if (
                        ((alignment_offset + cigar_oplen) > alignment_length - MIN_BP_OVERLAP)
                        && (alignment_offset <= alignment_length - MIN_BP_OVERLAP)
                ) {
                    // last valid value for intron_start
                    // don't forget +1 for 1-indexed (read_pos) to 1-indexed intron_start
                    last_pos = (1 + read_pos) + genomic_offset
                        + (advance_reference ? alignment_length - MIN_BP_OVERLAP - alignment_offset : 0);
                }
                alignment_offset += cigar_oplen;
            }
            if (advance_reference) {
                genomic_offset += cigar_oplen;
            }
            if (advance_reference != advance_read) {
                // only need update if one advanced but not the other
                const int current_genomic_pos = read_pos + genomic_offset;
                const int current_offset = alignment_offset;
                genomic_alignment_offsets.push_back({current_genomic_pos, current_offset});
            }
        }

        const char read_strand = _get_strand(read);

        // loop over introns
        for (; low != intronVec_[chrom].end(); ++low) {
            Intron * intron = *low;  // current intron from the iterator

            // is intron outside of admissible coordinates?
            if (intron->get_end() < first_pos) {
                continue;  // intron cannot overlap enough
            }
            if (intron->get_start() > last_pos) {
                break;  // this and all subsequent introns cannot overlap enough
            }

            // is read compatible with gene strand?
            const char gstrand = intron->get_gene()->get_strand() ;
            if (read_strand != '.' && gstrand != read_strand) {
                continue;  // read strand not compatible with intron gene strand
            }

            // does read include junctions overlapping intron?
            bool junc_found = false;
            for (const auto &j: junc_record) {
                if (j.first <= intron->get_end() && intron->get_start() <= j.second) {
                    // there is a junction that starts or ends within intronic
                    // coordinates, or there is a junction that goes over this
                    // intron entirely
                    junc_found = true;
                    break;
                } else if (j.second > intron->get_end()) {
                    // none of the remaining junctions (sorted) can intersect intron
                    break;
                }
            }
            if (junc_found) {
                continue;  // read not compatible with intron due to junctions
            }

            // since read includes intron within acceptable bounds, add intron
            // with appropriate offset
            // calculate relative to last acceptable position in alignment, which has offset alignment_length - MIN_BP_OVERLAP
            int relative_offset = 0;
            // don't forget 1 in front of read_pos, genomic offsets, which are
            // zero-indexed vs 1-indexed intron starts
            int read_start_vs_intron_start = (1 + read_pos) - intron->get_start();
            if (read_start_vs_intron_start >= 0) {
                // read start at or after intron start, further away from end of read
                relative_offset = read_start_vs_intron_start;
            } else {
                unsigned int i = 1;  // get index of first offset past intron start
                for (; i < genomic_alignment_offsets.size(); ++i) {
                    // could replace with something like std::lower_bound for log(n) search
                    // but premature optimization unless many cigar operations
                    if (genomic_alignment_offsets[i].first >= intron->get_start()) {
                        break;
                    }
                }
                relative_offset =
                    // how far last genomic coordinate before intron is from start
                    (1 + genomic_alignment_offsets[i - 1].first)  // last genomic coordinate as 1-indexed position
                    - intron->get_start()  // intron start (1-indexed position)
                    // adjust for offset on alignment
                    - genomic_alignment_offsets[i - 1].second;
            }
            // get intron offset relative to read (not alignment), which
            // requires adjusting for soft clipping on *right* (this is in
            // contrast to junctions -- intron coordinates are sort of inverted
            // relative to junctions for historical reasons)
            int intron_offset = alignment_length - MIN_BP_OVERLAP
                + relative_offset + clipping_lengths.second;
            #pragma omp critical
            {
                intron->add_read(intron_offset, eff_len_, 1);
            }
        }
        return 0 ;
    }


    /**
     * Set up file for reading with given number of threads
     */
    inline int _init_samfile(
            const char* bam_path,
            samFile* &in, bam_hdr_t* &header, bam1_t* &aln,
            htsThreadPool &p, const unsigned int nthreads
    ) {
        int ignore_sam_err = 0;
        int extra_hdr_nuls = 0;

        in = sam_open(bam_path, "rb") ;
        if (NULL == in) {
            fprintf(stderr, "Error opening \"%s\"\n", bam_path);
            return EXIT_FAILURE;
        }
        header = sam_hdr_read(in);
        if (NULL == header) {
            fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_path);
            sam_close(in);
            return EXIT_FAILURE;
        }
        header->ignore_sam_err = ignore_sam_err;
        if (extra_hdr_nuls) {
            char *new_text = (char*) realloc(header->text, header->l_text + extra_hdr_nuls);
            if (new_text == NULL) {
                fprintf(stderr, "Error reallocing header text\n");
                sam_close(in);
                bam_hdr_destroy(header);
                return EXIT_FAILURE;
            }
            header->text = new_text;
            memset(&header->text[header->l_text], 0, extra_hdr_nuls);
            header->l_text += extra_hdr_nuls;
        }

        aln = bam_init1();
        if (nthreads > 0) {
            p.pool = hts_tpool_init(nthreads);
            if (!p.pool) {
                fprintf(stderr, "Error creating thread pool\n");
            } else {
                hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
            }
        }
        return 0;
    }


    /**
     * Close file for reading
     */
    inline int _close_samfile(
            samFile* &in, bam_hdr_t* &header, bam1_t* &aln, htsThreadPool &p
    ) {
        int exit_code = 0;
        // close the read/header/input file, check if there was an error closing input
        bam_destroy1(aln);
        bam_hdr_destroy(header);
        if (sam_close(in) < 0) {
            fprintf(stderr, "Error closing input.\n");
            exit_code = 1 ;
        }
        if (p.pool)
            hts_tpool_destroy(p.pool);

        return exit_code;
    }


    /**
     * Obtain lower bound on eff_len from first num_reads, update eff_len
     */
    void IOBam::EstimateEffLenFromFile(int num_reads) {
        samFile *in;
        bam_hdr_t *header ;
        bam1_t *aln;
        htsThreadPool p = {NULL, 0};
        // open the sam file
        int exit_code = _init_samfile(bam_.c_str(), in, header, aln, p, nthreads_);
        if (exit_code) {
            cerr << "Error parsing input\n";
            return;
        }
        int read_ct = 0;
        int max_read_length = 0;
        int r = 0;  // error code from htslib
        // loop over the first num_reads reads
        while (read_ct < num_reads && (r = sam_read1(in, header, aln)) >= 0) {
            ++read_ct;  // increment the number of reads proceessed
            const int read_length = aln->core.l_qseq;  // current read length
            if (read_length > max_read_length) {
                max_read_length = read_length;
            }
        }
        if (r < -1) {
            cerr << "Error parsing input\n";
            return;
        }
        // close the sam file
        exit_code |= _close_samfile(in, header, aln, p);
        if (exit_code) {
            cerr << "Error closing input\n";
        }
        // update eff_len_ if greater than current value
        const unsigned int est_eff_len = eff_len_from_read_length(max_read_length);
        if (est_eff_len > eff_len_) {
            update_eff_len(est_eff_len);
        }
    }


    int IOBam::ParseJunctionsFromFile(bool ir_func) {
        samFile *in;
        bam_hdr_t *header ;
        bam1_t *aln;
        htsThreadPool p = {NULL, 0};
        // open the sam file
        int exit_code = _init_samfile(bam_.c_str(), in, header, aln, p, nthreads_);
        if (exit_code) {
            return exit_code;  // there was an error, so end with error code
        }

        int (IOBam::*parse_func)(bam_hdr_t *, bam1_t *) ;
        if(ir_func) parse_func = &IOBam::parse_read_for_ir ;
        else parse_func = &IOBam::parse_read_into_junctions ;

        int r = 0;  // error code from htslib
        while ((r = sam_read1(in, header, aln)) >= 0) {
            (this->*parse_func)(header, aln) ;
        }
        if (r < -1) {
            fprintf(stderr, "Error parsing input.\n");
            exit_code = 1;
        }
        if(!ir_func)
            junc_limit_index_ = junc_vec.size() ;

        // close the sam file
        exit_code |= _close_samfile(in, header, aln, p);

        return exit_code ;
    }

    /**
     * sort vec, identify positions 0:numpos that are not stacks, return numpos
     *
     * @param &vec vector of nonzero coverage values, sorted as side effect
     * @param pvalue_limit positions with coverage with right-tailed p-value
     * less than this limit will be marked as stacks
     *
     * @return the number of positions that are not stacks. The coverage at
     * these positions are found in sorted order as the corresponding first
     * values in vec
     *
     * @note Sorts vec as side effect to test only the most extreme values of
     * vec
     * @note Does not remove the stacks from vec. The return value gives enough
     * information to ignore these values in subsequent computation
     * @note Stacks are computed with respect to leave-one-out mean of all
     * other positions before any stacks removed. This means that previously
     * removed stacks will contribute to the distribution mean. It would be
     * trivial to modify this so that stack removal takes place against the
     * running leave-one-out mean after each stack removal
     */
    unsigned int IOBam::normalize_stacks(
            vector<float> &vec, const float pvalue_limit) {
        using boost::math::poisson_distribution;
        using boost::math::cdf;
        using boost::math::complement;

        // sort coverage so we only need to test extremes
        sort(vec.begin(), vec.end());
        // get sum of reads
        float sreads = std::accumulate(vec.begin(), vec.end(), float{});
        // number of positions that remain/for denominator
        unsigned int numpos = vec.size();
        const unsigned int other_npos = numpos - 1;  // denominator for leave-one-out mean
        // starting from highest coverage positions, identify stacks
        for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
            // coverage at current position
            const float vec_i = *it;
            // leave-one-out mean coverage at other positions
            const float mean_reads = (other_npos == 0) ? 0.5: (sreads - vec_i) / other_npos;
            // p-value at current position
            const poisson_distribution<float> stack_dist(mean_reads);
            const float pvalue = cdf(complement(stack_dist, vec_i));
            // decrement number of valid positions if outlier
            if (pvalue < pvalue_limit) {
                --numpos;  // decrement numpos for return value
            } else {
                break;  // don't need to test remaining positions
            }
        }
        // finally: values vec[:numpos] are not stacks, vec[numpos:] are stacks
        return numpos;
    }

    int IOBam::bootstrap_samples(int msamples, float* boots, float pvalue_limit) {
        const int njunc = junc_map.size();

        #pragma omp parallel for num_threads(nthreads_)
        for(int jidx=0; jidx < njunc; jidx++){

            vector<float> &coverage = *(junc_vec[jidx]);
            // copy nonzero coverage values to vec
            std::vector<float> vec;
            std::copy_if(coverage.begin(), coverage.end(),
                std::back_inserter(vec), [](float x) { return x > 0; });
            // if there are no nonzero positions, skip...
            if (vec.empty()) continue;

            // get number of positions to bootstrap over after stack removal
            const unsigned int numpos = pvalue_limit <= 0 ?
                    vec.size() : normalize_stacks(vec, pvalue_limit);
            // handle cases for getting bootstrap replicates
            if (numpos == 0) {
                // can't bootstrap from 0 positions (everything is default 0)
                continue;
            }
            // we perform parametric vs nonparametric bootstrap sampling
            // depending on whether variance of bootstrap distribution is
            // less than mean
            // obviously parametric if constant
            bool parametric_bootstrap = (numpos == 1 || vec[0] == vec[numpos - 1]);
            float bootstrap_mean;
            if (parametric_bootstrap) {
                bootstrap_mean = vec[0];
            } else {
                bootstrap_mean = majiq::mean(vec.begin(), vec.begin() + numpos);
                // variance<0> for ddof=0 --> population variance
                const float bootstrap_variance = majiq::variance<0>(
                    vec.begin(), vec.begin() + numpos, bootstrap_mean);
                parametric_bootstrap = bootstrap_variance < bootstrap_mean;
            }

            std::mt19937 &generator = generators_[omp_get_thread_num()];
            if (parametric_bootstrap) {
                // variance < mean if we bootstrapped nonparametrically, so use
                // Poisson distribution instead as conservative estimate of
                // overdispersed variance
                poisson_distribution<int> distribution(bootstrap_mean * numpos);
                for (int m = 0; m < msamples; ++m) {
                    const int idx2d = (jidx * msamples) + m;
                    boots[idx2d] = static_cast<float>(distribution(generator));
                }
            } else {  // performing nonparametric sampling where overdispersed
                // by sampling numpos - 1 positions to sum, bootstrap samples will
                // have correct variance (variance of difference between
                // two draws from bootstrap distribution equals difference
                // between two draws from generating distribution under
                // appropriate assumptions)
                const unsigned int ksamples = numpos - 1;
                // rescale sum to be comparable to sum over all nonzero positions
                const float scale_sum = static_cast<float>(numpos) / ksamples;
                // generate random numbers
                uniform_int_distribution<unsigned int> distribution(0, numpos - 1);
                for (int m = 0; m < msamples; ++m) {
                    // index to update
                    const int idx2d = (jidx * msamples) + m;
                    // accumulate ksamples reads from nonzero positions
                    float lambda = 0;
                    for (unsigned int k = 0; k < ksamples; ++k) {
                        lambda += vec[distribution(generator)];
                    }
                    // save rescaled sum over bootstrapped positions to output buffer
                    boots[idx2d] = lambda * scale_sum;
                }
            }
        }
        return 0 ;
    }

    int IOBam::get_njuncs(){
        return junc_vec.size() ;
    }

    const map<string, unsigned int>& IOBam::get_junc_map(){
        return junc_map ;
    }

    int * IOBam::get_junc_vec_summary(){
        int njunc = junc_vec.size() ;
        int * res = (int *) calloc(2*njunc, sizeof(int)) ;

        #pragma omp parallel for num_threads(nthreads_)
        for(int jidx=0; jidx < njunc; ++jidx){
            vector<float> &coverage = *(junc_vec[jidx]);
            float ss = 0 ;
            int np = 0 ;
            for(unsigned int i=0; i<eff_len_; ++i){
                ss += coverage[i];
                np += (coverage[i] > 0) ? 1 : 0;
            }
            res[jidx] = (int) ss ;
            res[jidx+njunc] = np ;
        }
        return res ;
    }


    void IOBam::parseJuncEntry(map<string, vector<overGene*>> & glist, string gid, string chrom, char strand,
                               int start, int end, unsigned int sreads, unsigned int minreads_t, unsigned int numpos,
                               unsigned int minpos_t, unsigned int denovo_t, bool denovo, bool denovo_ir, vector<Gene*>& oGeneList,
                               bool ir, vector<float>& ircov, float min_intron_cov, float min_bins, int minexp,
                               bool reset){

        vector<overGene*>::iterator low = lower_bound (glist[chrom].begin(), glist[chrom].end(),
                                                       start, _Region::func_comp ) ;
        if (low == glist[chrom].end())
            return ;
        if (ir){
            for (const auto &gObj: (*low)->glist){
                if (gObj->get_id() != gid) continue ;
                Intron * irptr = new Intron(start, end, false, gObj, simpl_) ;
                irptr->add_read_rates_buff(ircov.size()) ;
                irptr->initReadCovFromVector(ircov) ;

                const string key = irptr->get_key(gObj) ;
                omp_set_lock(&map_lck_) ;
                {
                    if (junc_map.count(key) == 0) {
                        junc_map[key] = junc_vec.size() ;
                        junc_vec.push_back(irptr->read_rates_ptr_) ;
                        gObj->add_intron(irptr, min_intron_cov, minexp, min_bins, reset, denovo_ir) ;
                    }
                }
                omp_unset_lock(&map_lck_) ;
            }
        } else {
            coord_key_t key = std::make_pair(start, end);
            add_junction(chrom, strand, start, end, 0, sreads);
            for (const auto &gObj: (*low)->glist){
                gObj->updateFlagsFromJunc(key, sreads, minreads_t, numpos, minpos_t, denovo_t, denovo, minexp, reset) ;
            }

        }
        return ;

    }

    void IOBam::detect_introns(float min_intron_cov, unsigned int min_experiments, float min_bins, bool reset, bool denovo){
        for (const auto & it: glist_){
            if (intronVec_.count(it.first)==0){
                intronVec_[it.first] = vector<Intron*>() ;
            }
            const int n = (it.second).size() ;
            #pragma omp parallel for num_threads(nthreads_)
            for(int g_it = 0; g_it<n; g_it++){
                for (const auto &g: ((it.second)[g_it])->glist){
                    g->detect_introns(intronVec_[it.first], simpl_) ;
                }
            }
            // get introns per chromosome in sorted order
            sort(intronVec_[it.first].begin(), intronVec_[it.first].end(), Intron::islowerRegion<Intron>) ;
            // get cumulative maximum for previous intron ends
            vector<int> cumulative_max;
            cumulative_max.reserve(intronVec_[it.first].size());
            int max_value = -1;  // initial maximum value
            for (const auto & intron : intronVec_[it.first]) {
                const int intron_end = intron->get_end();
                if (intron_end > max_value) {
                    max_value = intron_end;
                }
                cumulative_max.push_back(max_value);
            }
            intronVec_end_cummax_[it.first] = cumulative_max;
        }
        ParseJunctionsFromFile(true) ;
        for (const auto & it: intronVec_){
            const int n = (it.second).size() ;

            #pragma omp parallel for num_threads(nthreads_)
            for (int idx = 0; idx < n; idx++) {
                Intron * intrn_it = (it.second)[idx] ;
                const bool pass = intrn_it->is_reliable(min_bins, eff_len_);
                if (pass) {
                    intrn_it->normalize_readrates();  // normalize readrates if it passed
                    const string key = intrn_it->get_key(intrn_it->get_gene()) ;
                    #pragma omp critical
                    {
                        if (junc_map.count(key) == 0) {
                            junc_map[key] = junc_vec.size() ;
                            junc_vec.push_back(intrn_it->read_rates_ptr_) ;
                            (intrn_it->get_gene())->add_intron(intrn_it, min_intron_cov, min_experiments, min_bins,
                                                               reset, denovo) ;
                        }
                    }
                } else {
                    intrn_it->free_nreads() ;
                    delete intrn_it ;
                }
            }
        }
    }


    void IOBam::get_intron_raw_cov(float* out_cov){
        const unsigned int all_junc = junc_vec.size() ;
        #pragma omp parallel for num_threads(nthreads_)
        for(unsigned int idx = junc_limit_index_; idx<all_junc; idx++){
            const vector<float> &coverage = *(junc_vec[idx]);
            const unsigned int i = idx - junc_limit_index_ ;
            for (unsigned int j=0; j<eff_len_; j++){
                const unsigned int indx_2d = i*eff_len_ + j ;
                out_cov[indx_2d] = coverage[j] ;
            }
        }
    }

    void IOBam::get_junction_raw_cov(float* out_cov) {
        #pragma omp parallel for num_threads(nthreads_)
        for(unsigned int i = 0; i < junc_limit_index_; ++i) {
            const vector<float> &coverage = *(junc_vec[i]);
            for (unsigned int j = 0; j < eff_len_; ++j) {
                const unsigned int indx_2d = i * eff_len_ + j;
                out_cov[indx_2d] = coverage[j];
            }
        }
    }

    /**
     * Partition genes into chromosomes and groups of intersecting genes (overGenes) in sorted order
     */
    void prepare_genelist(map<string, Gene*>& gene_map, map<string, vector<overGene*>> & geneList){
        // construct lists of genes per chromosome
        map<string, vector<Gene*>> gmp;  // genes per chromosome
        for(const auto & p: gene_map){
            const string chrom = (p.second)->get_chromosome() ;
            if (gmp.count(chrom) == 0 )
                gmp[chrom] = vector<Gene*>() ;
            gmp[chrom].push_back(p.second)  ;
        }
        // for each chromosome/list of genes
        for (auto &gl: gmp){
            // iterate over genes in sorted order
            sort((gl.second).begin(), (gl.second).end(), Gene::islowerRegion<Gene>) ;
            // group genes (into overGene) that directly/indirectly overlap
            int gstart = 0, gend = 0 ;
            overGene * ov = nullptr ;

            for (const auto &g: (gl.second)){

                int st = g->get_utr_start() ;
                int nd = g->get_utr_end() ;
                // if current gene overlaps previous overgene, add it
                if (st <= gend && nd >= gstart){
                    // update overgene boundaries to include new gene
                    gstart = (gstart == 0) ? st : min(gstart, st) ;
                    gend = (gend == 0) ? nd : max(gend, nd) ;
                    ov->set_start(gstart) ;
                    ov->set_end(gend) ;
                }
                else{
                    // otherwise, done constructing current overgene, start new
                    // one for current gene
                    if (ov != nullptr )
                        geneList[gl.first].push_back(ov) ;
                    ov = new overGene(st, nd) ;
                    gstart = st ;
                    gend = nd ;
                }
                (ov->glist).push_back(g) ;
            }
            if(ov != nullptr)
                geneList[gl.first].push_back(ov) ;
        }
    }


    void free_genelist(map<string, vector<overGene*>> & geneList){
        for (auto &ov_vec: geneList){
            for(auto &ov: ov_vec.second){
                delete (ov) ;
            }
            (ov_vec.second).clear() ;
        }
    }

}
