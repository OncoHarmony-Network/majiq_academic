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
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "scythestat/distributions.h"

using namespace std;
namespace io_bam {

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
    inline bool _cigar_advance_read(const char cigar_op) {
        return (
            cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF
            || cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP
        );
    }

    // returns true if associated cigar operation offsets position on reference
    inline bool _cigar_advance_reference(const char cigar_op) {
        return (
            cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF
            || cigar_op == BAM_CDEL || cigar_op == BAM_CREF_SKIP
        );
    }

    /**
     * adjust cigar_ptr, n_cigar, read_length to ignore clipping
     *
     * @param n_cigar number of cigar operations passed by reference (adjusted by call)
     * @param cigar_ptr pointer to array of cigar operations (adjusted by call)
     * @param read_length length of read passed by reference (adjusted by call)
     *
     */
    void _length_adjust_cigar_soft_clipping(
            int &n_cigar, uint32_t** cigar_ptr, int &read_length
    ) {
        uint32_t* base_cigar = *cigar_ptr;  // original cigar array
        int length_removed = 0;  // keep track of how many bases removed
        // remove clipping on the right
        for (int i = n_cigar - 1; i >= 0; --i) {  // iterate backwards
            const char cigar_op = bam_cigar_op(base_cigar[i]);
            // ignore clipping operations and contribution to read length
            if (cigar_op == BAM_CHARD_CLIP) {
                --n_cigar;
            } else if (cigar_op == BAM_CSOFT_CLIP) {
                --n_cigar;
                read_length -= bam_cigar_oplen(base_cigar[i]);
            } else {
                break;
            }
        }
        // remove clipping on the left
        int lhs_clipping = 0;  // offset to apply to *cigar_ptr
        for (int i = 0; i < n_cigar; ++i) {
            const char cigar_op = bam_cigar_op(base_cigar[i]);
            // ignore clipping operations and contribution to read length
            if (cigar_op == BAM_CHARD_CLIP) {
                ++lhs_clipping;
            } else if (cigar_op == BAM_CSOFT_CLIP) {
                ++lhs_clipping;
                read_length -= bam_cigar_oplen(base_cigar[i]);
            } else {
                break;
            }
        }
        if (lhs_clipping > 0) {
            n_cigar -= lhs_clipping;
            *cigar_ptr = base_cigar + lhs_clipping;
        }
        return;
    }

    void IOBam::find_junction_genes(string chrom, char strand, int start, int end,
                                    float* nreads_ptr){
//        const int n = glist_[chrom].size() ;
        bool found_stage1 = false ;
        bool found_stage2 = false ;
        vector<Gene*> temp_vec1 ;
        vector<Gene*> temp_vec2 ;
        Junction * junc = new Junction(start, end, false, simpl_) ;
        const string key = junc->get_key() ;

        vector<overGene*>::iterator low = lower_bound (glist_[chrom].begin(), glist_[chrom].end(), start, _Region::func_comp ) ;
        if (low == glist_[chrom].end())
            return ;

        for (const auto &gObj: (*low)->glist){
//            if (gObj->get_start() >= end) break ;
            if (gObj->get_start() >= end) continue ;
            if (gObj->get_end() < start) continue ;
            if (start< gObj->get_start() || end> gObj->get_end()) continue ;
            if (strand == '.' || strand == gObj->get_strand()) {
                if(gObj->junc_map_.count(key) >0 ){
                    found_stage1 = true ;
                    gObj->initialize_junction(key, start, end, nreads_ptr, simpl_) ;
                } else if(found_stage1){
                    continue ;
                } else {
                    for(const auto &ex: gObj->exon_map_){
                        if (((ex.second)->get_start()-500) <= start && ((ex.second)->get_end()+500) >= end){
                            found_stage2 = true ;
                            temp_vec1.push_back(gObj) ;
                        }
                        else temp_vec2.push_back(gObj) ;
                    }
                }
            }
        }
        delete junc ;
        if (!found_stage1){
            if (found_stage2){
                for(const auto &g: temp_vec1){
                    g->initialize_junction(key, start, end, nreads_ptr, simpl_) ;
                }
            }else{
                for(const auto &g: temp_vec2){
                    g->initialize_junction(key, start, end, nreads_ptr, simpl_) ;
                }
            }
        }
        return ;
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
    void IOBam::add_junction(string chrom, char strand, int start, int end, const int offset, int sreads) {
        if (offset >= eff_len_ || offset < 0) {
            // XXX -- should warn that we are losing reads (should handle offsets
            // earlier, this function should not need to check a second time
            return;
        }
        // get string-based key to map over junctions (index of 2d coverage vector)
        string key = chrom + ":" + strand + ":" + to_string(start) + "-" + to_string(end) ;

        // indicator if this is the first time we have seen this junction
        bool new_j = false ;
        // pointer to junction coverage per position/offset
        float * v ;
        // acquire lock for working with junc_map/junc_vec
        omp_set_lock(&map_lck_) ;
        {
            if (junc_map.count(key) == 0 ) {
                // we haven't seen this junction before
                junc_map[key] = junc_vec.size();  // index to add to junc_vec
                // allocate array of floating point zeros
                v = (float*) calloc(eff_len_, sizeof(float)) ;
                // add array to junc_vec at position junc_map[key]
                junc_vec.push_back(v) ;
                // mark as new junction
                new_j = true ;
            } else {
                // get previously created coverage array for the junction
                v = junc_vec[junc_map[key]] ;
            }
        }
        omp_unset_lock(&map_lck_) ;

        if (new_j) {
            // associate junction with genes (prioritizing annotated, other stages for denovo)
            find_junction_genes(chrom, strand, start, end, v) ;
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

        // get length of read
        int read_length = read->core.l_qseq;

        // get cigar operations, adjust them and read length to ignore clipping
        int n_cigar = read->core.n_cigar;
        uint32_t *cigar = bam_get_cigar(read) ;
        _length_adjust_cigar_soft_clipping(n_cigar, &cigar, read_length);

        // track offsets on genomic and read position for moving along alignment
        int genomic_offset = 0;
        int read_offset = 0;

        // iterate over cigar operations, look for split reads
        for (int i = 0; i < n_cigar; ++i) {
            // extract operation and offset for the operation from cigar[i]
            const char cigar_op = bam_cigar_op(cigar[i]);
            const int cigar_oplen = bam_cigar_oplen(cigar[i]);
            if (cigar_op == BAM_CREF_SKIP) {
                if (read_offset < MIN_BP_OVERLAP) {
                    // too early in read, not enough overlap yet
                    continue;
                }
                // if we make it here, this is in a valid position to add
                // get start and end of junction
                const int junction_start = read_pos + genomic_offset;
                const int junction_end = junction_start + cigar_oplen + 1;
                // get junction position
                const int junction_pos = read_offset - MIN_BP_OVERLAP;
                // add the junction at the specified position
                try {
                    add_junction(chrom, _get_strand(read), junction_start,
                                 junction_end, junction_pos, 1);
                } catch (const std::logic_error& e) {
                    cout << "ERROR" << e.what() << '\n';
                }
                // finally, update the genomic offset
                genomic_offset += cigar_oplen;
            } else {
                // handle all other operations
                if (_cigar_advance_read(cigar_op)) {
                    read_offset += cigar_oplen;
                    if (read_offset > read_length - MIN_BP_OVERLAP) {
                        break;  // future operations will not lead to valid junction
                    }
                }
                if (_cigar_advance_reference(cigar_op)) {
                    genomic_offset += cigar_oplen;
                }
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

        // get iterator over all introns with end position after alignment start
        vector<Intron *>::iterator low = lower_bound(intronVec_[chrom].begin(), intronVec_[chrom].end(),
                                                     read_pos, _Region::func_comp ) ;
        if (low == intronVec_[chrom].end()) {
            // this read can't contribute to any of the introns, so exit here
            return 0;
        }

        // fill records of genomic/read offsets to determine offset of introns
        vector<pair<int, int>> genomic_read_offsets;
        genomic_read_offsets.push_back({read_pos, 0});
        // fill records of junctions splitting read (ignore introns if intersect)
        vector<pair<int, int>> junc_record ;
        // first and last genomic positions an intron can include to be considered
        int first_pos = -10;  // -10 is placeholder value
        int last_pos = -10;

        // get length of read
        int read_length = read->core.l_qseq;

        // get cigar operations, adjust them and read length to ignore clipping
        int n_cigar = read->core.n_cigar ;
        uint32_t *cigar = bam_get_cigar(read) ;
        _length_adjust_cigar_soft_clipping(n_cigar, &cigar, read_length);

        // parse CIGAR operations for genomic_read_offsets, junc_record, {first,last}_pos
        int genomic_offset = 0;
        int read_offset = 0;
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
            const bool advance_read = _cigar_advance_read(cigar_op);
            const bool advance_reference = _cigar_advance_reference(cigar_op);
            if (advance_read) {
                // check if this operation gives us first/last position for valid intron
                if (
                        (read_offset <= MIN_BP_OVERLAP)
                        && ((read_offset + cigar_oplen) > MIN_BP_OVERLAP)
                ) {
                    // + 1 magic number for 0-->1-based indexing for intron start
                    first_pos = read_pos + genomic_offset + 1
                        + (advance_reference ? MIN_BP_OVERLAP - read_offset : 0);
                }
                if (
                        ((read_offset + cigar_oplen) >= read_length - MIN_BP_OVERLAP)
                        && (read_offset < read_length - MIN_BP_OVERLAP)
                ) {
                    last_pos = read_pos + genomic_offset
                        + (advance_reference ? read_length - MIN_BP_OVERLAP - read_offset : 0);
                }
                read_offset += cigar_oplen;
            }
            if (advance_reference) {
                genomic_offset += cigar_oplen;
            }
            if (advance_reference != advance_read) {
                // only need update if one advanced but not the other
                const int current_genomic_pos = read_pos + genomic_offset;
                const int current_offset = read_offset;
                genomic_read_offsets.push_back({current_genomic_pos, current_offset});
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
                if (j.first < intron->get_end() && intron->get_start() < j.second) {
                    junc_found = true;
                    break;
                } else if (j.second >= intron->get_end()) {
                    // none of the remaining junctions (sorted) can intersect intron
                    break;
                }
            }
            if (junc_found) {
                continue;  // read not compatible with intron due to junctions
            }

            // since read includes intron within acceptable bounds, add intron
            // with appropriate offset
            // calculate relative to last acceptable position in read, which has offset read_length - MIN_BP_OVERLAP
            int relative_offset = 0;
            int read_start_vs_intron_start = read_pos - intron->get_start();
            if (read_start_vs_intron_start >= 0) {
                // read start at or after intron start, further away from end of read
                relative_offset = read_start_vs_intron_start;
            } else {
                int i = 1;  // get index of first offset past intron start
                for (; i < genomic_read_offsets.size(); ++i) {
                    // could replace with something like std::lower_bound for log(n) search
                    // but premature optimization unless many cigar operations
                    if (genomic_read_offsets[i].first >= intron->get_start()) {
                        break;
                    }
                }
                relative_offset =
                    (genomic_read_offsets[i - 1].first - intron->get_start())  // how far this genomic coordinate was from intron start
                    - genomic_read_offsets[i - 1].second;  // adjust relative offset on read
            }
            int intron_offset = read_length - MIN_BP_OVERLAP + relative_offset;
            #pragma omp critical
            {
                intron->add_read(intron_offset, eff_len_, 1);
            }
        }
        return 0 ;
    }


    int IOBam::ParseJunctionsFromFile(bool ir_func){

        samFile *in;
        int ignore_sam_err = 0;
        bam_hdr_t *header ;
        bam1_t *aln;

        int r = 0, exit_code = 0;
//        hts_opt *in_opts = NULL;
        int extra_hdr_nuls = 0;

        in = sam_open(bam_.c_str(), "rb") ;
        if (NULL == in) {
            fprintf(stderr, "Error opening \"%s\"\n", bam_.c_str());
            return EXIT_FAILURE;
        }
        header = sam_hdr_read(in);
        if (NULL == header) {
            fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_.c_str());
            return EXIT_FAILURE;
        }
        header->ignore_sam_err = ignore_sam_err;
        if (extra_hdr_nuls) {
            char *new_text = (char*) realloc(header->text, header->l_text + extra_hdr_nuls);
            if (new_text == NULL) {
                fprintf(stderr, "Error reallocing header text\n");
                return EXIT_FAILURE;
            }
            header->text = new_text;
            memset(&header->text[header->l_text], 0, extra_hdr_nuls);
            header->l_text += extra_hdr_nuls;
        }

        aln = bam_init1();
        htsThreadPool p = {NULL, 0};
        if (nthreads_ > 0) {
            p.pool = hts_tpool_init(nthreads_);
            if (!p.pool) {
                fprintf(stderr, "Error creating thread pool\n");
                exit_code = 1;
            } else {
                hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
            }
        }

        int (IOBam::*parse_func)(bam_hdr_t *, bam1_t *) ;
        if(ir_func) parse_func = &IOBam::parse_read_for_ir ;
        else parse_func = &IOBam::parse_read_into_junctions ;


        while ((r = sam_read1(in, header, aln)) >= 0) {
            (this->*parse_func)(header, aln) ;
        }
        if (r < -1) {
            fprintf(stderr, "Error parsing input.\n");
            exit_code = 1;
        }
        bam_destroy1(aln);
        bam_hdr_destroy(header);

        r = sam_close(in);
        if(!ir_func)
            junc_limit_index_ = junc_vec.size() ;
        if (r < 0) {
            fprintf(stderr, "Error closing input.\n");
            exit_code = 1 ;
        }
        if (p.pool)
            hts_tpool_destroy(p.pool);

        return exit_code ;
    }

    int IOBam::normalize_stacks(vector<float> vec, float sreads, int npos, float fitfunc_r, float pvalue_limit){
        if (fitfunc_r == 0.0){
            for (int i=0; i<(int)vec.size(); i++){
                const float mean_reads = (npos == 1) ? 0.5: (sreads-vec[i]) / (npos - 1) ;
                const float pvalue = 1 - scythe::ppois(vec[i], mean_reads) ;
                if (pvalue< pvalue_limit){
                    vec.erase(vec.begin() + i) ;
                    npos -- ;
                }
            }
        }else{
            for (int i=0; i<(int) vec.size(); i++){
                const float mean_reads = (npos == 1) ? 0.5: (sreads-vec[i]) / (npos - 1) ;
                const float r = 1 / fitfunc_r ;
                const float p = r/(mean_reads + r) ;
                const float pvalue = 1 - scythe::pnbinom(vec[i], r, p) ;
                if (pvalue< pvalue_limit){
                    vec.erase(vec.begin() + i);
                    npos -- ;
                }
            }
        }
        return npos ;
    }

    int IOBam::boostrap_samples(int msamples, int ksamples, float* boots, float fitfunc_r, float pvalue_limit){
        const int njunc = junc_map.size();

        #pragma omp parallel for num_threads(nthreads_)
        for(int jidx=0; jidx < njunc; jidx++){

            vector<float> vec ;
            int npos = 0 ;
            float sreads = 0 ;
            for(unsigned int i=0; i<eff_len_; ++i){
                if (junc_vec[jidx][i]>0){
                    npos ++ ;
                    sreads += junc_vec[jidx][i] ;
                    vec.push_back(junc_vec[jidx][i]) ;
                }
            }
            if (npos == 0) continue ;
            if (pvalue_limit > 0) npos = normalize_stacks(vec, sreads, npos, fitfunc_r, pvalue_limit) ;
            if (npos == 0) continue ;
            default_random_engine generator;
            uniform_int_distribution<int> distribution(0, npos-1);
            for (int m=0; m<msamples; m++){
                float lambda = 0;
                for (int k=0; k<ksamples; k++)lambda += vec[distribution(generator)] ;
                lambda /= ksamples ;
                const int idx2d = (jidx*msamples) + m ;
                boots[idx2d] = (lambda * npos) ;
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
            float ss = 0 ;
            int np = 0 ;
            for(unsigned int i=0; i<eff_len_; ++i){
                ss += junc_vec[jidx][i] ;
                np += (junc_vec[jidx][i]>0)? 1 : 0 ;
            }
            res[jidx] = (int) ss ;
            res[jidx+njunc] = np ;
        }
        return res ;
    }


    void IOBam::parseJuncEntry(map<string, vector<overGene*>> & glist, string gid, string chrom, char strand,
                               int start, int end, unsigned int sreads, unsigned int minreads_t, unsigned int npos,
                               unsigned int minpos_t, unsigned int denovo_t, bool denovo, vector<Gene*>& oGeneList,
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
                        junc_vec.push_back(irptr->read_rates_) ;
                        gObj->add_intron(irptr, min_intron_cov, minexp, min_bins, reset) ;
                    }
                }
                omp_unset_lock(&map_lck_) ;
            }
        } else {
            string key = to_string(start) + "-" + to_string(end) ;
            add_junction(chrom, strand, start, end, 0, sreads) ;
            for (const auto &gObj: (*low)->glist){
                gObj->updateFlagsFromJunc(key, sreads, minreads_t, npos, minpos_t, denovo_t, denovo, minexp, reset) ;
            }

        }
        return ;

    }


    void IOBam::detect_introns(float min_intron_cov, unsigned int min_experiments, float min_bins, bool reset){
        for (const auto & it: glist_){
            if (intronVec_.count(it.first)==0){
                intronVec_[it.first] = vector<Intron*>() ;
            }
            const int n = (it.second).size() ;
            #pragma omp parallel for num_threads(nthreads_)
            for(int g_it = 0; g_it<n; g_it++){
                for (const auto &g: ((it.second)[g_it])->glist){
                    g->detect_exons() ;
                    g->detect_introns(intronVec_[it.first], simpl_) ;
                    g->reset_exons() ;
                }
            }
            sort(intronVec_[it.first].begin(), intronVec_[it.first].end(), Intron::islowerRegion<Intron>) ;

        }
        ParseJunctionsFromFile(true) ;
        for (const auto & it: intronVec_){
            const int n = (it.second).size() ;

            #pragma omp parallel for num_threads(nthreads_)
            for(int idx = 0; idx<n; idx++){
                Intron * intrn_it = (it.second)[idx] ;
                const bool pass = intrn_it->is_reliable(min_bins, eff_len_) ;
                if( pass ){
                    const string key = intrn_it->get_key(intrn_it->get_gene()) ;
                    #pragma omp critical
                    {
                        if (junc_map.count(key) == 0) {
                            junc_map[key] = junc_vec.size() ;
                            junc_vec.push_back(intrn_it->read_rates_) ;
                            (intrn_it->get_gene())->add_intron(intrn_it, min_intron_cov, min_experiments, min_bins, reset) ;
                        }
                    }
                }else{
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
            const unsigned int i = idx - junc_limit_index_ ;
            for (unsigned int j=0; j<eff_len_; j++){
                const unsigned int indx_2d = i*eff_len_ + j ;
                out_cov[indx_2d] = junc_vec[idx][j] ;
            }
        }
    }

    void prepare_genelist(map<string, Gene*>& gene_map, map<string, vector<overGene*>> & geneList){
        map<string, vector<Gene*>> gmp ;

        for(const auto & p: gene_map){
            const string chrom = (p.second)->get_chromosome() ;
            if (gmp.count(chrom) == 0 )
                gmp[chrom] = vector<Gene*>() ;
            gmp[chrom].push_back(p.second)  ;
        }

        for (auto &gl: gmp){
            sort((gl.second).begin(), (gl.second).end(), Gene::islowerRegion<Gene>) ;
            int gstart = 0, gend = 0 ;
            overGene * ov = nullptr ;

            for (const auto &g: (gl.second)){

                int st = g->get_start() ;
                int nd = g->get_end() ;
                if (st <= gend && nd >= gstart){
                    gstart = (gstart == 0) ? st : min(gstart, st) ;
                    gend = (gend == 0) ? nd : max(gend, nd) ;
                    ov->set_start(gstart) ;
                    ov->set_end(gend) ;
                }
                else{
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
