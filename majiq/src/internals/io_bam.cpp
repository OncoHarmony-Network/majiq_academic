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


// initialize static random state
vector<default_random_engine> io_bam::IOBam::generators_ = vector<default_random_engine>(0);

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

    void IOBam::find_junction_genes(string chrom, char strand, int start, int end,
                                    shared_ptr<vector<float>> nreads_ptr){
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


    void IOBam::add_junction(string chrom, char strand, int start, int end, int read_pos, int first_offpos, int sreads) {

        const unsigned int offset = (first_offpos == -1) ? (start - (read_pos+ MIN_BP_OVERLAP)) : first_offpos ;
        if (offset >= eff_len_) return ;
        string key = chrom + ":" + strand + ":" + to_string(start) + "-" + to_string(end) ;

        bool new_j = false ;
        omp_set_lock(&map_lck_) ;
        {
            if (junc_map.count(key) == 0 ) {
                junc_map[key] = junc_vec.size() ;
                junc_vec.push_back(make_shared<vector<float>>(buff_len_, 0));
                // this is a new junction
                new_j = true ;
            }
        }
        omp_unset_lock(&map_lck_) ;
        shared_ptr<vector<float>> v_ptr = junc_vec[junc_map[key]];
        vector<float> &v = *v_ptr;

        if (new_j) {
            find_junction_genes(chrom, strand, start, end, v_ptr) ;
        }
        #pragma omp atomic
            v[offset] += sreads ;

        return ;
    }


    int IOBam::parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) {
        int n_cigar = read->core.n_cigar;
        if (!_unique(read) || _unmapped(read)) // max one cigar operation exists(likely all matches)
            return 0;

        int read_pos = read->core.pos;
        string chrom(header->target_name[read->core.tid]);

        uint32_t *cigar = bam_get_cigar(read) ;

        int off = 0 ;
        int first_offpos = -1 ;
        for (int i = 0; i < n_cigar; ++i) {
            const char op = bam_cigar_op(cigar[i]);
            const int ol = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CDEL ||  op == BAM_CEQUAL || op == BAM_CDIFF){
                 off += ol ;
            }
            else if( op == BAM_CREF_SKIP ){
                const int rlen = read->core.l_qseq ;
                if (off >= MIN_BP_OVERLAP && off<= (rlen - MIN_BP_OVERLAP)){
                    const int j_end = read->core.pos + off + ol +1 ;
                    const int j_start =  read->core.pos + off ;
                    first_offpos  = (first_offpos == -1) ? (j_start - (read_pos+ MIN_BP_OVERLAP)) : first_offpos ;

                    try {
                        add_junction(chrom, _get_strand(read), j_start, j_end, read_pos, first_offpos, 1) ;
                    } catch (const std::logic_error& e) {
                        cout << "ERROR" << e.what() << '\n';
                    }
                }
                off += ol ;
            }
        }
        return 0;
    }


    int IOBam::parse_read_for_ir(bam_hdr_t *header, bam1_t *read) {
        int n_cigar = read->core.n_cigar ;
        if (!_unique(read) || _unmapped(read)) // max one cigar operation exists(likely all matches)
            return 0;
        const int read_pos = read->core.pos;

        const string chrom(header->target_name[read->core.tid]) ;
        if (intronVec_.count(chrom) == 0)
             return 0 ;
        const int rlen = read->core.l_qseq ;
        vector<Intron *>::iterator low = lower_bound (intronVec_[chrom].begin(), intronVec_[chrom].end(),
                                                      read_pos, _Region::func_comp ) ;
        if (low ==  intronVec_[chrom].end()) return 0 ;
        vector<pair<int, int>> junc_record ;

        const char read_strand = _get_strand(read) ;
        int off = 0 ;
        uint32_t *cigar = bam_get_cigar(read) ;
        for (int i = 0; i < n_cigar; ++i) {
            const char op = bam_cigar_op(cigar[i]);
            const int ol = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CDEL ||  op == BAM_CEQUAL || op == BAM_CDIFF){
                 off += ol;
            }
            else if( op == BAM_CREF_SKIP && off >= MIN_BP_OVERLAP){
                const int j_end = read->core.pos + off + ol +1;
                const int j_start =  read->core.pos + off;
                try {
                    junc_record.push_back({j_start, j_end}) ;
                } catch (const std::logic_error& e) {
                    cout << "ERROR" << e.what() << '\n';
                }
                off += ol ;
            }
        }

        for (; low != intronVec_[chrom].end() ; low++){

            bool junc_found = false ;
            Intron * intron = *low;
            const char gstrand = intron->get_gene()->get_strand() ;

            if(intron->get_start()> read_pos+rlen) break ;
            if (intron->get_end() <= read_pos) continue ;
            if (read_strand == '.' || gstrand == read_strand){
                for (const auto & j:junc_record){
                    if ((j.first>=intron->get_start() && j.first<= intron->get_end() )
                        || (j.second>=intron->get_start() && j.second<= intron->get_end())){
                        junc_found = true ;
                        break ;
                    }
                }
                if (!junc_found){
                  #pragma omp critical
                        intron->add_read(read_pos, eff_len_, 1) ;
                }
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

    int IOBam::bootstrap_samples(int msamples, int ksamples, float* boots, float fitfunc_r, float pvalue_limit) {
        const int njunc = junc_map.size();

        #pragma omp parallel for num_threads(nthreads_)
        for(int jidx=0; jidx < njunc; jidx++){

            vector<float> &coverage = *(junc_vec[jidx]);
            vector<float> vec ;
            int npos = 0 ;
            float sreads = 0 ;
            for(unsigned int i=0; i<eff_len_; ++i){
                if (coverage[i] > 0) {
                    npos ++ ;
                    sreads += coverage[i];
                    vec.push_back(coverage[i]);
                }
            }
            if (npos == 0) continue ;
            if (pvalue_limit > 0) npos = normalize_stacks(vec, sreads, npos, fitfunc_r, pvalue_limit) ;
            if (npos == 0) continue ;
            default_random_engine &generator = generators_[omp_get_thread_num()];
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
                        junc_vec.push_back(irptr->read_rates_ptr_) ;
                        gObj->add_intron(irptr, min_intron_cov, minexp, min_bins, reset) ;
                    }
                }
                omp_unset_lock(&map_lck_) ;
            }
        } else {
            string key = to_string(start) + "-" + to_string(end) ;
            add_junction(chrom, strand, start, end, 0, 0, sreads) ;
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
                            (intrn_it->get_gene())->add_intron(intrn_it, min_intron_cov, min_experiments, min_bins, reset) ;
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
