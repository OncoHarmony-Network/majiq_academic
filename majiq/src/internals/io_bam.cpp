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

//    inline bool IOBam::new_junc_values(const string key){
//        bool res = false ;
//        if (junc_map.count(key) == 0 ) {
//            junc_map[key] = junc_vec.size() ;
//            float * v = (float*) calloc(eff_len_, sizeof(float)) ;
//            junc_vec.push_back(v) ;
//        }
//        return res;
//    }

    char IOBam::_get_strand(bam1_t * read){
        char strn = '.';

        if (strandness_ == FWD_STRANDED){
            strn = (((read->core.flag & 0x10) == (0x00 & is_read1(read)))
                    || ((read->core.flag & 0x10) == (0x10 & is_read2(read)))) ? '+' : '-' ;

        } else if (strandness_ == REV_STRANDED){
            strn = (((read->core.flag & 0x10) ==(0x10 & is_read1(read)))
                    || ((read->core.flag & 0x10) == (0x00 & is_read2(read)))) ? '+' : '-' ;
        }
//cout << strn << " " << read->core.flag << " " << read->id<<"\n" ;

        return (strn);
    }

    inline int _unique(bam1_t * read){
        return ((read->core.flag & 0x100) != 0x100) ;
    }

    void IOBam::find_junction_genes(string chrom, char strand, int start, int end,
                                    float* nreads_ptr){
        const int n = glist_[chrom].size();
        bool found_stage1 = false ;
        bool found_stage2 = false ;
        vector<Gene*> temp_vec1 ;
        vector<Gene*> temp_vec2 ;
        Junction * junc = new Junction(start, end, false) ;
        const string key = junc->get_key() ;
        int i = Gene::RegionSearch(glist_[chrom], n, end) ;

        if(i<0) return ;
        while(i< n){
            Gene * gObj  = glist_[chrom][i] ;
            ++i ;
            if (gObj->get_start() >= end) break ;
            if (gObj->get_end() < start) continue ;
//            if (start< gObj->get_start() || end> gObj->get_end()) continue ;
            if (strand == '.' || strand == gObj->get_strand()) {
                if(gObj->junc_map_.count(key) >0 ){
                    found_stage1 = true ;
                    gObj->initialize_junction(key, start, end, nreads_ptr) ;
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
                    g->initialize_junction(key, start, end, nreads_ptr) ;
                }
            }else{
                for(const auto &g: temp_vec2){
                    g->initialize_junction(key, start, end, nreads_ptr) ;
                }
            }
        }
        return ;
    }


    void IOBam::add_junction(string chrom, char strand, int start, int end, int read_pos, int first_offpos) {

        const unsigned int offset = (first_offpos == -1) ? (start - (read_pos+ MIN_BP_OVERLAP)) : first_offpos ;
        if (offset >= eff_len_) return ;
        string key = chrom + ":" + strand + ":" + to_string(start) + "-" + to_string(end) ;

        bool new_j = false ;
        float * v ;
        #pragma omp critical
        {
            if (junc_map.count(key) == 0 ) {
                junc_map[key] = junc_vec.size() ;
                v = (float*) calloc(eff_len_, sizeof(float)) ;
                junc_vec.push_back(v) ;
                new_j = true ;
            } else {
                v = junc_vec[junc_map[key]] ;
            }
        }
        if (new_j) {
            find_junction_genes(chrom, strand, start, end, v) ;
        }
        #pragma omp atomic
            junc_vec[junc_map[key]][offset] += 1 ;

        return ;
    }


    int IOBam::parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) {
        int n_cigar = read->core.n_cigar;
        if (n_cigar <= 1 || !_unique(read)) // max one cigar operation exists(likely all matches)
            return 0;

        int read_pos = read->core.pos;
        string chrom(header->target_name[read->core.tid]);

        uint32_t *cigar = bam_get_cigar(read) ;

        int off = 0 ;
        int first_offpos = -1 ;
//cout << "read id: " << bam_get_qname(read) << "\n" ;
        for (int i = 0; i < n_cigar; ++i) {
            const char op = bam_cigar_op(cigar[i]);
            const int ol = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CDEL ||  op == BAM_CEQUAL || op == BAM_CDIFF){
                 off += ol ;
            }
            else if( op == BAM_CREF_SKIP ){
                if (off >= MIN_BP_OVERLAP){
                    const int j_end = read->core.pos + off + ol +1 ;
                    const int j_start =  read->core.pos + off ;
                    first_offpos  = (first_offpos == -1) ? (j_start - (read_pos+ MIN_BP_OVERLAP)) : first_offpos ;

//cout << "add junction " << read_pos << ":: "<< j_start << "-" << j_end << "\n" ;
                    try {
                        add_junction(chrom, _get_strand(read), j_start, j_end, read_pos, first_offpos) ;
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
//cerr << "R1: read id: " << bam_get_qname(read) << " : " << n_cigar<< "\n" ;
        if (n_cigar <= 1 || !_unique(read)) // max one cigar operation exists(likely all matches)
            return 0;
        const int read_pos = read->core.pos;
        const string chrom(header->target_name[read->core.tid]) ;
//cerr << "R2 " << read_pos << " :: " << chrom << "\n" ;

        if (intronVec_.count(chrom) == 0)
             return 0 ;
//cerr << "R3\n" ;
        const int  nintrons = intronVec_[chrom].size() ;
        uint32_t *cigar = bam_get_cigar(read) ;
//cerr << "R4\n" ;

        int idx = _Region::RegionSearch(intronVec_[chrom], nintrons, read_pos) ;
//cerr << "R5" << idx << "\n" ;
        if (idx<0) return 0 ;
//cerr << "R6 \n" ;
        vector<pair<int, int>> junc_record ;

        int off = 0;
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
        for(int i=idx; i<nintrons; ++i){
            bool junc_found = false ;
            Intron * intron = intronVec_[chrom][i] ;
            if(intron->get_start()> read_pos) break ;
            for (const auto & j:junc_record){
                if ((j.first>=intron->get_start() && j.first<= intron->get_end() )
                    || (j.second>=intron->get_start() && j.second<= intron->get_end())){
                    junc_found = true ;
                    break ;
                }
            }
            if (!junc_found){
                #pragma omp critical
                {
                    if (intron->read_rates_ == nullptr){
                        intron->add_read_rates_buff(eff_len_) ;
                    }
                    const int offset = (int)(read_pos - intron->get_start()) % intron->nbins_ ;
                    intron->read_rates_[offset] += 1 ;
                }
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
        if (r < 0) {
            fprintf(stderr, "Error closing input.\n");
            exit_code = 1;
        }
        if (p.pool)
            hts_tpool_destroy(p.pool);

        return exit_code;
    }

    int IOBam::normalize_stacks(vector<float> vec, float sreads, int npos, float fitfunc_r, float pvalue_limit){
//cout << "iniT norm stacks " << sreads<< ": " << npos <<": " << pvalue_limit<< "\n" ;
        const float mean_reads = sreads/npos ;
        if (fitfunc_r == 0.0){
            for (int i=0; i<(int)vec.size(); i++){
                const float pvalue = 1 - scythe::ppois(vec[i], mean_reads) ;
                if (pvalue< pvalue_limit){
                    vec.erase(vec.begin() + i) ;
                    npos -- ;
                }
            }

        }else{
            for (int i=0; i<(int) vec.size(); i++){
                const float r = 1 / fitfunc_r ;
                const float p = r/(mean_reads + r) ;
                const float pvalue = 1 - scythe::pnbinom(vec[i], r, p) ;
                if (pvalue< pvalue_limit){
                    vec.erase(vec.begin() + i);
                    npos -- ;
                }
            }
        }
//cout << "OUT norm stacks " << npos << "\n" ;
        return npos ;
    }

    int IOBam::boostrap_samples(int msamples, int ksamples, float* boots, float fitfunc_r, float pvalue_limit){

//        float * p = boots ;
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
//                *p = (lambda * npos) ;
//                p++ ;
            }
        }
//cout << "OUT BOOTSTRAP\n" ;
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

    void IOBam::detect_introns(float min_intron_cov, unsigned int min_experiments, float min_bins, bool reset){

        junc_limit_index_ = junc_vec.size() ;

        for (const auto & it: glist_){
            if (intronVec_.count(it.first)==0){
                intronVec_[it.first] = vector<Intron*>() ;
            }
            const int n = (it.second).size() ;
            #pragma omp parallel for num_threads(nthreads_)
            for(int g_it = 0; g_it<n; g_it++){
                (it.second)[g_it]->detect_introns(intronVec_[it.first]) ;
            }
            sort(intronVec_[it.first].begin(), intronVec_[it.first].end(), Intron::islowerRegion<Intron>) ;
        }

        ParseJunctionsFromFile(true) ;
//cout << "Num junctions: " << junc_vec.size() << "\n" ;
        for (const auto & it: intronVec_){
            const int n = (it.second).size() ;
            #pragma omp parallel for num_threads(nthreads_)
            for(int idx = 0; idx<n; idx++){
                Intron * intrn_it = (it.second)[idx] ;
                const bool pass = intrn_it->is_reliable(min_bins) ;
                if( pass ){
                    const string key = intrn_it->get_key(intrn_it->get_gene()) ;
                    #pragma omp critical
                    {
                        if (junc_map.count(key) == 0) {
                            junc_map[key] = junc_vec.size() ;
                            junc_vec.push_back(intrn_it->read_rates_) ;
                            (intrn_it->get_gene())->add_intron(intrn_it, min_intron_cov, min_experiments, reset) ;
                        }
                    }
                }else{
                    intrn_it->free_nreads() ;
                    delete intrn_it ;
                }
            }
        }
    }

//    void IOBam::fit_nb(vector<float*> jlist, ){
//
//        vector<float> pvalues ;
//
//
//    }


}
