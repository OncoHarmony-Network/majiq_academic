#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <random>
#include <string>
#include <cmath>
#include "io_bam.hpp"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"

using namespace std;
//using namespace interval;
namespace io_bam {

//    struct CompareIntervals
//    {
//       pair<int, int> asInterval(const Gene *g) const // or static
//       {
//          return {g->start, g->end};
//       }
//
//       pair<int, int> asInterval(Junction * j) const // or static
//       {
//            return {j->start, j->end};
//       }
//
//       template< typename T1, typename T2 >
//       bool operator()( T1 const t1, T2 const t2 ) const
//       {
//           pair<int,int> p1 = asInterval(t1) ;
//           pair<int,int> p2 = asInterval(t2) ;
//           return ((p1.first <= p2.second) && (p1.second >= p2.first)) ;
//       }
//    };

    bool juncinGene(Gene* t1, Junction* t2){
        return ( (t1->get_start() <= t2->get_end()) && (t1->get_end() >= t2->get_start()) ) ;
    }

    int juncGeneSearch(vector<Gene *> a, int n, Junction * j) {
        int l = 0 ;
        int h = n ; // Not n - 1
        while (l < h) {
            int mid = (l + h) / 2 ;
            if (juncinGene(a[mid], j)) {
                h = mid ;
            } else {
                l = mid + 1 ;
            }
        }
        return l;
    }

    char IOBam::_get_strand(bam1_t * read){
        char strn = '.';

        if (strandness_ == FWD_STRANDED){
            strn = (((read->core.flag & 0x10) == (0x00 & is_read1(read)))
                    || ((read->core.flag & 0x10) == (0x10 & is_read2(read)))) ? '+' : '-';

        } else if (strandness_ == FWD_STRANDED){
            strn = (((read->core.flag & 0x10) ==(0x10 & is_read1(read)))
                    || ((read->core.flag & 0x10) == (0x00 & is_read2(read)))) ? '+' : '-';
        }
        return (strn);
    }

    inline int _unique(bam1_t * read){
        return ((read->core.flag & 0x100) != 0x100);
    }

    void IOBam::find_junction_genes(string chrom, char strand, int start, int end,
                                    int* nreads_ptr ){
        const int n = glist_[chrom].size();
        bool found_stage1 = false ;
        bool found_stage2 = false ;
        vector<Gene*> temp_vec1 ;
        vector<Gene*> temp_vec2 ;
        Junction * junc = new Junction(start, end, false) ;
        const string key = junc->get_key() ;
        int i = juncGeneSearch(glist_[chrom], n, junc) ;
        while(i< n ){
            Gene * gObj  = glist_[chrom][i] ;
            ++i ;
            if (!juncinGene(gObj, junc)) break ;
            if (strand == '.' || strand == gObj->get_strand()) {
                if(gObj->junc_map_.count(key) >0 ){
                    found_stage1 = true ;
                    gObj->junc_map_[key]->set_nreads_ptr(nreads_ptr) ;
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
                    g->junc_map_[key] = new Junction(start, end, false) ;
                    g->junc_map_[key]->set_nreads_ptr(nreads_ptr) ;
                }
            }else{
                for(const auto &g: temp_vec2){
                    g->junc_map_[key] = new Junction(start, end, false) ;
                    g->junc_map_[key]->set_nreads_ptr(nreads_ptr) ;
                }
            }
        }
        return ;
    }

    inline void IOBam::update_junction_read(string key, int offset, int count) {
        int* vec = junc_vec[junc_map[key]] ;
        vec[offset] += count ;
        return ;
    }


    void IOBam::add_junction(string chrom, char strand, int start, int end, int read_pos) {

        const unsigned int offset = start - (read_pos+ MIN_BP_OVERLAP) ;
        if (offset >= eff_len_) return ;

        string key = chrom + ":" + strand + ":" + to_string(start) + "-" + to_string(end) ;
        #pragma omp critical
        {
            if(junc_map.count(key) == 0) {
                junc_map[key] = junc_vec.size() ;
                int* v = (int*) calloc(eff_len_, sizeof(int)) ;
                junc_vec.push_back(v) ;
                find_junction_genes(chrom, strand, start, end, v) ;
            }
            update_junction_read(key, offset, 1) ;
        }
        return ;
    }


    int IOBam::parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) {
        int n_cigar = read->core.n_cigar;
        if (n_cigar <= 1 || !_unique(read)) // max one cigar operation exists(likely all matches)
            return 0;

        int read_pos = read->core.pos;
        string chrom(header->target_name[read->core.tid]);

        uint32_t *cigar = bam_get_cigar(read);

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
                    add_junction(chrom, _get_strand(read), j_start, j_end, read_pos);
                } catch (const std::logic_error& e) {
                    cout << "KKK" << e.what() << '\n';
                }
                off += ol ;
            }
        }
        return 0;
    }


    int IOBam::ParseJunctionsFromFile(){

        samFile *in;
        int flag = 0, ignore_sam_err = 0;
        char moder[8];
        bam_hdr_t *header ;
        bam1_t *aln;

        int r = 0, exit_code = 0;
        hts_opt *in_opts = NULL;
        int nreads = 0;
        int extra_hdr_nuls = 0;
        int multi_reg = 0;


    //    strcpy(moder, "r");
    //    if (flag & READ_CRAM) strcat(moder, "c");
    //    else if ((flag & READ_COMPRESSED) == 0) strcat(moder, "b");
    //    in = sam_open(filename, moder) ;
        cout << "FILE:" << bam_ << "\n";
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

        while ((r = sam_read1(in, header, aln)) >= 0) {
//            cout<<"read" << aln->core.tid << " :: " << aln->core.pos << "\n" ;
            parse_read_into_junctions(header, aln) ;
            if (nreads && --nreads == 0)
                break;
        }
        cout<< "FINISH FILE\n" ;
        // Finishing

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


    int IOBam::boostrap_samples(int msamples, int ksamples, float* boots){

        float * p = boots ;
        const int njunc = junc_map.size();

        #pragma omp parallel for num_threads(nthreads_)
        for(int jidx=0; jidx < njunc; jidx++){

            int* vec = junc_vec[jidx] ;
            int npos = 0 ;
            for(int i=0; i<eff_len_; ++i){
                npos = vec[i]? 1 : 0 ;
            }
            if (npos == 0) continue ;
            default_random_engine generator;
            uniform_int_distribution<int> distribution(0, npos-1);
            for (int m=0; m<msamples; m++){
                float lambda = 0;
                for (int k=0; k<ksamples; k++)lambda += vec[distribution(generator)] ;
                lambda /= ksamples ;
                *p = lambda * npos ;
                p++ ;
            }
        }
        return 0 ;
    }

    int IOBam::get_njuncs(){
        return junc_map.size() ;
    }

    const map<string, unsigned int>& IOBam::get_junc_map(){
        return junc_map ;
    }

    int * IOBam::get_junc_vec_summary(){
        int njunc = junc_map.size() ;
        int * res = (int *) calloc(2*njunc, sizeof(int)) ;

        #pragma omp parallel for num_threads(nthreads_)
        for(int jidx=0; jidx < njunc; ++jidx){
            int ss = 0 ;
            int np = 0 ;
            for(int i=0; i<eff_len_; ++i){
                ss += junc_vec[jidx][i] ;
                np += junc_vec[jidx][i]? 1 : 0 ;
            }
            res[jidx] = ss ;
            res[jidx+njunc] = np ;
        }
        return res ;
    }

}