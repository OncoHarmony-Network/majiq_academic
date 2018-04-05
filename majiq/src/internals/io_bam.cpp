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

    struct CompareIntervals
    {
       pair<int, int> asInterval(const Gene *g) const // or static
       {
          return {g->start, g->end};
       }

       pair<int, int> asInterval(Junction * j) const // or static
       {
            return {j->start, j->end};
       }

       template< typename T1, typename T2 >
       bool operator()( T1 const t1, T2 const t2 ) const
       {
           pair<int,int> p1 = asInterval(t1) ;
           pair<int,int> p2 = asInterval(t2) ;
           return ((p1.first <= p2.second) && (p1.second >= p2.first)) ;
       }
    };


    bool juncGeneSearch(Gene* t1, Junction* t2){
        return ((t1->start <= t2->end) && (t1->end >= t2->start)) ;
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



    void IOBam::find_junction_gene(string chrom, char strand, Junction * junc){
        vector<Gene *>::iterator gObjIt ;
        gObjIt = lower_bound(glist_[chrom].begin(), glist_[chrom].end(), junc, juncGeneSearch) ;
        bool found = false ;
        const string key = junc->get_key() ;
        vector<Gene*> temp_vec ;

        while(gObjIt != glist_[chrom].end()){
            if (!juncGeneSearch(*gObjIt, junc)) break ;
            if (strand == '.' || strand == (*gObjIt)->strand) {
                for(const auto &ex: (*gObjIt)->exon_map){

                    if (((ex.second)->start-500) <= junc->start && ((ex.second)->end+500) >= junc->start){
                        found = true ;
                        #pragma omp critical
                        (*gObjIt)->junc_map[key] = junc ;

                    }
                    else temp_vec.push_back((*gObjIt)) ;
                }
            }
            gObjIt++ ;
        }

        if (found) return ;
        for(const auto &g: temp_vec){
            #pragma omp critical
            g->junc_map[key] = junc ;
        }

        return ;
    }



    void IOBam::add_junction(string chrom, char strand, int start, int end, int read_pos) {

        const unsigned int offset = start - (read_pos+ MIN_BP_OVERLAP) ;
        bool newjunc = false ;

        if (offset >= eff_len_) {
            return ;
        }
        string key =  chrom + ":" + to_string(start) + "-" + to_string(end) ;
        #pragma omp critical
        {
            if(junc_map.count(key) == 0) {
                junc_map[key] = junc_vec.size() ;

                Junction * v = new Junction(start, end, false) ;
                junc_vec.push_back(v) ;
                newjunc = true ;
            }

            junc_vec[junc_map[key]]->update_junction_read(read_pos, 1) ;
        }

        find_junction_gene( chrom,  strand, junc_vec[junc_map[key]]) ;

        return;
    }


    int IOBam::parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) {
        int n_cigar = read->core.n_cigar;
        if (n_cigar <= 1 || !_unique(read)) // max one cigar operation exists(likely all matches)
            return 0;

        int chr_id = read->core.tid;
        int read_pos = read->core.pos;
        string chrom(header->target_name[chr_id]);

        uint32_t *cigar = bam_get_cigar(read);

        int off = 0;
        for (int i = 0; i < n_cigar; ++i) {
            char op = bam_cigar_op(cigar[i]);
            int ol = bam_cigar_oplen(cigar[i]);
//            cout<< "##" <<op<< read_pos <<"\n";
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF){
                 off += ol;
            }
            else if( op == BAM_CREF_SKIP || off >= MIN_BP_OVERLAP){
//                cout<< "@@@" << read_pos << "\n";
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
            const Junction * const jnc = junc_vec[jidx] ;
            const int npos = jnc->nreads_.size() ;
            default_random_engine generator;
            uniform_int_distribution<int> distribution(0,npos-1);
            for (int m=0; m<msamples; m++){
                float lambda = 0;
                for (int k=0; k<ksamples; k++) lambda += jnc->nreads_[distribution(generator)] ;
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

    const vector<Junction *>& IOBam::get_junc_vec(){
        return junc_vec ;
    }

}