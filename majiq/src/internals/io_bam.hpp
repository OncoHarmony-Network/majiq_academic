#ifndef IO_BAM_H
#define IO_BAM_H

//#include <algorithm>
#include <iomanip>
#include <iostream>
#include "grimoire.hpp"
#include "htslib/sam.h"
#include <map>
#include <set>
#include <omp.h>
//#include "interval.hpp"

#define MIN_INTRON_LEN 1000
#define MIN_BP_OVERLAP 8
#define MIN_JUNC_LENGTH 2
#define NUM_INTRON_BINS 10

#define UNSTRANDED 0
#define FWD_STRANDED 1
#define REV_STRANDED 2

#define is_read1(b)   (((b)->core.flag & BAM_FREAD1) != 0)
#define is_read2(b)   (((b)->core.flag & BAM_FREAD2) != 0)
#define is_read_reverse(b) (((b)->core.flag & BAM_FREVERSE) == 0x10)
////Sort a vector of junctions
//template <class CollectionType>
//inline void sort_junctions(CollectionType &junctions) {
//    sort(junctions.begin(), junctions.end(), compare_junctions);
//}

//The class that deals with creating the junctions
using namespace std;
using namespace grimoire;
//using namespace interval;

namespace io_bam{
    class IOBam {
        private:
            string bam_;
            int strandness_;
            unsigned int eff_len_;
            map<string, unsigned int> junc_map ;
            unsigned int nthreads_;
            map<string, vector<overGene*>> glist_ ;
            map<string, vector<Intron*>> intronVec_ ;
            unsigned int junc_limit_index_ ;

        public:
            vector<float *> junc_vec ;
            IOBam(){
                bam_ = string(".");
                strandness_ = 0 ;
                nthreads_ = 1 ;
            }

            IOBam(string bam1, int strandness1, unsigned int eff_len1, unsigned int nthreads1,
                  map<string, vector<overGene*>> glist1): strandness_(strandness1), eff_len_(eff_len1),
                                                                                nthreads_(nthreads1), glist_(glist1){
                bam_ = bam1 ;
            }

            IOBam(string bam1, int strandness1, unsigned int eff_len1): strandness_(strandness1), eff_len_(eff_len1){
                bam_ = string(bam1) ;
                nthreads_ = 1 ;
            }

            ~IOBam(){

                for(const auto &p1: junc_vec){
                    free(p1) ;
                }
                junc_vec.clear() ;
            }

            int parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) ;
            void add_junction(string chrom, char strand, int start, int end, int read_pos, int first_offpos, int sreads) ;
            int* get_junc_vec_summary() ;
            unsigned int get_junc_limit_index() { return junc_limit_index_ ; };
            int normalize_stacks(vector<float> vec, float sreads, int npos, float fitfunc_r, float pvalue_limit) ;
            int boostrap_samples(int msamples, int ksamples, float* boots, float fitfunc_r, float pvalue_limit) ;
            void detect_introns(float min_intron_cov, unsigned int min_experiments, float min_bins, bool reset) ;


            char _get_strand(bam1_t * read) ;
            void set_junction_strand(bam1_t  *aln, Junction& j1) ;
            void find_junction_genes(string chrom, char strand, int start, int end, float * nreads_ptr ) ;
            int  ParseJunctionsFromFile(bool ir_func) ;
            void parseJuncEntry(map<string, vector<overGene*>> & glist, string chrom, char strand, int start, int end,
                                int sreads, vector<Gene*>& oGeneList, bool ir) ;
            inline void update_junction_read(string key, int read_start, int count) ;
            inline bool new_junc_values(const string key) ;
            int parse_read_for_ir(bam_hdr_t *header, bam1_t *read) ;
            int get_njuncs() ;
            const map<string, unsigned int> &get_junc_map() ;
            const vector<Junction *>& get_junc_vec() ;
            void free_iobam() {
                int idx = 0 ;
                for(const auto &p1: junc_vec){
                    free(p1) ;
                    idx ++ ;
                }
                junc_vec.clear() ;
                junc_map.clear() ;
                intronVec_.clear() ;


            }



    };

    void prepare_genelist(map<string, Gene*>& gene_map, map<string, vector<overGene*>> & geneList) ;
    bool juncGeneSearch(Gene* t1, Junction* t2) ;
    void free_genelist(map<string, vector<overGene*>> & geneList) ;
}
#endif