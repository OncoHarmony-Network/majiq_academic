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

            vector<vector<unsigned int>*> junc_vec ;
            map<string, unsigned int> pos_map_ ;

            unsigned int nthreads_;
            map<string, vector<Gene*>> glist_ ;

        public:
            IOBam(){
                bam_ = string(".");
                strandness_ = 0 ;
                nthreads_ = 1 ;
            }

            IOBam(string bam1, int strandness1, unsigned int eff_len1, unsigned int nthreads1,
                            map<string, vector<Gene*>> &glist1): strandness_(strandness1), eff_len_(eff_len1),
                                                                                nthreads_(nthreads1), glist_(glist1){
                bam_ = bam1;
            }
            //Default constructor
            IOBam(string bam1, int strandness1, unsigned int eff_len1): strandness_(strandness1), eff_len_(eff_len1){
                bam_ = string(bam1);
                nthreads_ = 1 ;
            }
            ~IOBam(){
//                delete junc_map ;
//                delete junc_vec ;
//                delete pos_map_ ;
            }
//                for(const auto &p1: junc_vec){
//
//                    delete p1 ;
//                }
//                junc_map.clear() ;
////                junc_map.shrink_to_fit() ;
//
//            }


            int parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) ;
            void add_junction(string chrom, char strand, int start, int end, int read_pos) ;

            char _get_strand(bam1_t * read) ;
            void set_junction_strand(bam1_t  *aln, Junction& j1) ;
//            void find_junction_gene(string chrom, char strand, Junction * junc) ;
//            void find_junction_gene(string chrom, char strand, const Junction * junc) ;
            void find_junction_genes(string chrom, char strand, int start, int end,
                                    vector<unsigned int> * nreads_ptr ) ;
            int ParseJunctionsFromFile() ;
            inline void update_junction_read(string key, int read_start, int count) ;

            int boostrap_samples(int msamples, int ksamples, float* boots) ;
//             int boostrap_samples(int msamples, int ksamples, float* boots, char** out_id) ;
            int get_njuncs() ;
            const map<string, unsigned int> &get_junc_map() ;
            const vector<Junction *>& get_junc_vec() ;
            int* get_junc_vec_summary() ;
    };

    bool juncGeneSearch(Gene* t1, Junction* t2) ;
}
#endif