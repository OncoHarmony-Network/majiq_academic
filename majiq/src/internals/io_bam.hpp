#ifndef IO_BAM_H
#define IO_BAM_H

//#include <algorithm>
#include <iomanip>
#include <iostream>
#include "junction.hpp"
#include "htslib/sam.h"
#include <map>
#include <set>

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
namespace io_bam{
    class IOBam {
        private:
            string bam_;
            int strandness_;
//            map<string, Junction> junc_dict_;
            string region_;
            set<string> set_chroms_ ;
            set<string> set_prejuncs_ ;

        public:
            map<string, Junction> junc_dict_;
            IOBam(){
                bam_ = string(".");
                strandness_ = 0;
                region_ = string(".");

            }

            IOBam(char* bam1, int strandness1): strandness_(strandness1){
                bam_ = string(bam1);
                region_ = string(".");
            }
            //Default constructor
            IOBam(char* bam1, int strandness1, char* region1) : strandness_(strandness1){
                bam_ = string(bam1);
                region_ = string(region1);
            }

            void create_junctions_vector();
            int find_junctions();
            int parse_read_into_junctions(bam_hdr_t *header, bam1_t *read);
            void add_junction(string chrom, char strand, int start, int end, int read_pos);
            char _get_strand(bam1_t * read);
            map<string, Junction> get_dict();
            void set_filters(set<string> set_chroms1, set<string> set_prejuncs1);
            void set_junction_strand(bam1_t *aln, Junction& j1);

    };
}
#endif