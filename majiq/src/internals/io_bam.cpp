
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "io_bam.hpp"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

using namespace std;

namespace io_bam {
    char IOBam::_get_strand(bam1_t * read){
        char strn = '.';

        if (strandness_ == FWD_STRANDED){
            strn = (((read->core.flag & 0x10) == (0x00 && is_read1(read)))
                    || ((read->core.flag & 0x10) == (0x10 && is_read2(read)))) ? '+' : '-';

        } else if (strandness_ == FWD_STRANDED){
            strn = (((read->core.flag & 0x10) ==(0x10 && is_read1(read)))
                    || ((read->core.flag & 0x10) == (0x00 && is_read2(read)))) ? '+' : '-';
        }
        return (strn);
    }

    inline int _unique(bam1_t * read){
        return ((read->core.flag & 0x100) != 0x100);
    }

    void IOBam::add_junction(string chrom, char strand, int start, int end, int read_pos) {

        stringstream s1;
        stringstream s2;

        s1<< start;
        s2<< end;
        string key = chrom +  ":" + strand + string(":") + s1.str() + "-" + s2.str();

        if (set_prejuncs_.count(key) >0){
            return;
        }
        else if(junc_dict_.count(key) == 0) {
            Junction j1 = Junction(chrom, strand, start, end);
            j1.nreads = 1;
            junc_dict_[key] = j1;
        }
        else {
            Junction j0 = junc_dict_[key];
            j0.nreads +=  1;
        }
        return;
    }


    int IOBam::parse_read_into_junctions(bam_hdr_t *header, bam1_t *read) {
        int n_cigar = read->core.n_cigar;
        if (n_cigar <= 1 || !_unique(read)) // max one cigar operation exists(likely all matches)
            return 0;

        int chr_id = read->core.tid;
        int read_pos = read->core.pos;
        string chrom(header->target_name[chr_id]);
        if (set_chroms_.count(chrom)==0){
            return 0;
        }
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
                    cout << e.what() << '\n';
                }

                off += ol;

            }
        }
        return 0;
    }

    //The workhorse - identifies junctions from BAM
    int IOBam::find_junctions() {
        cout << "FILE:" << bam_ << "\n";
        if(!bam_.empty()) {

            samFile *in = sam_open(bam_.c_str(), "r");
            if(in == NULL) {
                string msg = "[ERROR]: Unable to open file";
                throw runtime_error(msg);
            }

            hts_idx_t *idx = sam_index_load(in, bam_.c_str());
            if(idx == NULL) {
                string msg = "[ERROR]: Unable to open index file from ";
                throw runtime_error(msg);
            }

            bam_hdr_t *header = sam_hdr_read(in);
            hts_itr_t *iter = NULL;
            iter  = sam_itr_querys(idx, header, region_.c_str());

            if(header == NULL || iter == NULL) {
                sam_close(in);
                string msg = "[ERROR]: INVALID Region ";// << region_ << "in " << bam_ << "not found.";
                throw runtime_error(msg);
            }
            bam1_t *aln = bam_init1();
            while(sam_itr_next(in, iter, aln) >= 0) {
                parse_read_into_junctions(header, aln);
            }
            hts_itr_destroy(iter);
            hts_idx_destroy(idx);
            bam_destroy1(aln);
            bam_hdr_destroy(header);
            sam_close(in);

//            for (const auto &p : junc_dict_) {
//                cout << "m[" << p.first << "] = \n";
//}

        }
        return 0;
    }

    map<string, Junction> IOBam::get_dict(){
        return (junc_dict_);

    }

    void IOBam::set_filters(set<string> set_chroms1, set<string> set_prejuncs1){
        set_chroms_ = set_chroms1;
        set_prejuncs_ = set_prejuncs1;
    }

//    //Create the junctions vector from the map
//    void IOBam::create_junctions_vector() {
//        for(map<string, Junction> :: iterator it = junctions_.begin(); it != junctions_.end(); it++) {
//            Junction j1 = it->second;
//            junctions_vector_.push_back(j1);
//        }
//    }
}