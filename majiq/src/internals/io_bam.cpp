
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

    void IOBam::add_junction(string gne_id, char strand, int start, int end, int read_pos) {

        const unsigned int offset = start - (read_pos+ MIN_BP_OVERLAP) ;
        if (offset >= eff_len_) {
            return ;
        }
        string key =  gne_id + string(":") + to_string(start) + "-" + to_string(end) ;
//        cout << "IN ADD JUNC " << nexps_<< " " << eff_len_ <<"\n"<< flush;
        if(denovo_juncs_.count(key) == 0) {
            denovo_juncs_[key]= new Junction(gne_id, start, end, nexps_, eff_len_) ;
        }
        denovo_juncs_[key]->update_junction_read(read_pos, exp_index,  1) ;
//        cout << "END ADD JUNC\n"<< flush;
        return;
    }


    int IOBam::parse_read_into_junctions(string gene_id, bam_hdr_t *header, bam1_t *read) {
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
                    add_junction(gene_id, _get_strand(read), j_start, j_end, read_pos);
                } catch (const std::logic_error& e) {
                    cout << "KKK" << e.what() << '\n';
                }

                off += ol;

            }
        }
        return 0;
    }


    int IOBam::find_junctions_from_region(Gene * gobj) {
//        cout << "FILE:" << bam_ << "\n";
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
            // for(int i; i<N ; i++){

                hts_itr_t *iter = sam_itr_querys(idx, header, gobj->get_region().c_str());
                if(header == NULL || iter == NULL) {
                    sam_close(in);
                    string msg = "[ERROR]: INVALID Region ";// << region_ << "in " << bam_ << "not found.";
                    return 0;
//                    throw runtime_error(msg);
                }
                bam1_t *aln = bam_init1();
                while(sam_itr_next(in, iter, aln) >= 0) {
                    parse_read_into_junctions(gobj->id, header, aln);
                }

            hts_itr_destroy(iter);
            bam_destroy1(aln);
            hts_idx_destroy(idx);
            bam_hdr_destroy(header);
            sam_close(in);

//            for (const auto &p : denovo_juncs_) {
//                cout << "m[" << p.first << "] = \n";
//            }

        }
        return 0;
    }



//
//    map<string, Junction> IOBam::get_dict(){
//        return (junc_dict_);
//
//    }
//
//    void IOBam::set_filters(map<string, set<string>>* prejuncs1){
//        set_prejuncs_ = prejuncs1;
//    }





}