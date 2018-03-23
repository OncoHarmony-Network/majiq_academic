
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
#include "htslib/thread_pool.h"

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

    void IOBam::add_junction(Gene* gobj, char strand, int start, int end, int read_pos) {

        const unsigned int offset = start - (read_pos+ MIN_BP_OVERLAP) ;
        if (offset >= eff_len_) {
            return ;
        }
        string key =  to_string(start) + "-" + to_string(end) ;
        if((gobj->junc_map).count(key) == 0) {
            (gobj->junc_map)[key]= new Junction(start, end, nexps_, eff_len_) ;
        }
        (gobj->junc_map)[key]->update_junction_read(read_pos, exp_index,  1) ;

        return;
    }


    int IOBam::parse_read_into_junctions(Gene* gobj, bam_hdr_t *header, bam1_t *read) {
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
                    add_junction(gobj, _get_strand(read), j_start, j_end, read_pos);
                } catch (const std::logic_error& e) {
                    cout << "KKK" << e.what() << '\n';
                }

                off += ol;

            }
        }
        return 0;
    }


    int IOBam::find_junctions_from_region(vector<Gene *> glist) {
        int exit_code = 0 ;
//        cout << "FILE:" << bam_ << "\n";
        if(!bam_.empty()) {
            int nthreads = 16 ;
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
            hts_itr_t *iter ;
            bam1_t *aln = bam_init1();

            htsThreadPool p = {NULL, 0};
            if (nthreads > 0) {
                p.pool = hts_tpool_init(nthreads);
                if (!p.pool) {
                    fprintf(stderr, "Error creating thread pool\n");
                    exit_code = 1;
                } else {
                    hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
                }
            }

            for(const auto &gg:glist){
                iter = sam_itr_querys(idx, header, gg->get_region().c_str());
                if(header == NULL || iter == NULL) {
//                    sam_close(in);
//                    bam_destroy1(aln);
                    string msg = "[ERROR]: INVALID Region ";// << region_ << "in " << bam_ << "not found.";
                    cout << msg << "\n" ;
                    continue ;
                }

                while(sam_itr_next(in, iter, aln) >= 0) {
                    parse_read_into_junctions(gg, header, aln);
                }
            }
            hts_itr_destroy(iter);
            bam_destroy1(aln);
            hts_idx_destroy(idx);
            bam_hdr_destroy(header);
            sam_close(in);

            if (p.pool)
                hts_tpool_destroy(p.pool);

//            for (const auto &p : denovo_juncs_) {
//                cout << "m[" << p.first << "] = \n";
//            }

        }
        return 0;
    }


//    int IOBam::parse_read_into_junctions2(Gene* gobj, bam_hdr_t *header, bam1_t *read) {
//        int n_cigar = read->core.n_cigar;
//        if (n_cigar <= 1 || !_unique(read)) // max one cigar operation exists(likely all matches)
//            return 0;
//
//        int chr_id = read->core.tid;
//        int read_pos = read->core.pos;
//        string chrom(header->target_name[chr_id]);
//
//        uint32_t *cigar = bam_get_cigar(read);
//
//        int off = 0;
//        for (int i = 0; i < n_cigar; ++i) {
//            char op = bam_cigar_op(cigar[i]);
//            int ol = bam_cigar_oplen(cigar[i]);
////            cout<< "##" <<op<< read_pos <<"\n";
//            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF){
//                 off += ol;
//            }
//            else if( op == BAM_CREF_SKIP || off >= MIN_BP_OVERLAP){
////                cout<< "@@@" << read_pos << "\n";
//                const int j_end = read->core.pos + off + ol +1;
//                const int j_start =  read->core.pos + off;
//                try {
//                    add_junction(gobj, _get_strand(read), j_start, j_end, read_pos);
//                } catch (const std::logic_error& e) {
//                    cout << "KKK" << e.what() << '\n';
//                }
//
//                off += ol;
//
//            }
//        }
//        return 0;
//    }


    int IOBam::ParseJunctionsFromFile(string filename, int nthreads){

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

        in = sam_open(filename.c_str(), "rb") ;
        if (NULL == in) {
            fprintf(stderr, "Error opening \"%s\"\n", filename.c_str());
            return EXIT_FAILURE;
        }
        header = sam_hdr_read(in);
        if (NULL == header) {
            fprintf(stderr, "Couldn't read header for \"%s\"\n", filename.c_str());
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
        if (nthreads > 0) {
            p.pool = hts_tpool_init(nthreads);
            if (!p.pool) {
                fprintf(stderr, "Error creating thread pool\n");
                exit_code = 1;
            } else {
                hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
            }
        }


    //    if (optind + 1 < argc && !(flag & READ_COMPRESSED)) { // BAM input and has a region
    //        int i;
    //        hts_idx_t *idx;
    //        if ((idx = sam_index_load(in, argv[optind])) == 0) {
    //            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
    //            return 1;
    //        }
    //        if (multi_reg) {
    //            int reg_count = 0;
    //            hts_reglist_t *reg_list = calloc(argc-(optind+1), sizeof(*reg_list));
    //            if (!reg_list)
    //                return 1;
    //
    //            // We need a public function somewhere to turn an array of region strings
    //            // into a region list, but for testing this will suffice for now.
    //            // Consider moving a derivation of this into htslib proper sometime.
    //            for (i = optind + 1; i < argc; ++i) {
    //                int j;
    //                uint32_t beg, end;
    //                char *cp = strrchr(argv[i], ':');
    //                if (cp) *cp = 0;
    //
    //                for (j = 0; j < reg_count; j++)
    //                    if (strcmp(reg_list[j].reg, argv[i]) == 0)
    //                        break;
    //                if (j == reg_count) {
    //                    reg_list[reg_count++].reg = argv[i];
    //                    if (strcmp(".", argv[i]) == 0) {
    //                        reg_list[j].tid = HTS_IDX_START;
    //
    //                    } else if (strcmp("*", argv[i]) == 0) {
    //                        reg_list[j].tid = HTS_IDX_NOCOOR;
    //
    //                    } else {
    //                        int k; // need the header API here!
    //                        for (k = 0; k < h->n_targets; k++)
    //                            if (strcmp(h->target_name[k], argv[i]) == 0)
    //                                break;
    //                        if (k == h->n_targets)
    //                            return 1;
    //                        reg_list[j].tid = k;
    //                        reg_list[j].min_beg = h->target_len[k];
    //                        reg_list[j].max_end = 0;
    //                    }
    //                }
    //
    //                hts_reglist_t *r = &reg_list[j];
    //                r->intervals = realloc(r->intervals, ++r->count * sizeof(*r->intervals));
    //                if (!r->intervals)
    //                    return 1;
    //                beg = 1;
    //                end = r->tid >= 0 ? h->target_len[r->tid] : 0;
    //                if (cp) {
    //                    *cp = 0;
    //                    // hts_parse_reg() is better, but awkward here
    //                    sscanf(cp+1, "%d-%d", &beg, &end);
    //                }
    //                r->intervals[r->count-1].beg = beg-1; // BED syntax
    //                r->intervals[r->count-1].end = end;
    //
    //                if (r->min_beg > beg)
    //                    r->min_beg = beg;
    //                if (r->max_end < end)
    //                    r->max_end = end;
    //            }
    //
    //            hts_itr_multi_t *iter = sam_itr_regions(idx, h, reg_list, reg_count);
    //            if (!iter)
    //                return 1;
    //            while ((r = sam_itr_multi_next(in, iter, b)) >= 0) {
    //                if (!benchmark && sam_write1(out, h, b) < 0) {
    //                    fprintf(stderr, "Error writing output.\n");
    //                    exit_code = 1;
    //                    break;
    //                }
    //                if (nreads && --nreads == 0)
    //                    break;
    //            }
    //            hts_itr_multi_destroy(iter);
    //        } else {
    //            for (i = optind + 1; i < argc; ++i) {
    //                hts_itr_t *iter;
    //                if ((iter = sam_itr_querys(idx, h, argv[i])) == 0) {
    //                    fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
    //                    continue;
    //                }
    //                while ((r = sam_itr_next(in, iter, b)) >= 0) {
    //                    if (!benchmark && sam_write1(out, h, b) < 0) {
    //                        fprintf(stderr, "Error writing output.\n");
    //                        exit_code = 1;
    //                        break;
    //                    }
    //                    if (nreads && --nreads == 0)
    //                        break;
    //                }
    //                hts_itr_destroy(iter);
    //            }
    //        }
    //        hts_idx_destroy(idx);
    //    } else
        while ((r = sam_read1(in, header, aln)) >= 0) {
            cout<<"read" << aln->core.tid << " :: " << aln->core.pos << "\n" ;

            if (nreads && --nreads == 0)
                break;
        }

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



}