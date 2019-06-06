#ifndef GRIMOIRE_H
#define GRIMOIRE_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <string>
#include <omp.h>
#include "psi.hpp"


#define EMPTY_COORD  -1
#define FIRST_LAST_JUNC  -2
#define MAX_DENOVO_DIFFERENCE 400
#define MIN_INTRON_BINSIZE 1000
#define MAX_TYPE_LENGTH 245
#define NA_LSV  "na"
using namespace std ;
//extern typedef vector<float> psi_distr_t ;
typedef vector<float> psi_distr_t ;

namespace grimoire{
    class Exon ;
    class Junction ;
    class Gene ;
    class Intron ;
    class LSV ;

    struct Ssite {
        int         coord ;
        bool        donor_ss ;
        Junction *  j ;
    };

    class Jinfo {
        public:
            unsigned int    index ;
            int             sreads ;
            int             npos ;

            Jinfo() {}
            Jinfo(unsigned int index1, int sreads1, int npos1): index(index1), sreads(sreads1), npos(npos1) {}
    } ;

    struct lsvtype {
        int         coord ;
        int         ref_coord ;
        Exon *      ex_ptr ;
        Junction *  jun_ptr ;
    } ;

    class _Region{
        protected:
            int start_ ;
            int end_ ;

            _Region(){}
            _Region(int start1, int end1): start_(start1), end_(end1) {}
            virtual ~_Region(){}
        public:
            int         get_start()             { return start_ ; }
            int         get_end()               { return end_ ; }
            void        set_start(int start1)   { start_ = start1 ; }
            void        set_end(int end1)       { end_ = end1 ; }
            inline int  length()                { return end_ - start_ ; }
            virtual void  set_simpl_fltr(bool val, bool in) { cerr << "NOT SHOW\n" ;}


            static bool func_comp (_Region* a, int coord){
                return a->get_end() < coord ;
            }

            template <class myRegion>
            static bool islowerRegion(myRegion * a, myRegion * b){
                return (a->get_start() < b->get_start()) ||
                        (a->get_start() == b->get_start() && a->get_end() < b->get_end());
            }

            template <class myRegion>
            static bool isgreaterRegion(myRegion * a, myRegion * b){
                return (a->get_start() > b->get_start()) ||
                        (a->get_start() == b->get_start() && a->get_end() > b->get_end());
            }

            template <class myRegion>
            static bool RegionsOverlap(myRegion* t1, myRegion* t2){
                return ( (t1->get_start() <= t2->get_end()) && (t1->get_end() >= t2->get_start()) ) ;
            }



    };


    class Junction: public _Region{
        private:
            bool annot_ ;
            bool intronic_ ;
            bool bld_fltr_ ;
            bool denovo_bl_ ;
            bool simpl_fltr_ ;
            unsigned int denovo_cnt_ ;
            unsigned int flter_cnt_ ;
            unsigned int simpl_cnt_in_ ;
            unsigned int simpl_cnt_out_ ;
            Exon * acceptor_;
            Exon * donor_;

        public:
            float* nreads_ ;

            Junction() {}
            Junction(int start1, int end1, bool annot1, bool simpl1): _Region(start1, end1),
                                                                      annot_(annot1), simpl_fltr_(simpl1){
                denovo_bl_ = annot1 ;
//                denovo_bl_ = false ;
                denovo_cnt_ = 0 ;
                bld_fltr_ = false ;
                flter_cnt_ = 0 ;
                intronic_ = false ;
                nreads_ = nullptr ;
                simpl_cnt_in_ = 0 ;
                simpl_cnt_out_ = 0 ;
            }
            ~Junction()             { clear_nreads(true) ; }

            string  get_key()       { return(to_string(start_) + "-" + to_string(end_)) ; }
            string  get_key(Gene * gObj) ;
            string  get_key(Gene * gObj, int strandness) ;
            bool    get_annot()     { return annot_ ; }
            bool    get_intronic()  { return intronic_ ; }
            bool    get_bld_fltr()  { return bld_fltr_ ; }
            bool    get_simpl_fltr(){ return simpl_fltr_ ; }
            bool    get_denovo_bl() { return denovo_bl_ ; }
            Exon*   get_acceptor()  { return acceptor_ ; }
            Exon*   get_donor()     { return donor_ ; }
            void  set_nreads_ptr(float * nreads1) { nreads_ = nreads1 ; }
            void  set_acceptor(Exon * acc) { acceptor_ = acc ; }
            void  set_donor(Exon * don) { donor_ = don ; }
            void  exonReset(){
                acceptor_ = nullptr ;
                donor_ = nullptr ;
            }

            void set_simpl_fltr(bool val, bool in){
                if (in)
                    simpl_cnt_in_ += val? 1 : 0 ;
                else
                    simpl_cnt_out_ += val? 1 : 0 ;
            }

            void update_simpl_flags(unsigned int min_experiments){
                simpl_fltr_ = simpl_fltr_ && simpl_cnt_in_ >= min_experiments && simpl_cnt_out_ >= min_experiments ;
                simpl_cnt_in_ = 0 ;
                simpl_cnt_out_ = 0 ;
            }

            inline void update_flags(unsigned int sreads, unsigned int minreads_t, unsigned int npos, unsigned int minpos_t,
                              unsigned int denovo_t, unsigned int min_experiments, bool denovo){
                if (( npos >= minpos_t) && (sreads >= minreads_t)){
                    ++ flter_cnt_ ;
                    bld_fltr_ = bld_fltr_ || (flter_cnt_ >= min_experiments) ;
                }
                if (sreads >= denovo_t){
                    ++ denovo_cnt_  ;
                    denovo_bl_ = denovo_bl_ || (denovo_cnt_ >= min_experiments) ;
                    if (!(denovo || annot_))
                        denovo_bl_ = false ;
                }

                return ;
            }

            void gen_and_update_flags(int efflen, unsigned int num_reads, unsigned int num_pos, unsigned int denovo_thresh,
                              unsigned int min_experiments, bool denovo){
                if (nreads_ == nullptr ){
                    return ;
                }
                unsigned int sum_reads = 0 ;
                unsigned int npos = 0 ;

                for(int i=0; i<efflen; ++i){
                    sum_reads += nreads_[i] ;
                    npos += nreads_[i]? 1 : 0 ;
                }
                update_flags(sum_reads, num_reads, npos, num_pos, denovo_thresh, min_experiments,  denovo) ;
//cerr << "UPDATE FLAGS " << get_key() << " bool:" << denovo_bl_<< " denovothresh: " << denovo_thresh << " sum_reads: " << sum_reads << "\n" ;
                return ;
            }


            void clear_nreads(bool reset_grp){
                nreads_ = nullptr ;
                flter_cnt_ = reset_grp ? 0: flter_cnt_ ;
                denovo_cnt_ = reset_grp ? 0: denovo_cnt_ ;
            }

            void update_junction_read(int read_start, unsigned int n) ;

    };
    class Exon: public _Region{

        public:
            bool            annot_ ;
            int             db_start_ ;
            int             db_end_ ;
            set<Junction *> ib ;
            set<Junction *> ob ;
            Intron *        ib_irptr ;
            Intron *        ob_irptr ;


            Exon(){}
            Exon(int start1, int end1): _Region(start1, end1), db_start_(start1), db_end_(end1){
                annot_  = true ;
                ib_irptr = nullptr ;
                ob_irptr = nullptr ;
            }
            Exon(int start1, int end1, bool annot1): _Region(start1, end1), annot_(annot1),
                                                                            db_start_(start1), db_end_(end1){
                ib_irptr = nullptr ;
                ob_irptr = nullptr ;
            }
            ~Exon(){}
            bool    is_lsv(bool ss) ;
            void    simplify(map<string, int>& junc_tlb, float simpl_percent, Gene* gObj, int strandness,
                        int denovo_simpl, int db_simple, int ir_simpl, bool last, unsigned int min_experiments) ;
            bool    has_out_intron()        { return ob_irptr != nullptr ; }
            string  get_key()       { return(to_string(start_) + "-" + to_string(end_)) ; }
            void set_simpl_fltr(bool val) {};
            void    revert_to_db(){
                set_start(db_start_) ;
                set_end(db_end_) ;
                for (auto const &j: ib){
                    j->exonReset() ;
                }
                ib.clear() ;

                for (auto const &j: ob){
                    j->exonReset() ;
                }
                ob.clear() ;
            }

    };

    class Intron: public _Region{
        private:
            bool    annot_ ;
            Gene *  gObj_ ;
            int     flt_count_ ;
            bool    ir_flag_ ;
            bool    markd_ ;

            int     nxbin_off_ ;
            int     nxbin_mod_ ;
            int     nxbin_ ;
            int     numbins_ ;
            bool    simpl_fltr_ ;
            unsigned int simpl_cnt_in_ ;
            unsigned int simpl_cnt_out_ ;

        public:
            float*  read_rates_ ;

            Intron () {}
            Intron (int start1, int end1, bool annot1, Gene* gObj1, bool simpl1): _Region(start1, end1),
                                                                     annot_(annot1), gObj_(gObj1), simpl_fltr_(simpl1){
                flt_count_  = 0 ;
                ir_flag_    = false ;
                read_rates_ = nullptr ;
                nxbin_      = 0 ;
                nxbin_off_  = 0 ;
                nxbin_mod_  = 0 ;
                markd_      = false ;
                numbins_    = 0 ;
                simpl_cnt_in_  = 0 ;
                simpl_cnt_out_  = 0 ;
            }

            bool    get_annot()             { return annot_ ; }
            Gene*   get_gene()              { return gObj_ ; }
            bool    get_ir_flag()           { return ir_flag_ ; }
            string  get_key()               { return (to_string(start_) + "-" + to_string(end_)) ; }
            bool    get_simpl_fltr()        { return simpl_fltr_ ; }
            int     get_nxbin()             { return nxbin_ ; }
            int     get_nxbin_off()         { return nxbin_off_ ; }
            int     get_nxbin_mod()         { return nxbin_mod_ ; }
            int     get_numbins()           { return numbins_ ; }
            int     get_fltcount()          { return flt_count_ ; }
            void    set_markd()             { markd_ = true ; }
            void    unset_markd()           { markd_ = false ; }
            bool    is_connected()          { return markd_ ; }
            string  get_key(Gene * gObj) ;
            void    calculate_lambda() ;
            bool    is_reliable(float min_bins, int eff_len) ;

            void set_simpl_fltr(bool val, bool in){
                if (in)
                    simpl_cnt_in_ += val? 1 : 0 ;
                else
                    simpl_cnt_out_ += val? 1 : 0 ;
            }

            void update_simpl_flags(unsigned int min_experiments){

                simpl_fltr_ = simpl_fltr_ && simpl_cnt_in_ >= min_experiments && simpl_cnt_out_ >= min_experiments ;
                simpl_cnt_in_  = 0 ;
                simpl_cnt_out_  = 0 ;

            }

            void add_read_rates_buff(int eff_len){
                nxbin_      = (int) ((length()+ eff_len) / eff_len) ;
                nxbin_mod_  = (length() % eff_len) ;
                nxbin_off_  = nxbin_mod_ * ( nxbin_+1 ) ;
                numbins_    = eff_len ;
                read_rates_ = (float*) calloc(numbins_, sizeof(float)) ;
            }

            void initReadCovFromVector(vector<float>& cov){
                for (unsigned int i=0; i<cov.size(); ++i) {
                    read_rates_[i] = cov[i] ;
                }

            }

            void  add_read(int read_pos, int eff_len, int s){
                int st = get_start() - eff_len ;
                if (read_rates_ == nullptr){
                    add_read_rates_buff(eff_len) ;
                }
                int offset = (read_pos - st) ;
                offset = (offset >=0) ? offset: 0 ;

                const int off1 = (int) ((offset - nxbin_off_) / nxbin_) + nxbin_mod_ ;
                const int off2 = (int) offset / (nxbin_+1) ;
                offset = (int)(offset < nxbin_off_) ? off2: off1 ;
                read_rates_[offset] += s ;
            }

            inline void  update_flags(float min_coverage, int min_exps, float min_bins) {
                int cnt = 0 ;
                const float pc_bins = numbins_ * min_bins ;
                for(int i =0 ; i< numbins_; i++){
                    cnt += (read_rates_[i]>= min_coverage) ? 1 : 0 ;
                }
                flt_count_ += (cnt >= pc_bins) ? 1 : 0 ;
                ir_flag_ = ir_flag_ || (flt_count_ >= min_exps) ;
                return ;
            }

            inline void update_boundaries(int start, int end){
                int st = get_start() ;
                int nd = get_end() ;

                st = (st > start)? st: start ;
                nd = (nd < end)? nd: end ;

                set_end(nd) ;
                set_start(st) ;
            }

            void overlaping_intron(Intron * inIR_ptr){
                start_ = max(start_, inIR_ptr->get_start()) ;
                end_ = min(end_, inIR_ptr->get_end()) ;

                nxbin_      = inIR_ptr->get_nxbin() ;
                nxbin_mod_  = inIR_ptr->get_nxbin_mod() ;
                nxbin_off_  = inIR_ptr->get_nxbin_off() ;
                numbins_    = inIR_ptr->get_numbins() ;
                read_rates_ = inIR_ptr->read_rates_ ;
                flt_count_ += inIR_ptr->get_fltcount() ;
                ir_flag_    = ir_flag_ || inIR_ptr->get_ir_flag() ;
                return ;
            }

            void clear_nreads(bool reset_grp){
                flt_count_ = reset_grp ? 0: flt_count_ ;
                return ;
            }

            void free_nreads(){
                if(read_rates_ != nullptr){
                    free(read_rates_);
                    read_rates_ = nullptr ;
                }
                return ;
            }
    };


    class Gene: public _Region{
        private:
            string  id_ ;
            string  name_ ;
            string  chromosome_ ;
            char    strand_ ;
            omp_lock_t map_lck_ ;

        public:
            map <string, Junction*> junc_map_ ;
            map <string, Exon*> exon_map_ ;
            vector <Intron*> intron_vec_ ;


            Gene (){}
            Gene(string id1, string name1, string chromosome1,
                 char strand1, unsigned int start1, unsigned int end1): _Region(start1, end1), id_(id1), name_(name1),
                                                                        chromosome_(chromosome1), strand_(strand1){
                omp_init_lock( &map_lck_ ) ;
            }

            ~Gene(){
                for(const auto &p: exon_map_){
                    delete p.second ;
                }
                for(const auto &p: junc_map_){
                    delete p.second ;
                }
                for(const auto &p: intron_vec_){
                    delete p ;
                }
            }
            string  get_id()        { return id_ ; }
            string  get_chromosome(){ return chromosome_ ;}
            char    get_strand()    { return strand_ ;}
            string  get_name()      { return name_ ;}
            void    set_simpl_fltr(bool val) {};
            void    create_annot_intron(int start_ir, int end_ir, bool simpl){
                Intron * ir = new Intron(start_ir +1, end_ir -1, true, this, simpl) ;
                intron_vec_.push_back(ir) ;
            }

            void reset_exons(){
                for (auto p  = exon_map_.begin(); p!= exon_map_.end();){
                    map <string, Exon*>::iterator pit = p ;
                    Exon * e = p->second ;
                    e->revert_to_db() ;
                    ++p ;
                    if (!e->annot_){
                        exon_map_.erase(pit) ;
                        delete e ;
                    }
                }

            }

            string  get_region() ;
            void    print_exons() ;
            void    detect_exons() ;
            void    connect_introns() ;
            void    detect_introns(vector<Intron*> &intronlist, bool simpl) ;
            void    add_intron(Intron * inIR_ptr, float min_coverage, unsigned int min_exps, float min_bins, bool reset) ;
            void    newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db) ;
            void    fill_junc_tlb(map<string, vector<string>> &tlb) ;
            int     detect_lsvs(vector<LSV*> &out_lsvlist) ;
            void    simplify(map<string, int>& junc_tlb, float simpl_percent, int strandness, int denovo_simpl,
                             int db_simple, int ir_simpl, bool last, unsigned int min_experiments) ;
            void    initialize_junction(string key, int start, int end, float* nreads_ptr, bool simpl) ;
            void    update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                      unsigned int denovo_thresh, unsigned int min_experiments, bool denovo) ;
            void    updateFlagsFromJunc(string key, unsigned int sreads, unsigned int minreads_t, unsigned int npos,
                                        unsigned int minpos_t, unsigned int denovo_t, bool denovo, int minexp,
                                        bool reset) ;

    };

    class LSV{
        private:
            string              id_ ;
            vector<Junction *>  junctions_ ;
            Intron *            ir_ptr_ ;
            string              type_ ;
            Gene *              gObj_ ;

        public:
            LSV(){}
            LSV(Gene* gObj1, Exon* ex, bool ss): gObj_(gObj1){
                const bool b = (gObj_->get_strand() == '+') ;

                const string t = (ss != b) ? "t" : "s" ;
                const string st =(ex->get_start() == -1) ? NA_LSV : to_string(ex->get_start()) ;
                const string end = (ex->get_end() == -1) ? NA_LSV : to_string(ex->get_end()) ;
                id_     = gObj_->get_id() + ":" + t + ":" + st + "-" + end ;

                ir_ptr_ = nullptr ;
                Intron * temp_ir = ss? ex->ob_irptr : ex->ib_irptr ;
                if (temp_ir != nullptr )
                    ir_ptr_ = (temp_ir->get_ir_flag() && !temp_ir->get_simpl_fltr()) ? temp_ir : nullptr ;

                type_   = set_type(ex, ss) ;
            }

            ~LSV(){}
            bool gather_lsv_info (float* source, float* target, list<Jinfo*> &info, map<string, Jinfo> &tlb,
                                  unsigned int msample) ;
            string      set_type (Exon* ex, bool ss) ;
            inline void get_variations(set<string> &t1) ;
            string      get_type ()                  { return type_ ; }
            string      get_id ()                    { return id_ ; }
            Gene*       get_gene()                   { return gObj_ ; }
            int         get_num_junctions()          { return junctions_.size() ; }
            vector<Junction *>& get_junctions()      { return junctions_ ; }
            Intron *    get_intron()                 {
                Intron* irptr = (ir_ptr_ != nullptr) ? ir_ptr_ : (Intron * ) 0 ;
                return irptr ;

            }
            int  get_num_variations() {
                const int n = junctions_.size() ;
                const int c = (ir_ptr_ != nullptr) ? 1 : 0 ;
                return (n+c) ;
            }
    };

    class overGene: public _Region{
        public:
            vector<Gene *> glist ;
        overGene(unsigned int start1, unsigned int end1):_Region(start1, end1) { }
        overGene(){}

        ~overGene(){
            for (auto &g: glist){
                delete g ;
            }
            glist.clear() ;
            glist.shrink_to_fit() ;
        }
    } ;


    void sortGeneList(vector<Gene*> &glist) ;
    vector<Intron *> find_intron_retention(Gene * gObj, int start, int end);
//    vector<Gene*> find_gene_from_junc(const map<string, vector<overGene*>>& glist, string chrom, int start, int end, bool ir) ;

    void find_gene_from_junc(map<string, vector<overGene*>> & glist, string chrom, char strand, int start, int end,
                             vector<Gene*>& oGeneList, bool ir, bool simpl) ;
    void fill_junc_tlb(vector<LSV*>& lsv_list, map<string, int>& tlb) ;
    bool isNullJinfo(Jinfo* x) ;
    void free_JinfoVec(vector<Jinfo*>& jvec);
    string key_format(string gid, int coord1, int coord2, bool ir) ;
    void free_lsvlist(vector<LSV*> & lsvList) ;
}

#endif

