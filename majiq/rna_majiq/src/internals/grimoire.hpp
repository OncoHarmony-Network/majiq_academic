#ifndef GRIMOIRE_H
#define GRIMOIRE_H

#include <iostream>
#include <map>
#include <vector>
#include <memory>
#include <set>
#include <list>
#include <string>
#include <omp.h>
#include "psi.hpp"


#define EMPTY_COORD  -1
#define FIRST_LAST_JUNC  -2
#define MAX_DENOVO_DIFFERENCE 400
#define MAX_TYPE_LENGTH 245
#define NA_LSV  "na"
using namespace std ;
typedef vector<float> psi_distr_t ;
typedef std::pair<int, int> coord_key_t;
typedef string lsv_id_t;

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
            int             numpos ;

            Jinfo() {}
            Jinfo(unsigned int index1, int sreads1, int npos1): index(index1), sreads(sreads1), numpos(npos1) {}
    } ;

    struct lsvtype {
        int         coord ;
        int         ref_coord ;
        Exon *      ex_ptr ;
        Junction *  jun_ptr ;
    } ;

    /**
     * mix-in class for objects that can be simplified
     */
    class _Simplifiable {
        protected:
            // is this object simplified?
            bool simpl_fltr_;

        private:
            // within current build group, number of experiments that give
            // evidence against simplification
            unsigned int simpl_cnt_in_;
            unsigned int simpl_cnt_out_;

        public:
            /**
             * virtual destructor
             */
            virtual ~_Simplifiable() {}
            /**
             * empty constructor
             */
            _Simplifiable() {}
            /**
             * Initialize simplifiable details
             *
             * @param simpl1 Does the object start simplified or not?
             */
            _Simplifiable(bool simpl1) : simpl_fltr_(simpl1),
                                         simpl_cnt_in_(0),
                                         simpl_cnt_out_(0) {
            }
            /**
             * Is the object simplified?
             */
            bool get_simpl_fltr() const { return simpl_fltr_; }
            /**
             * Count experiment as having evidence for or against
             * simplification in specified direction
             *
             * @param val experiment has evidence for keeping it simplified
             * @param in evidence is for junction as part of event going in vs
             * out of exon
             */
            void set_simpl_fltr(bool val, bool in) {
                if (!val) {
                    // we only count experiments where evidence is against
                    // simplification
                    if (in) {
                        ++simpl_cnt_in_;
                    } else {
                        ++simpl_cnt_out_;
                    }
                }
            }
            /**
             * Update if simplified using aggregated build group and reset counts
             *
             * @param min_experiments Minimum number of experiments with
             * evidence against simplification, making it unsimplified
             *
             * @note object stays simplified only if both directions had
             * insufficient experiments to make not simplified
             */
            void update_simpl_flags(unsigned int min_experiments) {
                // we only need to do anything if currently simplified
                if (simpl_fltr_) {
                    // object stays simplified only if both directions had
                    // insufficient experiments to make not simplified
                    simpl_fltr_ = simpl_cnt_in_ < min_experiments
                        && simpl_cnt_out_ < min_experiments;
                }
                // reset counts
                simpl_cnt_in_ = 0;
                simpl_cnt_out_ = 0;
            }
    };

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
            inline int  length()                { return 1 + end_ - start_ ; }


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


    class Junction: public _Region, public _Simplifiable {
        private:
            bool annot_ ;
            bool intronic_ ;
            bool bld_fltr_ ;
            bool denovo_bl_ ;
            bool constitutive_ ;
            unsigned int denovo_cnt_ ;
            unsigned int flter_cnt_ ;
            Exon * acceptor_;
            Exon * donor_;

        public:
            shared_ptr<vector<float>> nreads_ptr_;

            Junction() {}
            Junction(int start1, int end1, bool annot1, bool simpl1)
                    : _Region(start1, end1),
                      _Simplifiable(simpl1),
                      annot_(annot1) {
                denovo_bl_ = annot1 ;
                denovo_cnt_ = 0 ;
                bld_fltr_ = false ;
                flter_cnt_ = 0 ;
                intronic_ = false ;
                nreads_ptr_ = nullptr ;
                constitutive_  = false ;
            }
            ~Junction()             { clear_nreads(true) ; }

            coord_key_t get_key() { return std::make_pair(start_, end_); }
            string  get_key(Gene * gObj) ;
            string  get_key(Gene * gObj, int strandness) ;
            bool    get_annot()     { return annot_ ; }
            bool    get_intronic()  { return intronic_ ; }
            bool    get_bld_fltr()  { return bld_fltr_ ; }
            bool    get_denovo_bl() { return denovo_bl_ ; }
            Exon*   get_acceptor()  { return acceptor_ ; }
            Exon*   get_donor()     { return donor_ ; }
            bool    is_exitron()    { return get_donor() == get_acceptor() ; }

            bool  get_constitutive() { return constitutive_ ; }
            void  set_constitutive() { constitutive_ = true ; }
            void  set_nreads_ptr(shared_ptr<vector<float>> nreads1) { nreads_ptr_ = nreads1 ; }
            void  set_acceptor(Exon * acc) { acceptor_ = acc ; }
            void  set_donor(Exon * don) { donor_ = don ; }
            void  exonReset(){
                acceptor_ = nullptr ;
                donor_ = nullptr ;
            }

            /**
             * Identify if junction passes current experiments to update build/denovo filters
             *
             * @param sreads the number of reads for junction in current experiment
             * @param minreads_t the minimum number of reads to pass build filter
             * @param numpos the number of nonzero positions in current experiment
             * @param minpos_t the minimum number of nonzero positions to pass
             * build/denovo filters
             * @param denovo_t the minimum number of reads to pass denovo filter
             * @param min_experiments the total number of experiments in a
             * build group required to pass build/denovo filters
             * @param denovo are we accepting denovo junctions?
             *
             * @note Update counts of experiments passing build or denovo
             * filters, and update flag indicating if junction passed
             * corresponding filter if enough experiments passed in this build
             * group. Only requires one build group to pass
             */
            inline void update_flags(unsigned int sreads, unsigned int minreads_t, unsigned int numpos, unsigned int minpos_t,
                              unsigned int denovo_t, unsigned int min_experiments, bool denovo){
                // only update flags if experiment has enough nonzero positions
                if (numpos >= minpos_t) {
                    // only try updating build filter if hasn't passed and enough reads
                    if (!bld_fltr_ && sreads >= minreads_t) {
                        // increment number of experiments passing build filters
                        ++flter_cnt_;
                        // update flag for junction if enough experiments passed
                        if (flter_cnt_ >= min_experiments) {
                            bld_fltr_ = true;
                        }
                    }
                    // only try updating denovo filter if processing denovo,
                    // hasn't passed, and enough reads
                    if (denovo && !denovo_bl_ && sreads >= denovo_t) {
                        // increment number of experiments passing build filters
                        ++denovo_cnt_;
                        // update flag for junction if enough experiments passed
                        if (denovo_cnt_ >= min_experiments) {
                            denovo_bl_ = true;
                        }
                    }
                }
                return ;
            }

            void gen_and_update_flags(int efflen, unsigned int num_reads, unsigned int num_pos, unsigned int denovo_thresh,
                              unsigned int min_experiments, bool denovo){
                if (nreads_ptr_ == nullptr) {
                    return ;
                }
                const vector<float> nreads = *nreads_ptr_;
                unsigned int sum_reads = 0 ;
                unsigned int numpos = 0 ;

                for(int i=0; i<efflen; ++i){
                    sum_reads += nreads[i] ;
                    numpos += nreads[i]? 1 : 0 ;
                }
                update_flags(sum_reads, num_reads, numpos, num_pos, denovo_thresh, min_experiments,  denovo) ;
                return ;
            }


            void clear_nreads(bool reset_grp){
                nreads_ptr_ = nullptr ;
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
            coord_key_t get_key() { return std::make_pair(start_, end_); }
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

    class Intron: public _Region, public _Simplifiable {
        private:

            bool    annot_ ;
            Gene *  gObj_ ;
            bool    ir_flag_ ;

            int     flt_count_ ;
            bool    markd_ ;
            bool    const_donor_ ;
            bool    const_acceptor_ ;
            int     nxbin_off_ ;
            int     nxbin_mod_ ;
            int     nxbin_ ;
            int     numbins_ ;

        public:
            shared_ptr<vector<float>> read_rates_ptr_;

            Intron () {}
            Intron (int start1, int end1, bool annot1, Gene* gObj1, bool simpl1)
                    : _Region(start1, end1),
                      _Simplifiable(simpl1),
                      annot_(annot1),
                      gObj_(gObj1) {
                flt_count_  = 0 ;
                ir_flag_    = false ;
                read_rates_ptr_ = nullptr;
                nxbin_      = 0 ;
                nxbin_off_  = 0 ;
                nxbin_mod_  = 0 ;
                markd_      = false ;
                numbins_    = 0 ;
                const_donor_    = false ;
                const_acceptor_ = false ;
            }

             Intron (int start1, int end1, bool annot1, Gene* gObj1, bool simpl1, bool enable_annot)
                    : _Region(start1, end1),
                      _Simplifiable(simpl1),
                      annot_(annot1),
                      gObj_(gObj1),
                      ir_flag_(enable_annot) {
                flt_count_  = 0 ;
                read_rates_ptr_ = nullptr;
                nxbin_      = 0 ;
                nxbin_off_  = 0 ;
                nxbin_mod_  = 0 ;
                markd_      = false ;
                numbins_    = 0 ;
                const_donor_    = false ;
                const_acceptor_ = false ;
            }

            bool    get_annot()             { return annot_ ; }
            Gene*   get_gene()              { return gObj_ ; }
            bool    get_ir_flag()           { return ir_flag_ ; }
            coord_key_t get_key() { return std::make_pair(start_, end_); }
            int     get_nxbin()             { return nxbin_ ; }
            int     get_nxbin_off()         { return nxbin_off_ ; }
            int     get_nxbin_mod()         { return nxbin_mod_ ; }
            int     get_numbins()           { return numbins_ ; }
            int     get_fltcount()          { return flt_count_ ; }
            bool    get_constitutive()      { return const_donor_ && const_acceptor_ && ir_flag_ && !simpl_fltr_ ; }

            void    set_const_donor()       { const_donor_ = true ; }
            void    set_const_acceptor()    { const_acceptor_ = true ; }
            void    set_markd()             { markd_ = true ; }
            void    unset_markd()           { markd_ = false ; }
            bool    is_connected()          { return markd_ ; }
            string  get_key(Gene * gObj) ;
            void    calculate_lambda() ;

            /**
             * Adjusts read rates relative to intron and (max) read length.
             *
             * @note calling more than once currently causes normalization to
             * be applied twice; that is, we do not maintain state as to
             * whether normalization has been performed previously
             */
            void normalize_readrates();

            /**
             * Return boolean indicating if intron is "reliable"
             *
             * @note an intron is considered reliable if it is annotated or
             * previously passed group filters, or if it has sufficient
             * positions with nonzero coverage
             *
             * @params min_bins the percentage of bins that need nonzero
             * coverage distributed over them (in some way)
             * @param eff_len effective length of read to use in propagating
             * coverage to nearby bins used for short introns where only a
             * handful of reads could potentially cover the entire intron
             *
             * @return boolean indicating if intron is reliable
             */
            bool is_reliable(float min_bins, int eff_len);

            void add_read_rates_buff(int eff_len){
                // set up buffer and constants to bin raw positions to eff_len equivalent junction positions
                nxbin_      = (int) ((length() + eff_len) / eff_len);  // minimum number of raw positions per equivalent junction position
                nxbin_mod_  = length() % eff_len;  // number of equivalent positions that have 1 more raw position than nxbin_
                nxbin_off_  = nxbin_mod_ * (nxbin_ + 1);  // offset of raw position -- before: nxbin_+1 per position, after: nxbin_ per position
                numbins_    = eff_len;  // number of equivalent junction positions (keep value)
                read_rates_ptr_ = make_shared<vector<float>>(numbins_, 0);  // allocate per-position read rate buffer as zeroes
            }

            void initReadCovFromVector(vector<float>& cov){
                for (unsigned int i=0; i<cov.size(); ++i) {
                    (*read_rates_ptr_)[i] = cov[i];
                }

            }

            /**
             * Add s reads to intron at the given offset with respect to intron start
             *
             * @param intron_offset genomic/read offset of intron start to last
             * admissible position on read
             * @param eff_len number of bins for different positions
             * @param s number of reads to add
             */
            void add_read(int intron_offset, int eff_len, int s){
                if (intron_offset < 0 || intron_offset >= length() + eff_len) {
                    cerr << "Intron offset " << intron_offset << " miscalculated\n";  // XXX
                    return;
                }
                if (read_rates_ptr_ == nullptr) {
                    add_read_rates_buff(eff_len);
                }
                // get offset over the eff_len bins using intron_offset
                int offset = 0;
                if (intron_offset < nxbin_off_) {
                    // nxbin_ + 1 positions per bin
                    offset = (int) intron_offset / (nxbin_ + 1);
                } else {
                    // nxbin_ positions per bin after nxbin_mod_
                    offset = (int) ((intron_offset - nxbin_off_) / nxbin_) + nxbin_mod_;
                }
                (*read_rates_ptr_)[offset] += s;
            }

            inline void  update_flags(float min_coverage, int min_exps, float min_bins) {
                int cnt = 0 ;
                const float pc_bins = numbins_ * min_bins ;
                for(int i =0 ; i< numbins_; i++){
                    cnt += ((*read_rates_ptr_)[i] >= min_coverage) ? 1 : 0;
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
                read_rates_ptr_ = inIR_ptr->read_rates_ptr_;
                flt_count_ += inIR_ptr->get_fltcount() ;
                ir_flag_    = ir_flag_ || inIR_ptr->get_ir_flag() ;
                return ;
            }

            void clear_nreads(bool reset_grp){
                flt_count_ = reset_grp ? 0: flt_count_ ;
                return ;
            }

            void free_nreads(){
                read_rates_ptr_ = nullptr;
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
            map <coord_key_t, Junction*> junc_map_;
            map <coord_key_t, Exon*> exon_map_;
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
            void    create_annot_intron(int start_ir, int end_ir, bool simpl, bool enable_anot_ir){
                Intron * ir = new Intron(start_ir +1, end_ir -1, true, this, simpl, enable_anot_ir) ;
                intron_vec_.push_back(ir) ;
            }

            void reset_flags(){
                // clear reads, minreads, minpos for next group
                for (const auto &p: junc_map_){
                    (p.second)->clear_nreads(true) ;
                }
                for (const auto &ir: intron_vec_){
                    ir->clear_nreads(true) ;
                }
            }
            void reset_exons(){
                for (auto p  = exon_map_.begin(); p!= exon_map_.end();){
                    map <coord_key_t, Exon*>::iterator pit = p ;
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
            void    detect_exons() ;
            void    connect_introns() ;
            /**
             * Fill intronlist with potential introns between exons
             *
             * @note uses annotated exon boundaries to define potential introns
             */
            void    detect_introns(vector<Intron*> &intronlist, bool simpl) ;
            void    add_intron(Intron * inIR_ptr, float min_coverage, unsigned int min_exps, float min_bins, bool reset,
                               bool denovo) ;
            void    newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db) ;
            void    fill_junc_tlb(map<string, vector<string>> &tlb) ;

            /**
             * Output nonredundant (strict) vs non-mutually-redundant LSVs to
             * out_lsvlist, return number of LSVS added
             */
            int detect_lsvs(vector<LSV*> &out_lsvlist, bool lsv_strict, bool only_source, bool only_target);
            int     get_constitutive_junctions(vector<string>& v) ;
            void    simplify(map<string, int>& junc_tlb, float simpl_percent, int strandness, int denovo_simpl,
                             int db_simple, int ir_simpl, bool last, unsigned int min_experiments) ;
            void    initialize_junction(coord_key_t key, int start, int end, shared_ptr<vector<float>> nreads_ptr, bool simpl) ;
            void    update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                      unsigned int denovo_thresh, unsigned int min_experiments, bool denovo) ;
            void    updateFlagsFromJunc(coord_key_t key, unsigned int sreads, unsigned int minreads_t, unsigned int numpos,
                                        unsigned int minpos_t, unsigned int denovo_t, bool denovo, int minexp,
                                        bool reset) ;

    };

    class LSV{
        private:
            lsv_id_t id_;
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
            inline void get_variations(set<coord_key_t> &t1) ;
            string      get_type ()                  { return type_ ; }
            lsv_id_t get_id() { return id_; }
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

    /**
     * Group of genes over a region. Region is union over genes. Typically
     * constructed to partition genes on chromosome into directly/indirectly
     * overlapping gene lists
     */
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
    void find_gene_from_junc(map<string, vector<overGene*>> & glist, string chrom, char strand, int start, int end,
                             vector<Gene*>& oGeneList, bool ir, bool simpl) ;
    void fill_junc_tlb(vector<LSV*>& lsv_list, map<string, int>& tlb) ;
    bool isNullJinfo(Jinfo* x) ;
    void free_JinfoVec(vector<Jinfo*>& jvec);
    string key_format(string gid, int coord1, int coord2, bool ir) ;
    void free_lsvlist(vector<LSV*> & lsvList) ;
}

#endif

