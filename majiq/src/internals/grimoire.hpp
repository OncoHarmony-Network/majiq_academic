#ifndef GRIMOIRE_H
#define GRIMOIRE_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <string>

#define EMPTY_COORD  -1
#define FIRST_LAST_JUNC  -2
#define MAX_DENOVO_DIFFERENCE 500
#define MIN_INTRON_BINSIZE 1000

using namespace std ;

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

    struct Jinfo {
        unsigned int    index ;
        unsigned int    start ;
        unsigned int    end ;
        int             sreads ;
        int             npos ;
    };

    struct lsvtype {
        int         coord ;
        int         ref_coord ;
        Exon *      ex_ptr ;
        Junction *  jun_ptr ;
    };


    class _Region{
        protected:
            int start_ ;
            int end_ ;

            _Region(){}
            _Region(int start1, int end1): start_(start1), end_(end1) {}
        public:
            int         get_start()             { return start_ ; }
            int         get_end()               { return end_ ; }
            void        set_start(int start1)   { start_ = start1 ; }
            void        set_end(int end1)       { end_ = end1 ; }
            inline int  length()                { return end_ - start_ ; }

            template <class myRegion>
            static int RegionSearch(vector<myRegion *>  &a, int n, int coord) {
                int l = 0 ;
                int h = n ; // Not n - 1
                while (l < h) {
                    int mid = (l + h) / 2 ;
                    if (a[mid]->get_start() >= coord) {
                        h = mid ;
                    } else {
                        l = mid +1 ;
                    }
                }
                return l-1;
            }

            template <class myRegion>
            static bool islowerRegion(myRegion * a, myRegion * b){
                return (a->get_start() < b->get_start()) ||
                        (a->get_start() == b->get_start() && a->get_end() < b->get_end());
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
            unsigned int denovo_cnt_ ;
            unsigned int flter_cnt_ ;
            Exon * acceptor_;
            Exon * donor_;

        public:
            float* nreads_ ;

            Junction() {}
            Junction(int start1, int end1, bool annot1): _Region(start1, end1), annot_(annot1){
                denovo_bl_ = annot1 ;
                denovo_cnt_ = 0 ;
                bld_fltr_ = false ;
                flter_cnt_ = 0 ;
                intronic_ = false ;
                nreads_ = nullptr ;
            }
            ~Junction()             { clear_nreads(true) ; }

            string  get_key()       { return(to_string(start_) + "-" + to_string(end_)) ; }
            string  get_key(Gene * gObj) ;
            bool    get_annot()     { return annot_ ; }
            bool    get_intronic()  { return intronic_ ; }
            bool    get_bld_fltr()  { return bld_fltr_ ; }
            bool    get_denovo_bl() { return denovo_bl_ ; }
            Exon*   get_acceptor()  { return acceptor_ ; }
            Exon*   get_donor()     { return donor_ ; }

            void  set_nreads_ptr(float * nreads1) { nreads_ = nreads1 ; }
            void  set_acceptor(Exon * acc) { acceptor_ = acc ; }
            void  set_donor(Exon * don) { donor_ = don ; }

            void update_flags(int efflen, unsigned int num_reads, unsigned int num_pos, unsigned int denovo_thresh,
                              unsigned int min_experiments){
                if (nreads_ == nullptr ){
                    return ;
                }
                unsigned int sum_reads = 0 ;
                unsigned int npos = 0 ;

                for(int i=0; i<efflen; ++i){
                    sum_reads += nreads_[i] ;
                    npos += nreads_[i]? 1 : 0 ;
                }

                if (( npos >= num_pos) && (sum_reads >= num_reads)){
                    ++ flter_cnt_ ;
                    bld_fltr_ = bld_fltr_ || (flter_cnt_ >= min_experiments) ;
                }
                if (sum_reads >= denovo_thresh){
                    ++ denovo_cnt_  ;
                    denovo_bl_ = denovo_bl_ || (denovo_cnt_ >= min_experiments) ;
                }
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
            bool    has_out_intron()        { return ob_irptr != nullptr ; }

    };

    class Intron: public _Region{
        private:
//            string  id_ ;
            bool    annot_ ;
            Gene *  gObj_ ;
            int     flt_count_ ;
            bool    ir_flag_ ;

        public:
            float *   read_rates_ ;
            int     nbins_ ;

            Intron () {}
            Intron (int start1, int end1, bool annot1, Gene* gObj1): _Region(start1, end1),
                                                                     annot_(annot1), gObj_(gObj1) {
                flt_count_  = 0 ;
                ir_flag_    = false ;
                read_rates_ = nullptr ;
                nbins_      = 0 ;
            }

            bool    get_annot()             { return annot_ ; }
            Gene*   get_gene()              { return gObj_ ; }
            bool    get_ir_flag()           { return ir_flag_ ; }
            string  get_key()               { return(to_string(start_) + "-" + to_string(end_)) ; }
            string  get_key(Gene * gObj) ;
            void    calculate_lambda() ;
            void    add_read_rates_buff(int nbins1){

                if ( end_ - start_ >= nbins1 * MIN_INTRON_BINSIZE ) {
                    nbins_ = nbins1 ;
                } else {
                    nbins_ = (end_ - start_) / MIN_INTRON_BINSIZE ;
                    nbins_ = (nbins_ == 0) ? 1 : nbins_ ;

                }
//                read_rates_ = (int*) calloc(nbins_, sizeof(int)) ;
                read_rates_ = (float*) calloc(nbins1, sizeof(float)) ;
            }
            bool  is_reliable(float min_bins) ;
            inline void  update_flags(float min_coverage, int min_exps) {

                int cnt = 0 ;
                for(int i =0 ; i< nbins_; i++){
                    cnt += (read_rates_[i]>= min_coverage) ? 1 : 0 ;
                }
                flt_count_ += (cnt == nbins_) ? 1 : 0 ;
                ir_flag_ = ir_flag_ || (flt_count_ >= min_exps) ;
            }

            void overlaping_intron(Intron * inIR_ptr){
                start_ = max(start_, inIR_ptr->get_start()) ;
                end_ = min(end_, inIR_ptr->get_end()) ;
                read_rates_ = inIR_ptr->read_rates_ ;
            }
    };


    class Gene: public _Region{
        private:
            string  id_ ;
            string  name_ ;
            string  chromosome_ ;
            char    strand_ ;

        public:
            map <string, Junction*> junc_map_ ;
            map <string, Exon*> exon_map_ ;
            vector <Intron*> intron_vec_ ;


            Gene (){}
            Gene(string id1, string name1, string chromosome1,
                 char strand1, unsigned int start1, unsigned int end1): _Region(start1, end1), id_(id1), name_(name1),
                                                                        chromosome_(chromosome1), strand_(strand1){}
            ~Gene(){
                for(const auto &p2: exon_map_){
                    delete p2.second ;
                }
                for(const auto &p1: junc_map_){
                    delete p1.second ;
                }
            }
            string  get_id()        { return id_ ; }
            string  get_chromosome(){ return chromosome_ ;}
            char    get_strand()    { return strand_ ;}
            string  get_name()      { return name_ ;}

            void    create_annot_intron(int start_ir, int end_ir){
                Intron * ir = new Intron(start_ir, end_ir, true, this) ;
                intron_vec_.push_back(ir) ;
            }

            string  get_region() ;
            void    print_exons() ;
            void    detect_exons() ;
            void    detect_introns(vector<Intron*> &intronlist) ;
            void    add_intron(Intron * inIR_ptr, float min_coverage, unsigned int min_exps) ;
            void    connect_introns() ;
            void    newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db) ;
            void    fill_junc_tlb(map<string, vector<string>> &tlb) ;
            int     detect_lsvs(vector<LSV*> &out_lsvlist);
            void    update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                      unsigned int denovo_thresh, unsigned int min_experiments) ;

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
                id_     = gObj_->get_id() + ":" + t + ":" + to_string(ex->get_start()) + "-" + to_string(ex->get_end()) ;
                ir_ptr_ = ss? ex->ob_irptr : ex->ib_irptr ;
                type_   = set_type(ex, ss) ;
            }

            ~LSV(){}
            bool gather_lsv_info (float* source, float* target, list<Jinfo*> &info, map<string, Jinfo> &tlb, unsigned int msample) ;
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

    void sortGeneList(vector<Gene*> &glist) ;
    vector<Intron *> find_intron_retention(vector<Gene*> &gene_list, int start, int end) ;

}

#endif

