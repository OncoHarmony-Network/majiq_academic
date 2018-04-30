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


using namespace std ;

namespace grimoire{
    class Exon;
    class Gene;

    class Junction{
        private:
            int start_ ;
            int end_ ;
            bool annot_ ;
            bool intronic_ ;
            bool bld_fltr_ ;
            bool denovo_bl_ ;
            unsigned int denovo_cnt_ ;
            unsigned int flter_cnt_ ;
            Exon * acceptor_;
            Exon * donor_;

        public:
            int* nreads_ ;

            Junction() {}
            Junction(int start1, int end1, bool annot1): start_(start1), end_(end1), annot_(annot1){
                bld_fltr_ = false ;
                denovo_bl_ = false ;
                denovo_cnt_ = 0 ;
                flter_cnt_ = 0 ;
                intronic_ = false ;
                nreads_ = nullptr ;
            }
            ~Junction()             { clear_nreads(true) ; }

            string  get_key(Gene * gObj) ;
            string  get_key()       { return(to_string(start_) + "-" + to_string(end_)) ; }
            int     get_start()     { return start_ ; }
            int     get_end()       { return end_ ; }
            bool    get_annot()     { return annot_ ; }
            bool    get_intronic()  { return intronic_ ; }
            bool    get_bld_fltr()  { return bld_fltr_ ; }
            bool    get_denovo_bl() { return denovo_bl_ ; }
            Exon*   get_acceptor()  { return acceptor_ ; }
            Exon*   get_donor()     { return donor_ ; }

            void  set_nreads_ptr(int * nreads1) { nreads_ = nreads1 ; }
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
                if (sum_reads >= denovo_thresh || annot_){
                    ++ denovo_cnt_  ;
                    denovo_bl_ = denovo_bl_ || (denovo_cnt_ >= min_experiments) ;
                }

cout << get_key() << " : sumreads npos: " << sum_reads<< " : "<< npos << " denovo_thresh: "<<denovo_thresh

<< " : annot_: " <<  annot_ << "denovo::"<<denovo_bl_<<"\n" ;
                return ;
            }

            void clear_nreads(bool reset_grp){
                nreads_ = nullptr ;
                flter_cnt_ = reset_grp ? 0: flter_cnt_ ;
                denovo_cnt_ = reset_grp ? 0: denovo_cnt_ ;
            }

            void update_junction_read(int read_start, unsigned int n) ;
            int length() ;
    };
    class Exon{

        public:
            int     start_ ;
            int     end_ ;
            bool    annot_ ;
            bool    intron_ ;
            int     db_start_ ;
            int     db_end_ ;
            set<Junction *> ib ;
            set<Junction *> ob ;

            Exon(){}
            Exon(int start1, int end1): start_(start1), end_(end1), db_start_(start1), db_end_(end1){
                annot_ = true ;
                intron_ = false ;
            }
            Exon(int start1, int end1, bool annot1): start_(start1), end_(end1), annot_(annot1),
                                                                       db_start_(start1), db_end_(end1){
                intron_ = false;
            }
            Exon(int start1, int end1, bool annot1, bool intron1): start_(start1), end_(end1),
                                                                   annot_(annot1), intron_(intron1),
                                                                   db_start_(start1), db_end_(end1){}
            ~Exon(){}

            int     get_start()             { return start_ ; }
            int     get_end()               { return end_ ; }
            bool    get_intron()            { return intron_ ; }

            void    set_start(int start1)   { start_ = start1 ; }
            void    set_end(int end1)       { end_ = end1 ; }
    };

    class Gene {
        private:
            string id_ ;
            string name_ ;
            string chromosome_ ;
            char strand_ ;
            int start_ ;
            int end_ ;

        public:
            map<string, Junction*> junc_map_ ;
            map<string, Exon*> exon_map_ ;


            Gene(){}
            Gene(string id1, string name1, string chromosome1,
                 char strand1, unsigned int start1, unsigned int end1): id_(id1), name_(name1), chromosome_(chromosome1),
                                                                             strand_(strand1), start_(start1), end_(end1){}
            ~Gene(){
                for(const auto &p2: exon_map_){
                    delete p2.second ;
                }
                for(const auto &p1: junc_map_){
                    delete p1.second ;
                }
            }
            int     get_start()     { return start_ ; }
            int     get_end()       { return end_ ; }
            string  get_id()        { return id_ ; }
            string  get_chromosome(){ return chromosome_ ;}
            char    get_strand()    { return strand_ ;}
            string  get_name()      { return name_ ;}
            string  get_region() ;
            void    print_exons() ;
            void    detect_exons() ;

            void    newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db) ;
            void    fill_junc_tlb(map<string, vector<string>> &tlb) ;
            void    update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                      unsigned int denovo_thresh, unsigned int min_experiments) ;
    };
//            void    print_gene()    { cout<< id_ << chromosome_ <<":"<< start_ <<"-" << end_ << "\n"; }



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


    class LSV{
        private:
            string id_ ;
            vector<Junction *> junctions_ ;
            string type_ ;
            Gene* gObj_ ;

        public:
            LSV(){}
            LSV(Gene* gObj1, Exon* ex, bool ss): gObj_(gObj1){
                const bool b = (gObj_->get_strand() == '+') ;
                string t = (ss != b) ? "t" : "s" ;
                id_ = gObj_->get_id() + ":" + t + ":" + to_string(ex->get_start()) + "-" + to_string(ex->get_end()) ;
                type_ = set_type(ex, ss) ;
            }
            ~LSV(){}
            bool gather_lsv_info (float* source, float* target, list<Jinfo*> &info, map<string, Jinfo> &tlb, unsigned int msample) ;
            string set_type (Exon* ex, bool ss) ;
            string get_type ()                  { return type_ ; }
            string get_id ()                    { return id_ ; }
            Gene*  get_gene()                   { return gObj_ ; }
            int    get_num_junctions()          { return junctions_.size() ; }
            vector<Junction *>& get_junctions() { return junctions_ ; }
    };

    bool is_lsv(set<Junction*> &juncSet, bool ss) ;
    int detect_lsvs(vector<LSV*> &out_lsvlist, Gene * gObj);
    void sortGeneList(vector<Gene*> &glist) ;
}
#endif
