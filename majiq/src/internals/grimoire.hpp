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

    class Junction{
        public:
            int start ;
            int end ;
            bool annot ;
            bool intronic ;
            unsigned int n_exps ;
            vector<unsigned int> nreads_ ;
            map<unsigned int, unsigned int> pos_map_ ;
//            vector<vector<unsigned int>> nreads_ ;

//            map<pair<unsigned int, unsigned int>, unsigned int> pos_map_ ;
            Exon * acceptor;
            Exon * donor;

            Junction() {
                start = 0 ;
                end = 0 ;
                intronic = false ;

            }

//            Junction(int start1, int end1, int nexp, int eff_len): start(start1), end(end1), n_exps(nexp){
//                intronic = false ;
//            }
            Junction(int start1, int end1): start(start1), end(end1){

                intronic = false ;
            }
            Junction(int start1, int end1, bool annot1): start(start1), end(end1), annot(annot1){} ;

            ~Junction(){
                clear_nreads() ;
            }

            string get_key(){
                return(to_string(start) + "-" + to_string(end)) ;
            }

            void clear_nreads(){
                nreads_.clear() ;
                pos_map_.clear() ;
            }

//            void update_junction_read(int read_start, unsigned int exp_index, unsigned int n) ;
            void update_junction_read(int read_start, unsigned int n) ;
            int length() ;
    };

    class Gene{
        public:
            string id ;
            string name ;
            string chromosome ;
            char strand ;
            int start ;
            int end ;
            map<string, Junction*> junc_map ;
            map<string, Exon*> exon_map ;

            Gene(){}

            Gene(string id1, string name1, string chromosome1,
                 char strand1, unsigned int start1, unsigned int end1): id(id1), name(name1), chromosome(chromosome1),
                                                                             strand(strand1), start(start1), end(end1){}
            void newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db);

            ~Gene(){
                for(const auto &p2: exon_map){
                    delete p2.second ;
                }

                for(const auto &p1: junc_map){
                    delete p1.second ;
                }
            }

            void print_gene(){
                cout<< id << chromosome <<":"<< start <<"-" << end << "\n";
            }

//            void add_elements(map<string, Junction*> &junc_map1, map<string, Exon*> &exon_map1){
//                junc_map = junc_map1 ;
//                exon_map = exon_map1 ;
//            }

            void detect_exons();
            string get_region();
    };

    class Exon{
        public:
            int start ;
            int end ;
            bool annot ;
            bool intron ;
            int db_start ;
            int db_end ;
            set<Junction *> ib ;
            set<Junction *> ob ;

            Exon(){}

            Exon(int start1, int end1): start(start1), end(end1), db_start(start1), db_end(end1){
                annot = true ;
                intron = false ;
            }

            Exon(int start1, int end1, bool annot1): start(start1), end(end1), annot(annot1),
                                                                       db_start(start1), db_end(end1){
                intron = false;
            }

            Exon(int start1, int end1, bool annot1, bool intron1): start(start1), end(end1),
                                                                                     annot(annot1), intron(intron1),
                                                                                     db_start(start1), db_end(end1){}

            ~Exon(){}
    };

    struct Ssite{

        int coord ;
        bool donor_ss ;
        Junction * j ;

    };

    class LSV{
        public:
            string id;
            set<Junction *> junctions;
            string type;

            LSV(){}
            LSV(string gene_id1, char strand, Exon* ex, bool ss){
                bool b = (strand == '+') ;
                string t = (ss != b) ? "t" : "s" ;
                id = gene_id1 + ":" + t + ":" + to_string(ex->start) + "-" + to_string(ex->end) ;
                junctions = ss? ex->ob: ex->ib ;
            }
            ~LSV(){}
    };



    int detect_lsvs(list<LSV*> &out_lsvlist, Gene * gObj, unsigned int minexp, int minpos, int minreads);
    float* boostrap_samples(LSV * lsvObj, int msamples, int ksamples, int exp_idx, int eff_len) ;
    void sortGeneList(vector<Gene*> &glist) ;
}
#endif