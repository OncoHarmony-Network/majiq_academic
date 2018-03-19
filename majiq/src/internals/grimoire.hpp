#ifndef GRIMOIRE_H
#define GRIMOIRE_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>

#define EMPTY_COORD  -1
#define FIRST_LAST_JUNC  -2
#define MAX_DENOVO_DIFFERENCE 500


using namespace std ;

namespace grimoire{
    class Exon;

    class Junction{
        public:
            string gene_id ;
            int start ;
            int end ;
            bool annot ;
            bool intronic ;
            unsigned int n_exps ;
            unsigned int * nreads ;
            Exon * acceptor;
            Exon * donor;

            Junction() {
                start = 0 ;
                end = 0 ;
                intronic = false ;
                nreads = (unsigned int *) calloc(1, sizeof(unsigned int)) ;
            }
            Junction(string gene_id1, int start1, int end1, int nexp, int eff_len): gene_id(gene_id1), start(start1),
                                                                                               end(end1), n_exps(nexp){
                if (start>0 && end>0) {
//                    cout<< "allocating for " << nexp << " " << eff_len << "\n" ;
                    nreads = (unsigned int *) calloc(nexp * eff_len, sizeof(unsigned int)) ;
//                    cout<< "ALLOC DONE\n" ;
                }
                intronic = false ;

            }
            ~Junction(){
                cout<<"Trying to free junc:" << start << "-" << end<< "::"<<nreads << "\n" ;
                if (start>0 && end>0) {
//                     free(nreads) ;
                }
                cout<< "FREE DONE\n" ;
            }

            string get_key(){
                return(to_string(start) + "-" + to_string(end)) ;
            }

            void update_junction_read(int read_start, unsigned int exp_index, unsigned int n);
    };

    class Gene{
        public:
            string id ;
            string name ;
            string chromosome ;
            char strand ;
            unsigned int start ;
            unsigned int end ;

            Gene(){}

            Gene(string id1, string name1, string chromosome1,
                 char strand1, unsigned int start1, unsigned int end1): id(id1), name(name1), chromosome(chromosome1),
                                                                                                start(start1), end(end1){}
            ~Gene(){}

            string get_region();

    };

    class Exon{
        public:
            int start ;
            int end ;
            bool annot ;
            bool intron ;
            unsigned int db_start ;
            unsigned int db_end ;
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
    };

    struct Ssite{

        unsigned int coord ;
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
                cout << strand << "::"<< to_string(b) ;
                string t = (ss != b) ? "t" : "s" ;
                id = gene_id1 + ":" + t + ":" + to_string(ex->start) + "-" + to_string(ex->end) ;
                junctions = ss? ex->ob: ex->ib ;
            }
            ~LSV(){}
    };



    int detect_lsvs(list<LSV*> &out_lsvlist, map<string, Exon*> &exon_map, Gene * gObj, unsigned int nexp,
                    unsigned int eff_len, int minpos, int minreads);
    void detect_exons(map<string, Junction*> &junc_map, map<string, Exon*> &exon_map);

    void free_gene(Gene * gObj, map<string, Junction*> &junc_map, map<string, Exon*> &exon_map);

    float* boostrap_samples(LSV * lsvObj, int msamples, int ksamples, int exp_idx, int eff_len) ;

}
#endif