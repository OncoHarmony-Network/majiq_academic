#ifndef JUNCTION_H
#define JUNCTION_H

#include <iostream>


using namespace std;
struct Junction{
    public:
        unsigned int nreads;
        int index;
        int lsv_index;
//        bool annot;
//        string chrom;
//        char strand;
        string gene_id;
        unsigned int start;
        unsigned int end;
        unsigned int n_exps;

    public:
        Junction() {
            start = 0;
            end = 0;
            nreads = 0;

        }

        Junction(string gene_id1, unsigned int start1, unsigned int end1): gene_id(gene_id1),
                                                                           start(start1), end(end1){
            nreads = 0;
            n_exps = 0;

        }


        void update_junction_read(int n){
            nreads += n;
            n_exps = 0;
            return;
           }

        void reset(){
            nreads = 0;
            index = 0;
            lsv_index = 0;
            return;
        }

//    # @property
//    # def start(self):
//    #     return self.c_iobam.start
//    #
//    # @property
//    # def end(self):
//    #     return self.c_iobam.end
//    #
//    # @property
//    # def nreads(self):
//    #     return self.c_iobam.nreads
//    #
//    # def annot(self):
//    #     return self.c_iobam.annot
//
};

#endif