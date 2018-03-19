#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <list>
#include <set>
#include <random>
#include <algorithm>
#include "grimoire.hpp"

using namespace std;

namespace grimoire {

    float* boostrap_samples(LSV * lsvObj, int msamples, int ksamples, int exp_idx, int eff_len){
        float * boots = (float*) calloc(lsvObj->junctions.size() * eff_len, sizeof(double)) ;
        int jidx = 0 ;
        for (const auto &jnc: lsvObj->junctions){
            vector<int> jj ;
            for (int idx=0; idx<eff_len; idx){
                const int v = jnc->nreads[exp_idx, idx] ;
                if (v > 0) jj.push_back(v) ;
            }
            default_random_engine generator;
            uniform_int_distribution<int> distribution(0,jj.size());

            for (int m=0; m<msamples; m++){
                float lambda = 0;
                for (int k=0; k<ksamples; k++) lambda += distribution(generator) ;
                lambda /= ksamples ;
                boots[jidx, m] = lambda * jj.size() ;
            }
            jidx++ ;
        }

    }

    void Junction::update_junction_read(int read_start, unsigned int exp_idx, unsigned int n){
        int offs = start  - (read_start + 8) ;
//               cout << "vals " << exp_idx<< " " <<offs << "::" <<  nreads[exp_idx, offs] <<"\n" ;
        nreads[exp_idx, offs] += n ;
        return ;
    }

    bool sort_ss(const Ssite &a, const Ssite &b){
        return (a.coord<b.coord) || (a.coord == b.coord && a.donor_ss) ;
    }

    string Gene::get_region(){

        string reg ;
        stringstream s1, s2 ;

        s1<<start; s2<<end ;
        reg = chromosome + ":" + s1.str() + "-" + s2.str() ;
//        cout << "GENE REG  " << reg << "\n";
        return reg ;

    }

    Exon* addExon (map<string, Exon*> &exon_map, int start, int end, bool in_db){

        Exon * ex ;
        string key ;
        stringstream s1 ;
        stringstream s2 ;

        unsigned int start1 = (start == EMPTY_COORD) ? end - 10 : start ;
        unsigned int end1 = (end == EMPTY_COORD) ? start + 10 : end ;

        s1 << start1 ; s2 << end1 ;
        key = s1.str() + "-" + s2.str() ;

        if (exon_map.count(key) > 0){
            ex = exon_map[key];
        } else {
            ex = new Exon(start, end, in_db) ;
            exon_map[key] = ex ;
        }
        return (ex) ;
    }

    inline Exon* exonOverlap(map<string, Exon*> &exon_map, int start, int end){

        map<string, Exon*>::iterator exon_mapIt ;
        for (exon_mapIt = exon_map.begin(); exon_mapIt != exon_map.end(); ++exon_mapIt) {
            Exon *x = exon_mapIt->second ;
            if (x->start != EMPTY_COORD && x->end != EMPTY_COORD && start < x->end && end > x->start) {
                return x ;
            }
        }
        return nullptr;
    }



    void newExonDefinition(int start, int end, map<string, Exon*> &exon_map,
                           Junction *inbound_j, Junction *outbound_j, bool in_db){

        Exon *ex1 = nullptr, *ex2 = nullptr ;
        string key ;
        stringstream s1 ;
        stringstream s2 ;
        if ((end - start) < 1){
            return;
        }

        ex1 = exonOverlap(exon_map, start, end) ;

        if (inbound_j->intronic && outbound_j->intronic) {
            s1 << start ; s2 << end ;
            key = s1.str() + "-" + s2.str() ;
            if (exon_map.count(key) > 0){
                ex1 = exon_map[key] ;
                ex2 = ex1 ;
            } else {
                return ;
            }
        }else if (nullptr == ex1){

            if ((end - start) <= MAX_DENOVO_DIFFERENCE){
                s1 << start ; s2 << end ;
                key = s1.str() + "-" + s2.str() ;
                ex1 = new Exon(start, end, in_db) ;
                ex2 = ex1 ;
                exon_map[key] = ex1 ;

            } else {
                ex1 = addExon(exon_map, start, EMPTY_COORD, in_db) ;
                ex2 = addExon(exon_map, EMPTY_COORD, end, in_db) ;
            }

        } else {
            ex2 = ex1;
            if (start < (ex1->start - MAX_DENOVO_DIFFERENCE)){
                ex1 = addExon(exon_map, start, EMPTY_COORD, in_db) ;
            } else if (start < ex1->start){
                ex1->start = start ;
            }

            if (end > (ex2->end + MAX_DENOVO_DIFFERENCE)){
                ex2 = addExon(exon_map, EMPTY_COORD, end, in_db) ;
            } else if (end > ex2->end){
                ex2->end = end ;
            }
        }
        if (nullptr != inbound_j){
            (ex1->ib).insert(inbound_j) ;
            inbound_j->acceptor = ex1 ;
        }
        if (outbound_j != nullptr){
            (ex2->ob).insert(outbound_j) ;
            outbound_j->donor = ex2 ;
        }

        return ;
    }


    void detect_exons(map<string, Junction*> &junc_map, map<string, Exon*> &exon_map){
//        Ssite spl ;
        int opened = 0 ;
        vector<Ssite> ss_vec ;
        vector<Junction *> opened_exon ;
        Junction * last_5prime = nullptr, *first_3prime = nullptr ;
        vector<Ssite>::iterator ss_vecIt ;
        vector<Junction *>::iterator exonIterator ;

        for (const auto &j :junc_map){
            if ((j.second)->start > 0) {
                Ssite s = {(unsigned int)(j.second)->start, true, j.second} ;
                ss_vec.push_back(s) ;
            }
            if ((j.second)->end > 0) {
                Ssite s = {(unsigned int)(j.second)->end, false, j.second} ;
                ss_vec.push_back(s) ;
            }
        }
        sort(ss_vec.begin(), ss_vec.end(), sort_ss) ;

        for(ss_vecIt = ss_vec.begin(); ss_vecIt != ss_vec.end(); ss_vecIt++){
    //        #TODO: check this in advance
    //        if not (ss_listIterator->j.is_reliable() or spl.j.annot):
    //            continue

            if (ss_vecIt->donor_ss) {
                if (opened > 0){
                    newExonDefinition(opened_exon[opened-1]->end, ss_vecIt->coord, exon_map, opened_exon[-1],
                                      ss_vecIt->j, false) ;
                    opened_exon.pop_back() ;
                    opened -= 1 ;

                } else if ( 0 == opened) {
                    if (nullptr == first_3prime){
                        newExonDefinition((ss_vecIt->j)->start-10, (ss_vecIt->j)->start, exon_map, nullptr, ss_vecIt->j,
                                          false) ;
                    }else{
                        newExonDefinition(first_3prime->end, ss_vecIt->coord, exon_map, first_3prime,
                                          ss_vecIt->j, false) ;
                    }
                }
                last_5prime = ss_vecIt->j ;

            }else{
                if (opened > 0){
                    if (nullptr != last_5prime){
                        for(exonIterator = opened_exon.begin(); exonIterator != opened_exon.end(); exonIterator++){
                            newExonDefinition((*exonIterator)->end, last_5prime->start, exon_map, (*exonIterator),
                                              last_5prime, false) ;
                        }
                        last_5prime = nullptr ;
                        opened_exon.clear() ;
                        first_3prime = (ss_vecIt->j) ;
                    }
                } else {
                    last_5prime = nullptr ;
                    first_3prime = (ss_vecIt->j) ;
                }
                opened_exon.push_back(ss_vecIt->j) ;
                opened += 1 ;

            }
        }

        for(exonIterator = opened_exon.begin(); exonIterator != opened_exon.end(); exonIterator++){
            newExonDefinition((*exonIterator)->end, (*exonIterator)->end+10, exon_map, (*exonIterator), nullptr, false) ;
        }

        ss_vec.clear() ;

        return ;
    }

    inline bool is_lsv(set<Junction*> &juncSet, unsigned int nexp, unsigned int eff_len, int minpos, int minreads){
        set<Junction*>::iterator juncIt ;
        unsigned int njuncs = 0 ;
        for (juncIt = juncSet.begin(); juncIt != juncSet.end(); juncIt++){
            if (FIRST_LAST_JUNC != (*juncIt)->start) {
                unsigned int * cov = (*juncIt)->nreads ;
                bool pass = false ;
                njuncs += 1 ;
                for (unsigned int i=0; i<= nexp; i++){
                    int npos = 0 ;
                    int nreads = 0 ;

                    for (unsigned int j = 0; j<= eff_len; j++){
                        npos += (int) (cov[i,j] >0) ;
                        nreads += cov[i,j] ;
                    }
                    pass = pass || (npos >= minpos && nreads >= minreads) ;
                    if (njuncs >1 && pass){
                        return true ;
                    }
                }
            }
        }
        return false;
    }

    int detect_lsvs(list<LSV*> &out_lsvlist, map<string, Exon*> &exon_map, Gene* gObj, unsigned int nexp, unsigned int eff_len,
                    int minpos, int minreads){

        int count = 0 ;
        map<string, Exon*>::iterator exon_mapIt ;
        list<LSV*> lsv_list ;
        set<pair<set<string>, LSV*>> source ;
        LSV * lsvObj ;

        for (exon_mapIt = exon_map.begin(); exon_mapIt != exon_map.end(); ++exon_mapIt) {
            Exon * ex = exon_mapIt->second ;
            set<Junction*>::iterator juncIt ;

            if (is_lsv(ex->ob, nexp, eff_len, minpos, minreads)) {
                lsvObj = new LSV(gObj->id, gObj->strand, ex, true) ;
                set<Junction *>::iterator iter ;
                set<string> t1 ;
                for (iter = (lsvObj->junctions).begin(); iter != (lsvObj->junctions).end(); iter ++){
                    t1.insert((*iter)->get_key()) ;
                }
                pair<set<string>, LSV*> _p1 (t1, lsvObj) ;
                source.insert(_p1) ;
                lsv_list.push_back(lsvObj) ;
            }

            if (is_lsv(ex->ib, nexp, eff_len, minpos, minreads)) {
                lsvObj = new LSV(gObj->id, gObj->strand, ex, false) ;
                set<Junction *>::iterator iter ;
                set<string> t1 ;
                for (iter = (lsvObj->junctions).begin(); iter != (lsvObj->junctions).end(); iter ++){
                    t1.insert((*iter)->get_key()) ;
                }

                set<pair<set<string>, LSV*>> ::iterator setIter ;
                for (setIter = source.begin(); setIter != source.end(); setIter++){
                    set<string> d1 ;
                    set<string> d2 ;

                    set_difference((setIter->first).begin(), (setIter->first).end(), t1.begin(), t1.end(), inserter(d1, d1.begin())) ;
                    set_difference(t1.begin(), t1.end(), (setIter->first).begin(), (setIter->first).end(), inserter(d2, d2.begin())) ;

                    if (d1.size()>0 and d2.size()>0){
                        lsv_list.push_back(lsvObj) ;
                    } else if (d1.size() >0 or (d1.size()==0 and d2.size() == 0)){
                        delete lsvObj ;

                    } else if (d2.size() >0){
                        lsv_list.push_back(lsvObj) ;
                        lsv_list.remove(setIter->second) ;
                    }
                }
            }
        }

//        list<LSV*>::iterator lsvIter ;
//
//        for (lsvIter=lsv_list.begin(); lsvIter!=lsv_list.end(); lsvIter++){
//            (*lsvIter)->add_lsv() ;
//            count += 1 ;
////            delete lsvIter ;
//        }
        return count ;

    }

    void free_gene(Gene * gObj, map<string, Junction*> &junc_map, map<string, Exon*> &exon_map){
        for(const auto &p1 :junc_map){
            delete p1.second;
        }
        for(const auto &p2 : exon_map){
            delete p2.second;
        }
        junc_map.clear() ;
        exon_map.clear() ;
        delete gObj;
    }



}