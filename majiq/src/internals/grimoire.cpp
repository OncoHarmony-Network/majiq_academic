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

    void Junction::update_junction_read(int read_start, unsigned int n){

        if (pos_map_.count(read_start)>0){
            nreads_[pos_map_[read_start]] += n ;
        }else{

            pos_map_[read_start] =  nreads_.size() ;
            nreads_.push_back(n) ;
        }
        return ;
    }

    int Junction::length(){
        return end - start ;

    }

    bool sort_ss(const Ssite &a, const Ssite &b){


        bool crd =  (a.coord<b.coord) ;
        bool r = (a.coord == b.coord) && (a.donor_ss > b.donor_ss);
        bool same = (a.coord == b.coord) && (a.donor_ss == b.donor_ss) && ((a.j->length()) < b.j->length()) ;
        return  crd || r || same ;
    }

    string Gene::get_region(){

        string reg ;
        stringstream s1, s2 ;
        s1<<start; s2<<end ;
        reg = chromosome + ":" + s1.str() + "-" + s2.str() ;
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



    void Gene::newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db){

        Exon *ex1 = nullptr, *ex2 = nullptr ;
        string key ;
        stringstream s1 ;
        stringstream s2 ;
        if ((end - start) < 1){
            return;
        }

        ex1 = exonOverlap(exon_map, start, end) ;
        if (nullptr != inbound_j && nullptr != outbound_j && inbound_j->intronic && outbound_j->intronic) {
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


    void Gene::detect_exons(){
        vector<Ssite> ss_vec ;
        vector<Junction *> opened_exon ;
        Junction * last_5prime = nullptr, *first_3prime = nullptr ;
        vector<Ssite>::iterator ss_vecIt ;
        vector<Junction *>::iterator exonIterator ;
//        print_gene() ;
        for (const auto &jnc :this->junc_map){
//            cout << "JUNC: " << jnc.first << " : " << (jnc.second) << " :"<< (jnc.second)->start << "-" << (jnc.second)->end << "\n" ;

            if ((jnc.second)->start > 0) {
                Ssite s = {(jnc.second)->start, 1, jnc.second} ;
                ss_vec.push_back(s) ;
//                cout<< "SS1:" << "{" << s.coord << ", " << s.donor_ss << ", " << s.j << " :: " << s.j->length() << "}\n" ;
            }
            if ((jnc.second)->end > 0) {
                Ssite s = {(jnc.second)->end, 0, jnc.second} ;
//                cout<< "SS2:" << "{" << s.coord << ", " << s.donor_ss << ", " << s.j << " :: " << s.j->length() << "}\n" ;
                ss_vec.push_back(s) ;
            }

        }
        sort(ss_vec.begin(), ss_vec.end(), sort_ss) ;
        for(ss_vecIt = ss_vec.begin(); ss_vecIt != ss_vec.end(); ss_vecIt++){
    //        #TODO: check this in advance
    //        if not (ss_listIterator->j.is_reliable() or spl.j.annot):
    //            continue

            if (ss_vecIt->donor_ss) {
                if (opened_exon.size() > 0){
                    newExonDefinition(opened_exon.back()->end, ss_vecIt->coord, opened_exon.back(), ss_vecIt->j, false) ;
                    opened_exon.pop_back() ;


                } else if ( 0 == opened_exon.size()) {
                    if (nullptr == first_3prime){
                        newExonDefinition((ss_vecIt->j)->start-10, (ss_vecIt->j)->start, nullptr, ss_vecIt->j, false) ;
                    }else{
                        newExonDefinition(first_3prime->end, ss_vecIt->coord, first_3prime, ss_vecIt->j, false) ;
                    }
                }
                last_5prime = ss_vecIt->j ;
            }else{
                if (opened_exon.size() > 0){
                    if (nullptr != last_5prime){
                        for(exonIterator = opened_exon.begin(); exonIterator != opened_exon.end(); exonIterator++){
                            newExonDefinition((*exonIterator)->end, last_5prime->start, (*exonIterator),
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


            }
        }


        for(exonIterator = opened_exon.begin(); exonIterator != opened_exon.end(); exonIterator++){
            newExonDefinition((*exonIterator)->end, (*exonIterator)->end+10, (*exonIterator), nullptr, false) ;
        }
        ss_vec.clear() ;
        ss_vec.shrink_to_fit();
        return ;
    }

//    inline bool is_lsv(set<Junction*> &juncSet, unsigned int nexp, unsigned int eff_len, int minpos, int minreads){
//        set<Junction*>::iterator juncIt ;
//        unsigned int njuncs = 0 ;
//        for (juncIt = juncSet.begin(); juncIt != juncSet.end(); juncIt++){
//            if (FIRST_LAST_JUNC != (*juncIt)->start) {
//                bool pass = false ;
//                njuncs += 1 ;
//                for (unsigned int i=0; i<= nexp; i++){
//                    const int npos =   (*juncIt)->nreads_[i].size();
//                    int nreads = 0;
//                    for (auto& n : (*juncIt)->nreads_[i])
//                        nreads += n;
//
//                    pass = pass || (npos >= minpos && nreads >= minreads) ;
//                    if (njuncs >1 && pass){
//                        return true ;
//                    }
//                }
//            }
//        }
//        return false;
//    }

    int detect_lsvs(list<LSV*> &out_lsvlist, Gene* gObj, unsigned int nexp, unsigned int eff_len, int minpos, int minreads){

        int count = 0 ;
//        map<string, Exon*>::iterator exon_mapIt ;
//
//        list<LSV*> lsv_list ;
//        set<pair<set<string>, LSV*>> source ;
//        LSV * lsvObj ;
//        for (exon_mapIt = gObj->exon_map.begin(); exon_mapIt != gObj->exon_map.end(); ++exon_mapIt) {
//            Exon * ex = exon_mapIt->second ;
//            set<Junction*>::iterator juncIt ;
//
//            if (is_lsv(ex->ob, nexp, eff_len, minpos, minreads)) {
//                lsvObj = new LSV(gObj->id, gObj->strand, ex, true) ;
//                set<Junction *>::iterator iter ;
//                set<string> t1 ;
//                for (iter = (lsvObj->junctions).begin(); iter != (lsvObj->junctions).end(); iter ++){
//                    t1.insert((*iter)->get_key()) ;
//                }
//                pair<set<string>, LSV*> _p1 (t1, lsvObj) ;
//                source.insert(_p1) ;
//                lsv_list.push_back(lsvObj) ;
//            }
//
//            if (is_lsv(ex->ib, nexp, eff_len, minpos, minreads)) {
//                lsvObj = new LSV(gObj->id, gObj->strand, ex, false) ;
//                set<Junction *>::iterator iter ;
//                set<string> t1 ;
//                for (iter = (lsvObj->junctions).begin(); iter != (lsvObj->junctions).end(); iter ++){
//                    t1.insert((*iter)->get_key()) ;
//                }
//
//                set<pair<set<string>, LSV*>> ::iterator setIter ;
//                for (setIter = source.begin(); setIter != source.end(); setIter++){
//                    set<string> d1 ;
//                    set<string> d2 ;
//
//                    set_difference((setIter->first).begin(), (setIter->first).end(), t1.begin(), t1.end(), inserter(d1, d1.begin())) ;
//                    set_difference(t1.begin(), t1.end(), (setIter->first).begin(), (setIter->first).end(), inserter(d2, d2.begin())) ;
//
//                    if (d1.size()>0 and d2.size()>0){
//                        lsv_list.push_back(lsvObj) ;
//                    } else if (d1.size() >0 or (d1.size()==0 and d2.size() == 0)){
//                        delete lsvObj ;
//
//                    } else if (d2.size() >0){
//                        lsv_list.push_back(lsvObj) ;
//                        lsv_list.remove(setIter->second) ;
//                    }
//                }
//            }
//        }
//
////        list<LSV*>::iterator lsvIter ;
////
////        for (lsvIter=lsv_list.begin(); lsvIter!=lsv_list.end(); lsvIter++){
////            (*lsvIter)->add_lsv() ;
////            count += 1 ;
//////            delete lsvIter ;
////        }
        return count ;

    }


}