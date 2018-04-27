#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <list>
#include <set>
#include <random>
#include <algorithm>
#include <string>
#include "grimoire.hpp"

using namespace std;

namespace grimoire {

    bool sort_ss(const Ssite &a, const Ssite &b){

        bool crd =  (a.coord<b.coord) ;
        bool r = (a.coord == b.coord) && (a.donor_ss > b.donor_ss);
        bool same = (a.coord == b.coord) && (a.donor_ss == b.donor_ss) && ((a.j->length()) < b.j->length()) ;
        return  crd || r || same ;
    }

    string Junction::get_key(Gene * gObj) {
        return(gObj->get_chromosome() + ":" + to_string(start_) + "-" + to_string(end_)) ;
    }

    int Junction::length() { return end_ - start_ ; }


    Exon* addExon (map<string, Exon*> &exon_map, int start, int end, bool in_db){

        Exon * ex ;
        unsigned int start1 = (start == EMPTY_COORD) ? end - 10 : start ;
        unsigned int end1 = (end == EMPTY_COORD) ? start + 10 : end ;
        const string key = to_string(start1) + "-" + to_string(end1) ;

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
            if (   ( x->get_start() != EMPTY_COORD ) && ( x->get_end() != EMPTY_COORD )
                && ( start < x->get_end() ) && ( end > x->get_start()) ) {
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
            return ;
        }

//        cout << "INPUT:: " << inbound_j << "::" << outbound_j << "\n";

        ex1 = exonOverlap(exon_map_, start, end) ;
//        cout << "KKK0:: " << ex1 << ":: " << start << "-" << end << "\n" ;

        if (nullptr != inbound_j && nullptr != outbound_j && inbound_j->get_intronic() && outbound_j->get_intronic()) {
            s1 << start ; s2 << end ;
            key = s1.str() + "-" + s2.str() ;
            if (exon_map_.count(key) > 0){
                ex1 = exon_map_[key] ;
                ex2 = ex1 ;
            } else {
                return ;
            }
        } else if (nullptr == ex1){
//            cout << "KK1\n" ;
            if ((end - start) <= MAX_DENOVO_DIFFERENCE){
                s1 << start ; s2 << end ;
                key = s1.str() + "-" + s2.str() ;
                ex1 = new Exon(start, end, in_db) ;
                ex2 = ex1 ;
                exon_map_[key] = ex1 ;

            } else {
                ex1 = addExon(exon_map_, start, EMPTY_COORD, in_db) ;
                ex2 = addExon(exon_map_, EMPTY_COORD, end, in_db) ;
            }

        } else {

//            cout << "KK2\n" ;

            ex2 = ex1;
            if (start < (ex1->get_start() - MAX_DENOVO_DIFFERENCE)){
                ex1 = addExon(exon_map_, start, EMPTY_COORD, in_db) ;
            } else if (start < ex1->get_start()){
                ex1->set_start(start) ;
            }

            if (end > (ex2->get_end() + MAX_DENOVO_DIFFERENCE)){
                ex2 = addExon(exon_map_, EMPTY_COORD, end, in_db) ;
            } else if (end > ex2->get_end()){
                ex2->set_end(end) ;
            }
        }
        if (nullptr != inbound_j){
            (ex1->ib).insert(inbound_j) ;
            inbound_j->set_acceptor(ex1) ;
        }
        if (outbound_j != nullptr){
            (ex2->ob).insert(outbound_j) ;
            outbound_j->set_donor(ex2) ;
        }

        return ;
    }


    void Gene::detect_exons(){
        vector<Ssite> ss_vec ;
        vector<Junction *> opened_exon ;
        Junction * last_5prime = nullptr ;
        Junction * first_3prime = nullptr ;

        for (const auto &jnc : junc_map_){
            if (!(jnc.second)->get_denovo_bl()) continue ;
            if ((jnc.second)->get_start() > 0) {
                Ssite s = {(jnc.second)->get_start(), 1, jnc.second} ;
                ss_vec.push_back(s) ;
            }
            if ((jnc.second)->get_end() > 0) {
                Ssite s = {(jnc.second)->get_end(), 0, jnc.second} ;
                ss_vec.push_back(s) ;
            }

        }
        sort(ss_vec.begin(), ss_vec.end(), sort_ss) ;

//        for(const auto & p : ss_vec){
//            cout<< "##" << p.coord << ", " << p.donor_ss << ", " << p.j << "\n" ;
//        }

        for(const auto & ss : ss_vec){
            if (ss.donor_ss) {
                if (opened_exon.size() > 0){
                    newExonDefinition(opened_exon.back()->get_end(), ss.coord, opened_exon.back(), ss.j, false) ;
                    opened_exon.pop_back() ;

                } else if ( 0 == opened_exon.size()) {
                    if (nullptr == first_3prime){
                        newExonDefinition((ss.j)->get_start()-10, (ss.j)->get_start(), nullptr, ss.j, false) ;
                    }else{
                        newExonDefinition(first_3prime->get_end(), ss.coord, first_3prime, ss.j, false) ;
                    }
                }
                last_5prime = ss.j ;
            }else{
                if (opened_exon.size() > 0){
                    if (nullptr != last_5prime){
                        for (const auto &jj2: opened_exon){
                            newExonDefinition(jj2->get_end(), last_5prime->get_start(), jj2, last_5prime, false) ;
                        }
                        last_5prime = nullptr ;
                        opened_exon.clear() ;
                        first_3prime = (ss.j) ;
                    }
                } else {
                    last_5prime = nullptr ;
                    first_3prime = (ss.j) ;
                }
                opened_exon.push_back(ss.j) ;
            }
        }

        for (const auto &jj2: opened_exon){
            newExonDefinition(jj2->get_end(), jj2->get_end()+10, jj2, nullptr, false) ;
        }

        ss_vec.clear() ;
        ss_vec.shrink_to_fit();
        return ;
    }

    bool is_lsv(set<Junction*> &juncSet, bool ss){
        unsigned int c1 = 0 ;
        unsigned int c2 = 0 ;
        for(const auto &juncIt:  juncSet){
            const int coord = ss? juncIt->get_end() : juncIt->get_start() ;
            c2 += (FIRST_LAST_JUNC != coord) ? 1:0 ;
            if (FIRST_LAST_JUNC != coord && juncIt->get_bld_fltr()) {
                ++c1 ;
            }
        }

        return (c2>1 and c1>0);
    }

    int detect_lsvs(vector<LSV*> &lsv_list, Gene* gObj){

        map<string, Exon*>::iterator exon_mapIt ;
        vector<LSV*> lsvGenes ;
        set<pair<set<string>, LSV*>> source ;
        LSV * lsvObj ;
        set<string> remLsv ;
//cout << "DETECT LSVS PER GENE1 "<< gObj<<"\n" ;
        for(const auto &exon_mapIt: gObj->exon_map_){
//        for (exon_mapIt = gObj->exon_map_.begin(); exon_mapIt != gObj->exon_map_.end(); ++exon_mapIt) {
//            cout << "DETECT LSVS PER GENE2\n" ;
            Exon * ex = exon_mapIt.second ;
//            cout << "DETECT LSVS PER GENE2.1: " << ex <<"\n" ;
            if (is_lsv(ex->ob, true)) {
//cout << "DETECT LSVS PER GENE2.2: " << ex <<"\n" ;
                lsvObj = new LSV(gObj, ex, true) ;
//cout << "DETECT LSVS PER GENE2.5 "<< gObj<<"\n" ;
                set<string> t1 ;
//cout << "DETECT LSVS PER GENE2.6 "<< gObj<<"\n" ;
                for (const auto &jl1: lsvObj->get_junctions()){
                    t1.insert(jl1->get_key()) ;
                }
//cout << "DETECT LSVS PER GENE2.7 "<< gObj<<"\n" ;
                pair<set<string>, LSV*> _p1 (t1, lsvObj) ;
//cout << "DETECT LSVS PER GENE2.8 "<< gObj<<"\n" ;
                source.insert(_p1) ;
//cout << "DETECT LSVS PER GENE2.9 "<< gObj<<"\n" ;
                lsvGenes.push_back(lsvObj) ;
//cout << "DETECT LSVS PER GENE2.95 "<< gObj<<"\n" ;
            }
//cout << "DETECT LSVS PER GENE3 "<< ex->get_start()<< " :: "<<ex->get_end()<<"\n" ;
            if (is_lsv(ex->ib, false)) {
//cout << "DETECT LSVS PER GENE3.1\n" ;
                lsvObj = new LSV(gObj, ex, false) ;
//cout << "DETECT LSVS PER GENE3.5\n" ;
                set<string> t1 ;
//cout << "DETECT LSVS PER GENE3.6\n" ;
                for (const auto &jl1: lsvObj->get_junctions()){
                    t1.insert(jl1->get_key()) ;
                }
//cout << "DETECT LSVS PER GENE3.7\n" ;

                lsvGenes.push_back(lsvObj) ;
//                for (const auto &slvs: source){
//cout << "DETECT LSVS PER GENE6\n" ;
//                    set<string> d1 ;
//                    set<string> d2 ;
//
//                    set_difference((slvs.first).begin(), (slvs.first).end(), t1.begin(), t1.end(), inserter(d1, d1.begin())) ;
//                    set_difference(t1.begin(), t1.end(), (slvs.first).begin(), (slvs.first).end(), inserter(d2, d2.begin())) ;
//cout << "DETECT LSVS PER GEN75\n" ;
//                    if (d1.size()>0 and d2.size()>0){
//                        lsvGenes.push_back(lsvObj) ;
//                    } else if (d1.size() >0 or (d1.size()==0 and d2.size() == 0)){
//cout << "DELETE 1\n" ;
//                        delete lsvObj ;
//                        break ;
//
//                    } else if (d2.size() >0){
//                        lsvGenes.push_back(lsvObj) ;
//                        remLsv.insert(slvs.second->get_id()) ;
//                    }
//                }
            }
//            cout << "DETECT LSVS PER GENE3.1b\n" ;
        }

//        for(const auto &l: lsvGenes){
//            cout<< "##" << l->get_id() << "\n" ;
//        }
//cout << "DETECT LSVS PER GENE9\n" ;

        #pragma omp critical
        for(const auto &l: lsvGenes){
//cout << "DETECT LSVS PER GENE9.5"<< l->get_id() <<"\n" ;
            if (remLsv.count(l->get_id())==0){
//                cout << "DETECT LSVS PER GENE9.7"<< l->get_id() <<"\n" ;
                lsv_list.push_back(l) ;
            }else{
//                cout << "DETECT LSVS PER GENE9.8"<< l->get_id() <<"\n" ;
                delete l;
            }
        }
//cout << "DETECT LSVS PER GENE10\n" ;
        return  0;

    }

    void sortGeneList(vector<Gene*> &glist) {
        sort(glist.begin(), glist.end(),[](Gene* a, Gene* b) {
                                return (a->get_start() < b->get_start()) || (a->get_start() == b->get_start() && a->get_end() < b->get_end()) ;}) ;
    }

    bool LSV::gather_lsv_info(float* source, float* target, list<Jinfo*> &info, map<string, Jinfo> &tlb,
                                                                                  unsigned int msample){


//        target = calloc((junctions_.size*msample), sizeof(float)));

        float * s1 ;
        float * t1 = target ;
        unsigned int count = 0 ;
        cout<< id_ << " source @: " << source<< " target @: " << target<<" ##\n" ;
        for(const auto &j: junctions_){
            const string key = gObj_->get_chromosome() + ":" + to_string(j->get_start()) + "-" + to_string(j->get_end()) ;
            if(tlb.count(key)>0){
                const int idx = tlb[key].index ;
                s1 = source + (msample*idx) ;
                for(unsigned int m=0; m<msample; m++){
                    t1 = s1 ;
                    t1 ++ ;
                    s1 ++ ;
                }
                count ++ ;
                info.push_back(&tlb[key]) ;
            } else {
                t1 += msample ;

            }
        }
        return count>0 ;
    }

    bool positive(lsvtype a, lsvtype b){
        return ( a.coord<b.coord ||
                (a.coord == b.coord && a.ref_coord<b.ref_coord) ||
                (a.ex_ptr->get_intron() && !b.ex_ptr->get_intron()) ) ;
    }

    bool reverse(lsvtype a, lsvtype b){
        return ( a.coord>b.coord ||
                (a.coord == b.coord && a.ref_coord>b.ref_coord) ||
                (a.ex_ptr->get_intron() && !b.ex_ptr->get_intron()) ) ;
    }

    string LSV::set_type(Exon* ref_ex, bool ss){

        string ref_exon_id = to_string(ref_ex->get_start()) + "-" + to_string(ref_ex->get_end()) ;
        vector<lsvtype> sp_list ;
        set<unsigned int> ref_ss_set ;
        set<Junction *> junc_list = ss? ref_ex->ob: ref_ex->ib ;

        for (const auto &j: junc_list){

            Exon * ex = ss? j->get_acceptor() :  j->get_donor() ;
            if (ex==nullptr) continue ;
            const int coord = ss? j->get_end() : j->get_start() ;
            const int ref_coord = ss ? j->get_start() : j->get_end() ;
            lsvtype lsvtypeobj = {coord, ref_coord, ex, j} ;
            sp_list.push_back(lsvtypeobj) ;
        }

        bool b = (gObj_->get_strand() == '+') ;
        if (b) sort(sp_list.begin(), sp_list.end(), positive) ;
        else sort(sp_list.begin(), sp_list.end(), reverse) ;

        string ext_type = (ss != b) ? "t" : "s" ;

//        string prev_ex ;
//        for (const auto &ptr: sp_list){
//            if(ptr.ex_ptr == nullptr) continue ;
//            prev_ex = to_string(ptr.ex_ptr->get_start()) + "-" + to_string(ptr.ex_ptr->get_end()) ;
//            break ;
//        }

        string prev_ex = to_string(sp_list[0].ex_ptr->get_start()) + "-" + to_string(sp_list[0].ex_ptr->get_end()) ;
        unsigned int excount = 1 ;
        if (ss) for (const auto &j: ref_ex->ob) ref_ss_set.insert(j->get_start()) ;
        else for (const auto &j: ref_ex->ib) ref_ss_set.insert(j->get_end()) ;

        unsigned int jidx = 0 ;

//cout << " ## " << gObj_->get_id() << "strand:"<< gObj_->get_strand()<<" SS: "<< ss << " and " << b << "\n" ;


        for (const auto &ptr: sp_list){
            jidx ++ ;
            const string exid = to_string((ptr.ex_ptr)->get_start()) + "-" + to_string((ptr.ex_ptr)->get_end()) ;
            if (exid == ref_exon_id) continue ;
            if ((ptr.jun_ptr)->get_intronic()){
                ext_type = ext_type + "|i" ;
                junctions_.push_back(ptr.jun_ptr) ;
                continue ;
            }
            unsigned int total = 0 ;
            unsigned int pos = 0 ;
            if (ss) {
                 set<unsigned int> ss_set ;
                 for (const auto &j: (ptr.ex_ptr)->ib){
                    ss_set.insert(j->get_end()) ;
                 }
                 total = ss_set.size() ;
                 pos = distance(ss_set.begin(), ss_set.find((ptr.jun_ptr)->get_end()))+1 ;
            }else{
                 set<unsigned int> ss_set ;
                 for (const auto &j: (ptr.ex_ptr)->ob){
                    ss_set.insert(j->get_start()) ;
                 }
                 total = ss_set.size() ;
                 pos = distance(ss_set.begin(), ss_set.find((ptr.jun_ptr)->get_start()))+1 ;

            }
            if(prev_ex != exid){
                prev_ex = exid ;
                excount += 1 ;
            }

            ext_type = ext_type + "|" + to_string(jidx) + "e" + to_string(excount) + "."
                                + to_string(pos) + "o" +  to_string(total) ;
            junctions_.push_back(ptr.jun_ptr) ;
        }
        return ext_type ;
    }


    void Gene::update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                  unsigned int denovo_thresh, unsigned int min_experiments){
        for(const auto &p: junc_map_){
            (p.second)->update_flags(efflen, minreads, minpos, denovo_thresh, min_experiments) ;
            (p.second)->clear_nreads(is_last_exp) ;
        }
        return ;
    }


    void Gene::print_exons(){

        for(const auto & ex: exon_map_ ){
            cout << "EXON:: "<< ex.first << "\n" ;
            for (const auto & j1: (ex.second)->ib){
                cout << "<<< " << j1->get_start() << "-" << j1->get_end() << "\n";
            }
            for (const auto & j1: (ex.second)->ob){
                cout << ">>>" << j1->get_start() << "-" << j1->get_end() << "\n";
            }
        }
    }
}


