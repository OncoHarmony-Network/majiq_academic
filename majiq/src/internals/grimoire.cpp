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

#include <omp.h>

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

/*
    Gene functions
*/

    bool islowergene(Gene * a, Gene * b){
        return (a->get_start() < b->get_start()) || (a->get_start() == b->get_start() && a->get_end() < b->get_end());
    }

    void sortGeneList(vector<Gene*> &glist) {
        sort(glist.begin(), glist.end(), islowergene) ;
        return ;
    }

    bool islowerintron(Intron * a, Intron * b){
        return (a->get_start() < b->get_start()) || (a->get_start() == b->get_start() && a->get_end() < b->get_end());
    }

    void Gene::newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db){

        Exon *ex1 = nullptr, *ex2 = nullptr ;
        string key ;
        stringstream s1 ;
        stringstream s2 ;
        if ((end - start) < 1) return ;

        ex1 = exonOverlap(exon_map_, start, end) ;
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

    void Gene::fill_junc_tlb(map<string, vector<string>> &tlb){

        for(const auto &j: junc_map_){
            if (!(j.second)->get_denovo_bl()) continue ;
            const string key = chromosome_ + ":" + strand_ + ":" + j.first ;
            const string key2 = chromosome_ + ":.:" + j.first ;
            if(tlb.count(key) == 0){
                tlb[key] = vector<string>() ;
            }
            if(tlb.count(key2) == 0){
                tlb[key2] = vector<string>() ;
            }
            tlb[key].push_back(id_) ;
            tlb[key2].push_back(id_) ;
        }
        return ;
    }

    void Gene::update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                  unsigned int denovo_thresh, unsigned int min_experiments){
        for(const auto &p: junc_map_){
            (p.second)->update_flags(efflen, minreads, minpos, denovo_thresh, min_experiments) ;
            (p.second)->clear_nreads(is_last_exp) ;
        }
        return ;
    }

    void Gene::detect_introns(vector<Intron*> intronlist){

//int tid = omp_get_thread_num();
//cout << tid<< " GEN ss_vec for gene:"<< id_ << "\n" ;
        vector<Ssite> ss_vec ;
        vector<Junction *> opened_exon ;
        int start_ir = 0 ;
        int end_ir = 0 ;



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
//cout << tid<< " Analize ss_vec for gene:"<< id_ << "\n" ;
        for(const auto & ss : ss_vec){
            if (ss.donor_ss) {
                start_ir = ss.coord ;

            } else {
                if (start_ir <= 0) {
                    continue ;
                } else {
                    end_ir = ss.coord ;
                    // new intron
                    // who creates the intron int_ptrs ?
                    // check if the introns are created in the io_bam functions for bam
                    #pragma omp critical
                    {
                        Intron * irObj = new Intron(start_ir, end_ir, false, this) ;
                        intronlist.push_back(irObj) ;
                    }

                    start_ir = 0 ;
                }
            }
        }
//cout << tid<< " END detect_introns for  gene:"<< id_ << "\n" ;
    }

    void Gene::add_intron(Intron * inIR_ptr, unsigned int min_exps){
        bool found = false ;

        for (const auto &ir: intron_vec_){
            if (ir->get_end() < inIR_ptr->get_start() || ir->get_start() > inIR_ptr->get_end()) continue ;
            if (ir->get_end() >= inIR_ptr->get_start() && ir->get_start() <= inIR_ptr->get_end()){
                ir->overlaping_intron(inIR_ptr) ;
                ir->update_flags(min_exps) ;
                found = true ;
            }
        }

        if(!found){
            inIR_ptr->update_flags(min_exps) ;
            intron_vec_.push_back(inIR_ptr) ;
        }else{
            delete(inIR_ptr) ;
        }
    }

    void Gene::connect_introns(){

        sort(intron_vec_.begin(), intron_vec_.end(), islowerintron) ;

        const int n_itrons = intron_vec_.size() ;
        int intron_index = 0 ;
        Exon * prev_ex = nullptr ;
        int ir_start = 0 ;

        for(const auto & ex: exon_map_ ){
            if (ir_start != 0) {
                const int ir_end = (ex.second)->get_start() - 1 ;
    cout << "CONNECT INTRONS: prev_ex: " << prev_ex << " ir_coord: "<< ir_start << " :: " << ir_end << "\n" ;
                while(intron_index< n_itrons){
                    Intron * ir_ptr = intron_vec_[intron_index];

                    if (ir_start <= ir_ptr->get_end() && ir_end >= ir_ptr->get_start() && ir_ptr->get_ir_flag()){
                        prev_ex->ob_irptr = ir_ptr ;
                        (ex.second)->ib_irptr = ir_ptr ;
                    }
                    ++intron_index ;
                }
            }
            ir_start = (ex.second)->get_end() + 1;
            prev_ex = (ex.second) ;
        }

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

    string Intron::get_key(Gene * gObj) {
        return(gObj->get_chromosome() + ":" + to_string(start_) + "-" + to_string(end_)) ;
    }

    bool Intron::is_reliable(float min_intron_cov){
        return true ;
    }

/*
    LSV fuctions
*/

    bool is_lsv(Exon * ex, bool ss){
//cout << "in is_LSVS\n" ;
        unsigned int c1 = 0 ;
        unsigned int c2 = 0 ;
        set<Junction*> &juncSet = ss? ex->ob : ex->ib ;
        for(const auto &juncIt:  juncSet){
            const int coord = ss? juncIt->get_end() : juncIt->get_start() ;
            c2 += (FIRST_LAST_JUNC != coord) ? 1:0 ;
            if (FIRST_LAST_JUNC != coord && juncIt->get_bld_fltr()) {
                ++c1 ;
            }
        }
//cout << "intron is_LSVS: "<< ss<<"\n" ;
        Intron * ir_ptr = ss? ex->ob_irptr : ex->ib_irptr ;
//cout << "intron is_LSVS ptr: "<< ir_ptr<<"\n" ;
        if (ir_ptr != nullptr){
//cout << " is_LSVS in if\n" ;
            ++c2 ;
            const int c = ir_ptr->get_ir_flag()? 1 : 0 ;
            c1 = c1 + c ;
        }
//cout << "out is_LSVS\n" ;
        return (c2>1 and c1>0);
    }

    inline void LSV::get_variations(set<string> &t1){
        for (const auto &jl1: junctions_){
            t1.insert(jl1->get_key()) ;
        }
        if(ir_ptr_ != nullptr){
            t1.insert(ir_ptr_->get_key()) ;
        }
    }

    int detect_lsvs(vector<LSV*> &lsv_list, Gene* gObj){

        map<string, Exon*>::iterator exon_mapIt ;
        vector<LSV*> lsvGenes ;
        set<pair<set<string>, LSV*>> source ;
        LSV * lsvObj ;
        set<string> remLsv ;

        for(const auto &exon_mapIt: gObj->exon_map_){
            Exon * ex = exon_mapIt.second ;
//            cout<< "IN 0 :: \n" ;
            if (is_lsv(ex, true)) {
                lsvObj = new LSV(gObj, ex, true) ;
                set<string> t1 ;
//cout<< "IN 0a :: \n" ;
                lsvObj->get_variations(t1) ;
//cout<< "IN 1a :: \n" ;
                pair<set<string>, LSV*> _p1 (t1, lsvObj) ;
                source.insert(_p1) ;
                lsvGenes.push_back(lsvObj) ;
            }
            if (is_lsv(ex, false)) {
                lsvObj = new LSV(gObj, ex, false) ;
                set<string> t1 ;
//cout<< "IN 0b :: \n" ;
                lsvObj->get_variations(t1) ;
//cout<< "IN 1b :: \n" ;
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
//cout << "NEXT EXON::" << ex->start_ << "-" << ex->end_ << "\n" ;
        }
//cout<< "IN 2 :: \n" ;
        int nlsv = 0 ;
        #pragma omp critical
        for(const auto &l: lsvGenes){
            if (remLsv.count(l->get_id())==0){
                ++nlsv ;
                lsv_list.push_back(l) ;
            }else{
                delete l;
            }
        }
//cout<< "OUT LSV\n" ;
        return  nlsv;

    }

    bool positive(lsvtype a, lsvtype b){
        return ( a.coord<b.coord ||
                (a.coord == b.coord && a.ref_coord<b.ref_coord)) ;
    }

    bool reverse(lsvtype a, lsvtype b){
        return ( a.coord>b.coord ||
                (a.coord == b.coord && a.ref_coord>b.ref_coord)) ;
    }



    bool LSV::gather_lsv_info(float* source, float* target, list<Jinfo*> &info, map<string, Jinfo> &tlb,
                                                                                  unsigned int msample){
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

    string LSV::set_type(Exon* ref_ex, bool ss){
        string ref_exon_id = to_string(ref_ex->get_start()) + "-" + to_string(ref_ex->get_end()) ;
        vector<lsvtype> sp_list ;
        set<unsigned int> ref_ss_set ;
        set<Junction *> junc_list = ss? ref_ex->ob: ref_ex->ib ;

        for (const auto &j: junc_list){
            Exon * ex = ss? j->get_acceptor() :  j->get_donor() ;
            const int coord = ss? j->get_end() : j->get_start() ;
            const int ref_coord = ss ? j->get_start() : j->get_end() ;
            if (ex==nullptr || coord < 0 || ref_coord < 0) continue ;
            lsvtype lsvtypeobj = {coord, ref_coord, ex, j} ;
            sp_list.push_back(lsvtypeobj) ;
        }

        bool b = (gObj_->get_strand() == '+') ;
        if (b) sort(sp_list.begin(), sp_list.end(), positive) ;
        else sort(sp_list.begin(), sp_list.end(), reverse) ;
        string ext_type = (ss != b) ? "t" : "s" ;

        string prev_ex = to_string(sp_list[0].ex_ptr->get_start()) + "-" + to_string(sp_list[0].ex_ptr->get_end()) ;
        unsigned int excount = 1 ;
        if (ss) for (const auto &j: ref_ex->ob) ref_ss_set.insert(j->get_start()) ;
        else for (const auto &j: ref_ex->ib) ref_ss_set.insert(j->get_end()) ;

        unsigned int jidx = 0 ;

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
//                 cout << "LSV_TYPE 4.5.1 "<<ptr.ex_ptr<<"\n";
                 for (const auto &j: (ptr.ex_ptr)->ob){
                    ss_set.insert(j->get_start()) ;
                 }
                 total = ss_set.size() ;
                 pos = distance(ss_set.begin(), ss_set.find((ptr.jun_ptr)->get_start()))+1 ;
//                 cout << "LSV_TYPE 4.5.5\n";

            }
//cout << "LSV_TYPE 4.6\n";
            if(prev_ex != exid){
                prev_ex = exid ;
                excount += 1 ;
            }
//cout << "LSV_TYPE 4.7\n";
            ext_type = ext_type + "|" + to_string(jidx) + "e" + to_string(excount) + "."
                                + to_string(pos) + "o" +  to_string(total) ;
//cout << "LSV_TYPE 4.8\n";
            junctions_.push_back(ptr.jun_ptr) ;
//cout << "LSV_TYPE 4.9\n";
        }

        if (ir_ptr_ != nullptr) ext_type += "|i" ;
        return ext_type ;
    }

    int GeneSearch(vector<Gene *> & a, int n, int coord) {
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

    int intronSearch(vector<Intron *> & a, int n, const int start_pos) {
        int l = 0 ;
        int h = n ; // Not n - 1
        while (l < h) {
            int mid = (l + h) / 2 ;
            if (a[mid]->get_start()>=start_pos) {
                h = mid ;
            } else {
                l = mid +1 ;
            }
        }
        return l-1;
    }


    vector<Intron *>& find_intron_retention(vector<Gene*> & gene_list, int start, int end){
        vector<Intron*> ir_vec ;
        const int n = gene_list.size() ;
        const int gIdx = GeneSearch(gene_list, n, end) ;
        if (gIdx <0){ return ir_vec;}
        if (gene_list[gIdx]->get_start()<= end && gene_list[gIdx]->get_end()>= start){
            const int nir = (gene_list[gIdx]->intron_vec_).size() ;
            int irIdx = intronSearch(gene_list[gIdx]->intron_vec_, nir, start);
            for(int i = irIdx; i<nir; ++i){
                Intron * irp = gene_list[gIdx]->intron_vec_[i] ;
                if(irp->get_start()> end){
                    break ;
                }
                if(irp->get_end() < start){
                    continue ;
                }
                ir_vec.push_back(irp) ;
            }
        }
        return ir_vec;
    }

}


