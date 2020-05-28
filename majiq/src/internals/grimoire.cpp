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
#include <functional>
#include "grimoire.hpp"
#include "io_bam.hpp"

#include <omp.h>

using namespace std;

namespace grimoire {

    bool positive(lsvtype a, lsvtype b){
        return ( a.ref_coord<b.ref_coord ||
                (a.ref_coord == b.ref_coord && a.coord<b.coord)) ;
    }

    bool reverse(lsvtype a, lsvtype b){
        return ( a.ref_coord>b.ref_coord ||
                (a.ref_coord == b.ref_coord && a.coord>b.coord)) ;
    }

    bool fless(unsigned int lhs, unsigned int rhs) { return lhs<rhs ; }
    bool fgrt(unsigned int lhs, unsigned int rhs) { return lhs>rhs ; }


    bool sort_ss(const Ssite &a, const Ssite &b){

        bool crd =  (a.coord<b.coord) ;
        bool r = (a.coord == b.coord) && (a.donor_ss < b.donor_ss);  // acceptors before donors if given same coordinate
        bool same = (a.coord == b.coord) && (a.donor_ss == b.donor_ss) && ((a.j->length()) < b.j->length()) ;
        return  crd || r || same ;
    }

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
                && ( start <= x->get_end() ) && ( end >= x->get_start()) ) {
                return x ;
            }
        }
        return nullptr;
    }


    string  Junction::get_key(Gene * gObj, int strandness) {

        string strand = (UNSTRANDED == strandness) ? "." : to_string(gObj->get_strand()) ;
        return(gObj->get_chromosome() + ":" + strand + ":" + to_string(start_) + "-" + to_string(end_)) ;
    }

    string  Junction::get_key(Gene * gObj) {
        return(gObj->get_id() + ":" + to_string(start_) + "-" + to_string(end_)) ;
    }

/*
    Gene functions
*/

    void sortGeneList(vector<Gene*> &glist) {
        sort(glist.begin(), glist.end(), Gene::islowerRegion<Gene>) ;
        return ;
    }

    void Gene::newExonDefinition(int start, int end, Junction *inbound_j, Junction *outbound_j, bool in_db){

        Exon *ex1 = nullptr, *ex2 = nullptr ;
        string key ;
        stringstream s1 ;
        stringstream s2 ;
        if ((end - start) < 0) return ;
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

                const int st1 = (inbound_j == nullptr) ? EMPTY_COORD : start ;
                const int nd1 = (outbound_j == nullptr) ? EMPTY_COORD : end ;
                ex1 = addExon(exon_map_, st1, nd1, in_db) ;
                ex2 = ex1 ;

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
            int coord = last_5prime == nullptr ? jj2->get_end()+10 : last_5prime->get_start() ; 
            newExonDefinition(jj2->get_end(), coord, jj2, last_5prime, false) ;
        }
        ss_vec.clear() ;
        ss_vec.shrink_to_fit();
        return ;
    }


    void Gene::update_junc_flags(int efflen, bool is_last_exp, unsigned int minreads, unsigned int minpos,
                                  unsigned int denovo_thresh, unsigned int min_experiments, bool denovo){
        for(const auto &p: junc_map_){
            (p.second)->gen_and_update_flags(efflen, minreads, minpos, denovo_thresh, min_experiments, denovo) ;
            (p.second)->clear_nreads(is_last_exp) ;
        }
        return ;
    }

    void Gene::updateFlagsFromJunc(string key, unsigned int sreads, unsigned int minreads_t, unsigned int npos,
                               unsigned int minpos_t, unsigned int denovo_t, bool denovo, int minexp, bool reset){
        if (junc_map_.count(key) > 0){
            omp_set_lock(&map_lck_) ;
            Junction * jnc = junc_map_[key] ;
            jnc->update_flags(sreads, minreads_t, npos, minpos_t, denovo_t, minexp, denovo) ;
            jnc->clear_nreads(reset) ;
            omp_unset_lock(&map_lck_) ;
        }
    }

    void Gene::initialize_junction(string key, int start, int end, shared_ptr<vector<float>> nreads_ptr, bool simpl){

        omp_set_lock(&map_lck_) ;
        if (junc_map_.count(key) == 0){
            junc_map_[key] = new Junction(start, end, false, simpl) ;
        }
        junc_map_[key]->set_nreads_ptr(nreads_ptr) ;
        omp_unset_lock(&map_lck_) ;
        return ;
    }


    void Gene::detect_introns(vector<Intron*> &intronlist, bool simpl){

        vector<Ssite> ss_vec ;
        vector<Junction *> opened_exon ;
        int start_ir = 0 ;
        int end_ir = 0 ;

        for (const auto &jnc : junc_map_){
            if (!(jnc.second)->get_denovo_bl()) continue ;


            if ((jnc.second)->get_start() > 0 && (jnc.second)->get_donor()->get_start()>0) {
                Ssite s = {(jnc.second)->get_start(), true, jnc.second} ;
                ss_vec.push_back(s) ;
            }
            if ((jnc.second)->get_end() > 0 && (jnc.second)->get_acceptor()->get_end()>0) {
                Ssite s = {(jnc.second)->get_end(), false, jnc.second} ;
                ss_vec.push_back(s) ;
            }
        }
        sort(ss_vec.begin(), ss_vec.end(), sort_ss) ;

        for(const auto & ss : ss_vec){
            if (ss.donor_ss) {
                start_ir = ss.coord ;
            } else {
                if (start_ir <= 0) {
                    continue ;
                } else {
                    // are we going to close a valid intron?
                    if ((ss.coord - start_ir) >= 2) {
                        // the intron has positive length, so yes!
                        end_ir = ss.coord;
                        // create the intron
                        Intron * irObj = new Intron(
                            start_ir + 1, end_ir - 1, false, this, simpl
                        );
                        // critical region to add the intron to intronlist
                        #pragma omp critical
                        {
                            intronlist.push_back(irObj);
                        }
                    }
                    // either way, close open intron
                    start_ir = 0;
                }
            }
        }
    }

    void Gene::add_intron(Intron * inIR_ptr, float min_coverage, unsigned int min_exps, float min_bins, bool reset,
                          bool denovo){
        bool found = false ;
        for (const auto &ir: intron_vec_){
            if (ir->get_end() < inIR_ptr->get_start() || ir->get_start() > inIR_ptr->get_end()) continue ;
            if (ir->get_end() >= inIR_ptr->get_start() && ir->get_start() <= inIR_ptr->get_end()){
                ir->overlaping_intron(inIR_ptr) ;
//cerr << "found " << ir->get_key() << endl ;
                ir->update_flags(min_coverage, min_exps, min_bins) ;
                ir->clear_nreads(reset) ;
                found = true ;
            }
        }
        if(!found && denovo){
//cerr << "NOT found " << inIR_ptr->get_key() << endl ;

            inIR_ptr->update_flags(min_coverage, min_exps, min_bins) ;
            inIR_ptr->clear_nreads(reset) ;
            intron_vec_.push_back(inIR_ptr) ;
        }else{
            delete(inIR_ptr) ;
        }
    }

    void Gene::connect_introns(){

        sort(intron_vec_.begin(), intron_vec_.end(), Intron::islowerRegion<Intron>) ;

        const int n_itrons = intron_vec_.size() ;
        int intron_index = 0 ;
        Exon * prev_ex = nullptr ;
        int ir_start = 0 ;

        for(const auto & ex: exon_map_ ){

            if (((ex.second)->get_start() < 0) || (ex.second)->get_end()<0)
                continue ;
            const int ex_start = ((ex.second)->get_start() < 0) ? (ex.second)->get_end()-10 : (ex.second)->get_start() ;
            const int ex_end = ((ex.second)->get_end() < 0) ? (ex.second)->get_start()+10 : (ex.second)->get_end() ;

            if (ir_start == 0 ) {
                ir_start = ex_end + 1;
                prev_ex = (ex.second) ;
                continue ;
            }

            if( prev_ex->get_end()> 0 && (ex.second)->get_start()>0) {
                const int ir_end = ex_start - 1 ;

                while(intron_index< n_itrons){
                    Intron * ir_ptr = intron_vec_[intron_index];
                    if (ir_ptr->get_end() >= ex_end ) break ;
                    if (ir_start <= ir_ptr->get_end() && ir_end >= ir_ptr->get_start() && ir_ptr->get_ir_flag()){

                        if (prev_ex->ob_irptr != nullptr){
                            (prev_ex->ob_irptr)->unset_markd() ;
                        }
//    cerr << "#1 " << ir_start << "-" << ir_end<< " :: " << prev_ex->ob_irptr->get_gene() << ":" << prev_ex->ob_irptr->get_start() << "-" << prev_ex->ob_irptr->get_end()<< "\n" ;
//    cerr << "#2 " << ir_start << "-" << ir_end<< " :: " << ir_ptr->get_gene() << ":" << ir_ptr->get_start() << "-" << ir_ptr->get_end()<< "\n" ;
                        prev_ex->ob_irptr = ir_ptr ;
                        (ex.second)->ib_irptr = ir_ptr ;
                        ir_ptr->update_boundaries(prev_ex->get_end() +1, (ex.second)->get_start() -1) ;
                        ir_ptr->set_markd();
                    }

                    ++intron_index ;
                }
                ir_start = ex_end + 1;
                prev_ex = (ex.second) ;
            }

        }

    }

    int
    Gene::get_constitutive_junctions(vector<string>& v){
        vector<Exon*> ex_vector ;
        for(const auto &exon_mapIt: exon_map_){
            ex_vector.push_back(exon_mapIt.second) ;
        }
        sort(ex_vector.begin(), ex_vector.end(), Exon::islowerRegion<Exon>) ;

        for(const auto &ex: ex_vector){
            // first, we look at inbound connections to determine whether
            // inbound intron has constitutive acceptor (but only if we have a
            // valid inbound intron in the first place)
            if (ex->ib_irptr != nullptr && !(ex->ib_irptr)->get_simpl_fltr()) {
                bool has_junctions = false ;  // we haven't seen junctions yet
                for (const auto &juncIt: (ex->ib)) {
                    const int coord = juncIt->get_start() ;
                    if (FIRST_LAST_JUNC != coord && !juncIt->get_simpl_fltr() && !juncIt->is_exitron()) {
                        has_junctions = true ;  // we have seen a junction
                        break ;  // there's no point at looking for more junctions
                    }
                }
                // if we don't have any in-bound junctions, constitutive acceptor
                if ((ex->ib_irptr)->get_ir_flag() && !has_junctions) {
                    // if ir passed group filters and there were no junctions...
                    (ex->ib_irptr)->set_const_acceptor() ;
                }
            }

            // now, we look at outbound connections.
            // we want to determine whether we have constitutive donor intron
            // or constitutive junction
            unsigned int num_outbound_junctions = 0 ;  // counting valid junctions
            Junction * jnc = nullptr ;  // will be assigned an outbound junction
            for (const auto &juncIt: (ex->ob)) {
                const int coord = juncIt->get_end() ;
                if (FIRST_LAST_JUNC != coord && !juncIt->get_simpl_fltr() && !juncIt->is_exitron()) {
                    // this junction is valid for being part of an event
                    ++num_outbound_junctions ;
                    if (num_outbound_junctions > 1) {
                        // there are multiple junctions so, no constitutive
                        // donor intron or constitutive junction.
                        // There is no point in looping over more junctions
                        break ;
                    }
                    jnc = juncIt ;  // change from nullptr to last junction
                }
            }
            // if there are multiple junctions, go to next exon
            if (num_outbound_junctions > 1) {
                continue ;
            }
            // 3 scenarios -- constitutive junction(1)/intron(2) or 1 junction, 1 intron
            if (ex->ob_irptr != nullptr && !(ex->ob_irptr)->get_simpl_fltr()) {
                // this is a valid intron.
                if ((ex->ob_irptr)->get_ir_flag() && num_outbound_junctions < 1) {
                    // the intron passed build filters, and no outbound
                    // junctions, so the intron is a const_donor
                    (ex->ob_irptr)->set_const_donor() ;
                }  // otherwise, it has 2 valid connections
            } else {
                // there is no valid intron.
                // if there is a junction that passed previous conditions and
                // passed build filters, then it is is a constant donor
                if (jnc != nullptr && jnc->get_bld_fltr()) {
                    // is it also a constitutive acceptor?
                    Exon * accex = jnc->get_acceptor() ;
                    if(!accex->is_lsv(false)){
                        // it is a constitutive acceptor, so our junction is
                        // indeed constitutive
                        jnc->set_constitutive() ;
                        // get tab-delimited output for dump-constitutive
                        string str_ln = id_ + "\t" + chromosome_ + "\t" +
                                        to_string(jnc->get_start()) + "\t" + to_string(jnc->get_end()) + "\t" +
                                        to_string(ex->get_start()) + "\t" + to_string(ex->get_end()) + "\t" +
                                        to_string(accex->get_start()) + "\t" + to_string(accex->get_end()) ;
                        // push tab-delimited output to shared memory
                        #pragma omp critical
                            v.push_back(str_ln) ;
                    }
                }  // otherwise, no junctions/introns or 1 junction that didn't pass build filters
            }
            // done looking at all exons
        }
        return 0 ;
    }

    int Gene::detect_lsvs(vector<LSV*> &lsv_list, bool lsv_strict) {

        map<string, Exon*>::iterator exon_mapIt ;
        vector<LSV*> lsvGenes ;
        set<pair<set<string>, LSV*>> source ;
        LSV * lsvObj ;
        set<string> remLsv ;
        const bool ss = strand_ == '+' ;

        vector<Exon*> ex_vector ;
        for(const auto &exon_mapIt: exon_map_){
            ex_vector.push_back(exon_mapIt.second) ;
        }
        if (ss)
            sort(ex_vector.begin(), ex_vector.end(), Exon::islowerRegion<Exon>) ;
        else
            sort(ex_vector.begin(), ex_vector.end(), Exon::isgreaterRegion<Exon>) ;

        for(const auto &ex: ex_vector){

            if (ex->is_lsv(ss)) {
                lsvObj = new LSV(this, ex, ss) ;
                set<string> t1 ;
                lsvObj->get_variations(t1) ;
                pair<set<string>, LSV*> _p1 (t1, lsvObj) ;
                source.insert(_p1) ;
                lsvGenes.push_back(lsvObj) ;
            }

            if (ex->is_lsv(!ss)) {
                lsvObj = new LSV(this, ex, !ss) ;
                set<string> t1 ;
                lsvObj->get_variations(t1) ;
                bool rem_src = false ;

                for (const auto &slvs: source){
                    set<string> d1 ;  // fill with connections unique to source
                    set<string> d2 ;  // fill with connections unique to target
                    set_difference((slvs.first).begin(), (slvs.first).end(), t1.begin(), t1.end(), inserter(d1, d1.begin())) ;
                    set_difference(t1.begin(), t1.end(), (slvs.first).begin(), (slvs.first).end(), inserter(d2, d2.begin())) ;

                    // when do we remove target LSV from output?
                    if (d2.size() == 0 && (lsv_strict || d1.size() == 0)) {
                        rem_src = true ;  // remove target LSV
                        break ;
                    }

                    // when do we remove source LSV from output?
                    if (lsv_strict && (d1.size() == 0 && d2.size() > 0)) {
                        remLsv.insert(slvs.second->get_id()) ;
                    }
                }
                if (rem_src){
                    delete lsvObj ;
                }else{
                    lsvGenes.push_back(lsvObj) ;
                }
            }
        }
        int nlsv = 0 ;
        #pragma omp critical
        {
            for(const auto &l: lsvGenes){
                if (remLsv.count(l->get_id())==0){
                    ++nlsv ;
                    lsv_list.push_back(l) ;
                }else{
                    delete l ;
                }
            }
        }
        return  nlsv;

    }


    void Exon::simplify(map<string, int>& junc_tlb, float simpl_percent, Gene* gObj, int strandness,
                        int denovo_simpl, int db_simple, int ir_simpl, bool last, unsigned int min_experiments){
        float sumall = 0 ;
        {
            vector<float> sumj ;
            vector<_Simplifiable *> jnc_vec ;
            vector<int> thrshld_vect;
            unsigned int i = 0 ;
            for(const auto &juncIt: ib){
                if (!juncIt->get_denovo_bl())
                    continue ;
                string key = juncIt->get_key(gObj, strandness) ;
                jnc_vec.push_back(juncIt) ;
                float s = junc_tlb.count(key)>0 ? junc_tlb[key] : 0 ;
                sumj.push_back(s) ;
                int min_reads = juncIt->get_annot() ?  db_simple: denovo_simpl ;
                thrshld_vect.push_back(min_reads) ;
                sumall += s ;
                i++ ;
            }
            if (ib_irptr != nullptr && ib_irptr->get_ir_flag()){
                const int ir_strt =  ib_irptr->get_start() ;
                const int ir_end  =  ib_irptr->get_end() ;
                jnc_vec.push_back(ib_irptr) ;
                string key = key_format(gObj->get_id(), ir_strt, ir_end, true) ;
                float s = junc_tlb.count(key)>0 ? junc_tlb[key] : 0 ;
                thrshld_vect.push_back(ir_simpl) ;
                sumj.push_back(s) ;
                sumall += s ;
            }

            for(i=0; i<jnc_vec.size(); i++){
                float x = (sumall >0) ? sumj[i]/sumall : 0 ;
                jnc_vec[i]->set_simpl_fltr(x<simpl_percent || sumj[i]< thrshld_vect[i], true) ;
            }
        }
        sumall = 0 ;
        {
            vector<float> sumj ;
            vector<_Simplifiable *> jnc_vec ;
            vector<int> thrshld_vect;
            unsigned int i = 0 ;
            for(const auto &juncIt: ob){
                if (!juncIt->get_denovo_bl())
                    continue ;
                string key = juncIt->get_key(gObj, strandness) ;
                jnc_vec.push_back(juncIt) ;
                float s = junc_tlb.count(key)>0 ? junc_tlb[key] : 0 ;
                sumj.push_back(s) ;
                int min_reads = juncIt->get_annot() ?  db_simple: denovo_simpl ;
                thrshld_vect.push_back(min_reads) ;
                sumall += s ;
                i++ ;
            }
            if (ob_irptr != nullptr && ob_irptr->get_ir_flag()){
                const int ir_strt =  ob_irptr->get_start() ;
                const int ir_end  =  ob_irptr->get_end() ;
                jnc_vec.push_back(ob_irptr) ;
                string key = key_format(gObj->get_id(), ir_strt, ir_end, true) ;
                float s = junc_tlb.count(key)>0 ? junc_tlb[key] : 0 ;
                thrshld_vect.push_back(ir_simpl) ;
                sumj.push_back(s) ;
                sumall += s ;
            }

            for(i=0; i<jnc_vec.size(); i++){
                float x = (sumall >0) ? sumj[i]/sumall : 0 ;
                jnc_vec[i]->set_simpl_fltr(x<simpl_percent || sumj[i]< thrshld_vect[i], false) ;
            }
        }
    }


    void Gene::simplify(map<string, int>& junc_tlb, float simpl_percent, int strandness, int denovo_simpl,
                        int db_simple, int ir_simpl, bool last, unsigned int min_experiments){
        for(const auto &ex: exon_map_){
            (ex.second)->simplify(junc_tlb, simpl_percent, this, strandness, denovo_simpl, db_simple,
                                  ir_simpl, last, min_experiments) ;
        }
        if (last){
            for (const auto &j: junc_map_){
                (j.second)->update_simpl_flags(min_experiments) ;
            }
            for(const auto &ex: exon_map_){
                if ((ex.second)->ib_irptr != nullptr)
                    ((ex.second)->ib_irptr)->update_simpl_flags(min_experiments) ;
//                if ((ex.second)->ob_irptr != nullptr)
//                    ((ex.second)->ob_irptr)->update_simpl_flags(min_experiments, false) ;
            }
        }
    }


    void Gene::print_exons(){

        for(const auto & ex: exon_map_ ){
            cerr << "EXON:: "<< ex.first << "\n" ;
            for (const auto & j1: (ex.second)->ib){
                cerr << "<<< " << j1->get_start() << "-" << j1->get_end() << "\n" ;
            }
            for (const auto & j1: (ex.second)->ob){
                cerr << ">>>" << j1->get_start() << "-" << j1->get_end() << "\n" ;
            }
        }
    }

    string Intron::get_key(Gene * gObj) {
        return(gObj->get_chromosome() + ":" + gObj->get_strand() + ":" + to_string(start_) + "-" + to_string(end_)
               + ":" + gObj->get_id()) ;
    }

   void Intron::normalize_readrates() {
        // exit early if nothing to normalize
        if (read_rates_ptr_ == nullptr) {
            return;
        }
        // get read rates vector from shared pointer
        vector<float> &read_rates_ = *read_rates_ptr_;
        // normalize each position
        for (int i = 0; i < numbins_; ++i) {
            if (read_rates_[i] > 0) {
                // use the number of intronic positions assigned to the i-th
                // bin to normalize read rates
                const int positions_in_bin = i < nxbin_mod_ ? nxbin_ + 1 : nxbin_;
                read_rates_[i] = read_rates_[i] / positions_in_bin;
            }
        }
        return;
    }

    bool Intron::is_reliable(float min_bins, int eff_len) {
        // exit early if possible
        if (length() <= 0 || numbins_ <= 0 || read_rates_ptr_ == nullptr) {
            return false;
        }

        // no need to distribute coverage if intron previously known reliable
        if (ir_flag_) {
            return true;
        }
        // otherwise, need to look at percent positions with nonzero coverage

        // how many bins could a read with eff_len cover given that nxbin_ + 1
        // is the maximum number of positions per bin?
        const int ext = (int) (eff_len / (nxbin_ + 1));
        // initialize vector of read coverage where we share coverage with up
        // to ext adjacent bins
        vector<float> cov (numbins_);
        // get read rates vector from shared pointer
        vector<float> &read_rates_ = *read_rates_ptr_;
        // distribute read_rates per bin to coverage shared among adjacent bins
        for(int i = 0; i < numbins_; i++) {
            if (read_rates_[i] > 0) {
                // XXX this only distributes coverage to the right, which could
                // be problematic for extremely short introns
                for (int j = 0; j < ext && (i + j) < numbins_; j++) {
                    cov[i + j] += read_rates_[i];
                }
            }
        }
        // get the numer/percent positions covered
        float npos = 0;
        for (const auto &p: cov) {
            npos += (p > 0) ? 1 : 0;
        }
        const float pct_pos = (npos > 0) ? (npos / numbins_) : 0;
        return (pct_pos >= min_bins);
     }

/*
    LSV fuctions
*/

    bool Exon::is_lsv(bool ss){
        // True if we have multiple valid junctions/introns, with at least one
        // passing group filters
        unsigned int c1 = 0 ;  // count junctions/introns passing group filters
        unsigned int c2 = 0 ;  // count valid junctions/introns
        set<Junction*> &juncSet = ss? ob : ib ;
        for(const auto &juncIt:  juncSet){
            if (juncIt->is_exitron()){
                // we don't count exitrons towards junction count for LSV
                continue ;
            }
            const int coord = ss? juncIt->get_end() : juncIt->get_start() ;
            bool is_valid = FIRST_LAST_JUNC != coord && !juncIt->get_simpl_fltr() ;
            c2 += is_valid ? 1:0 ;
            if (is_valid && juncIt->get_bld_fltr()) {
                ++c1 ;
            }
        }
        Intron * ir_ptr = ss? ob_irptr : ib_irptr ;
        if (ir_ptr != nullptr && !ir_ptr->get_simpl_fltr()){
            ++c2 ;
            const int c = ir_ptr->get_ir_flag()? 1 : 0 ;
            c1 = c1 + c ;
        }
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

    bool LSV::gather_lsv_info(float* source, float* target, list<Jinfo*> &info, map<string, Jinfo> &tlb,
                                                                                  unsigned int msample){
        float * s1 ;
        float * t1 = target ;
        unsigned int count = 0 ;
//        cout<< id_ << " source @: " << source<< " target @: " << target<<" ##\n" ;
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
            if (ex==nullptr || coord < 0 || ref_coord < 0 || j->get_simpl_fltr()) continue ;
            lsvtype lsvtypeobj = {coord, ref_coord, ex, j} ;
            sp_list.push_back(lsvtypeobj) ;
        }

        bool b = (gObj_->get_strand() == '+') ;
        if (b) sort(sp_list.begin(), sp_list.end(), positive) ;
        else sort(sp_list.begin(), sp_list.end(), reverse) ;

        bool(*bfunc)(unsigned int, unsigned int) = b ? fless: fgrt ;

        string ext_type = (ss != b) ? "t" : "s" ;

        if (ss) for (const auto &j: ref_ex->ob) ref_ss_set.insert(j->get_start()) ;
        else for (const auto &j: ref_ex->ib) ref_ss_set.insert(j->get_end()) ;

        unsigned int jidx = 0 ;
        int prev_coord = 0 ;
        map<string, unsigned int > exDct;  // exon key -> exon count in strand order

        // fill exDct beforehand since it is not in same order as junctions
        set <pair<int, string>> other_exons;
        for (const auto &ptr: sp_list) {
            // don't add exon for exitrons
            if (ptr.ex_ptr == ref_ex) {
                continue;
            }
            // get exon id for exDct
            const string exid = to_string((ptr.ex_ptr)->get_start()) + "-" + to_string((ptr.ex_ptr)->get_end());
            // get coordinate for sorting
            // we use this coordinate because unique and handles half exons in order
            const int sort_coordinate = (b ? 1 : -1) * (ss ? ptr.ex_ptr->get_start() : ptr.ex_ptr->get_end());
            // should be unique pair, sorted by valid coordinate
            other_exons.insert(make_pair(sort_coordinate, exid));
        }
        for (const auto &ex_pair: other_exons) {
            const string exid = ex_pair.second;
            const unsigned int ex_ct = exDct.size() + 1;
            exDct[exid] = ex_ct;
        }

        // now loop over junctions
        for (const auto &ptr: sp_list){
            jidx = (prev_coord != ptr.ref_coord) ? jidx+1 : jidx ;
            prev_coord = ptr.ref_coord ;

            const string exid = to_string((ptr.ex_ptr)->get_start()) + "-" + to_string((ptr.ex_ptr)->get_end()) ;
            if (exid == ref_exon_id) continue ;

            unsigned int total = 0 ;
            unsigned int pos = 0 ;

            set<unsigned int, bool(*)(unsigned int,unsigned int)> ss_set (bfunc);

            if (ss) {
                 for (const auto &j: (ptr.ex_ptr)->ib){
                     if (
                         j->get_donor() != nullptr
                         && j->get_start() >= 0
                         && !(j->get_simpl_fltr())
                     ) {
                         // only count junctions that could make junc_list
                         // (i.e. junction starting/ending at exon(s), unsimplified)
                         ss_set.insert(j->get_end());
                     }
                 }
                 total = ss_set.size() ;
                 pos = distance(ss_set.begin(), ss_set.find((ptr.jun_ptr)->get_end()))+1 ;
            }else{
                 for (const auto &j: (ptr.ex_ptr)->ob){
                     if (
                         j->get_acceptor() != nullptr
                         && j->get_end() >= 0
                         && !(j->get_simpl_fltr())
                     ) {
                         // only count junctions that could make junc_list
                         // (i.e. junction starting/ending at exon(s), unsimplified)
                         ss_set.insert(j->get_start()) ;
                     }
                 }
                 total = ss_set.size() ;
                 pos = distance(ss_set.begin(), ss_set.find((ptr.jun_ptr)->get_start()))+1 ;

            }

            try {
                ext_type = (
                    ext_type
                    + "|" + to_string(jidx) + "e" + to_string(exDct.at(exid))
                    + "." + to_string(pos) + "o" + to_string(total)
                );
            } catch (const std::out_of_range& e) {
                // exid wasn't in exDct
                // this key should have been added during its construction...
                cerr << "Assertion error: exon (exid = '" << exid
                    << "') not found in exDct in LSV::set_type()\n";
                throw(std::out_of_range("Missing key " + exid + " in exDct"));
            }
            junctions_.push_back(ptr.jun_ptr) ;
        }

        if (ext_type.length() > MAX_TYPE_LENGTH){
            ext_type = (ss != b) ? "t" : "s" ;
            ext_type += "|na"  ;
        }
        if (ir_ptr_ != nullptr) ext_type += "|i" ;
        return ext_type ;
    }

    void fill_junc_tlb(vector<LSV*>& lsv_list, map<string, int>& tlb){

        for (const auto &l: lsv_list){
            const string gid = (l->get_gene())->get_id() ;
            for (const auto &j: l->get_junctions()){
                const string k = j->get_key(l->get_gene()) ;
                const int n  = tlb.size() ;
                if (tlb.count(k) == 0)
                    tlb[k] = n ;
            }

            Intron * ir_ptr = l->get_intron() ;
            if (ir_ptr != 0){
                const string k = "IR:" + gid + ":" + to_string(ir_ptr->get_start()) + "-" + to_string(ir_ptr->get_end()) ;
                const int n  = tlb.size() ;
                if (tlb.count(k) == 0)
                    tlb[k] = n ;
            }
        }
    }

    vector<Intron *> find_intron_retention(Gene * gObj, int start, int end){
        vector<Intron*> ir_vec ;
        vector<Intron *>::iterator low = lower_bound (gObj->intron_vec_.begin(), gObj->intron_vec_.end(),
                                                      start, _Region::func_comp ) ;
        if (low ==  gObj->intron_vec_.end()) return ir_vec ;
        for (; low != gObj->intron_vec_.end() ; low++){
            Intron * irp = *low;

            if(irp->get_start()> end){
                break ;
            }
            if(irp->get_end() < start){
                continue ;
            }
            if (irp->get_ir_flag() && irp->is_connected())
                ir_vec.push_back(irp) ;
        }
        return ir_vec ;
    }

    void find_gene_from_junc(map<string, vector<overGene*>> & glist, string chrom, char strand, int start, int end,
                             vector<Gene*>& oGeneList, bool ir, bool simpl){

        Junction * junc = new Junction(start, end, false, simpl) ;
        const string key = junc->get_key() ;
        vector<overGene*>::iterator low = lower_bound (glist[chrom].begin(), glist[chrom].end(),
                                                       start, _Region::func_comp ) ;
        if (low == glist[chrom].end())
            return ;
        if (ir){
            for (const auto &gObj: (*low)->glist){
                if(gObj->get_start() < start  && gObj->get_end() > end){
                     oGeneList.push_back(gObj) ;
                }
            }
        } else {
            for (const auto &gObj: (*low)->glist){
                const bool stbool = (strand == '.' || strand == gObj->get_strand()) ;
                if(gObj->junc_map_.count(key) >0 && (gObj->junc_map_[key])->get_denovo_bl() && stbool){
                    oGeneList.push_back(gObj) ;
                }
            }
        }
        delete junc ;
        return ;
    }

    bool isNullJinfo(Jinfo* x){
        return (x == nullptr) ;
    }

    void free_JinfoVec(vector<Jinfo*> & jvec){
        for (const auto &jobj_ptr: jvec){
            if (jobj_ptr != nullptr) 
                delete jobj_ptr ; 
        } 
        jvec.clear() ;
    }

    void free_lsvlist(vector<LSV*> & lsvList){
        for (auto &lsv: lsvList){
            delete lsv ;
        }
    }

    string key_format(string gid, int coord1, int coord2, bool ir){
        const string g_string = ir? "IR:"+ gid : gid ;
        return(g_string + ":" + to_string(coord1) + "-" + to_string(coord2)) ;

    }
}


