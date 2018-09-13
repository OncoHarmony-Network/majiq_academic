#ifndef qLSV_H
#define qLSV_H

#include <iostream>
#include <random>
#include <algorithm>
#include <string>
#include <math.h>
#include <omp.h>
#include <vector>

using namespace std ;
typedef vector<float> psi_distr_t ;

class qLSV{

    protected:
        int nways_ ;
        bool is_ir_ ;

    public:
        vector<psi_distr_t> samps ;
        qLSV(int nways1, bool is_ir1): nways_(nways1), is_ir_(is_ir1){

        }

        ~qLSV(){
            clear_samps() ;
        }

        void clear_samps(){
            samps.clear() ;
        }

        void reset_samps(){

            for(int xx=0; xx< nways_; xx++){
                fill(samps[xx].begin(), samps[xx].end(), (float)0.0) ;
            }
        }

        void add(float* coverage, int msamples){

            if (samps.size() ==0)
                samps = vector<psi_distr_t>(nways_, psi_distr_t(msamples)) ;

//cout << "samps[" << samps.size() << ", " << samps[0].size() << "]\n" ;
//cout << "loops[" << nways_ << ", " << msamples << "]\n" ;
            int cidx = 0 ;
            for(int xx=0; xx< nways_; xx++){
                for(int yy=0; yy< msamples; yy++){
                    samps[xx][yy] = samps[xx][yy] + coverage[cidx] ;
                    cidx ++ ;
                }
            }
        }

        int get_num_ways() {
            return nways_ ;
        }

        bool is_ir(){
            return is_ir_ ;
        }

};


class hetLSV: public qLSV{
    private:
        int joffset_ ;

    public:

//TODO: move vectors to private and create setter and adders functions

        vector<vector<psi_distr_t>> mu_psi ;
        vector<vector<psi_distr_t>> post_psi ;

        vector<vector<psi_distr_t>> cond_sample1 ;
        vector<vector<psi_distr_t>> cond_sample2 ;

        hetLSV(int nways1, int j_offset1, int max_nfiles1, int nbins1, bool is_ir, int ncond1):
                                                                    qLSV(nways1, is_ir), joffset_(j_offset1){
            mu_psi = vector<vector<psi_distr_t>>(ncond1, vector<psi_distr_t>(max_nfiles1, psi_distr_t(nways1))) ;
            post_psi = vector<vector<psi_distr_t>>(ncond1, vector<psi_distr_t>(nways1, psi_distr_t(nbins1))) ;

        }

        ~hetLSV() {}

        int get_junction_index(){
            return joffset_ ;
        }

        void create_condition_samples(int size1, int size2, int psi_samples){
            cond_sample1 = vector<vector<psi_distr_t>>(size1, vector<psi_distr_t>(nways_, psi_distr_t(psi_samples))) ;
            cond_sample2 = vector<vector<psi_distr_t>>(size2, vector<psi_distr_t>(nways_, psi_distr_t(psi_samples))) ;
        }

        void add_condition1(float * k, int x, int nways, int psi_samples){
            for (int i=0; i< nways; i++){
                for (int j=0; j< psi_samples; j++){
                    int  index = (i*psi_samples) + j ;
                    cond_sample1[x][i][j] = k[index] ;
                }
            }
        }

        void add_condition2(float * k, int x, int nways, int psi_samples){
            for (int i=0; i< nways; i++){
                for (int j=0; j< psi_samples; j++){
                    int  index = (i*psi_samples) + j ;
                    cond_sample2[x][i][j] = k[index] ;
                }
            }
        }


        void clear(){
            clear_samps() ;
            mu_psi.clear() ;
            mu_psi.shrink_to_fit() ;
            post_psi.clear() ;
            post_psi.shrink_to_fit() ;
            cond_sample1.clear() ;
            cond_sample1.shrink_to_fit() ;
            cond_sample2.clear() ;
            cond_sample2.shrink_to_fit() ;
        }

};

#endif

