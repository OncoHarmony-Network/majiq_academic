#include <random>
#include <algorithm>
#include <string>
#include "scythestat/distributions.h"
#include <math.h>
#include "psi.hpp"
#include "qLSV.hpp"
#include "stats/stats.hpp"

using namespace std ;

void get_psi_border(psi_distr_t& psi_border, int nbins){
    const float bsize = 1.0 / nbins ;

    for(int i=0; i<=nbins; i++){
        psi_border[i] = i*bsize ;
    }
    return ;
}

void prob_data_sample_given_psi(psi_distr_t& out_array, float sample, float all_sample, psi_distr_t & psi_border,
                                int nbins, float alpha_prior, float beta_prior){

    const float a = sample + alpha_prior ;
    const float b = (all_sample - sample) + beta_prior ;

//cout << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " betap="<< beta_prior <<"\n" ;
    float prev = scythe::pbeta(psi_border[0], a, b) ;
    for (int i=0; i<nbins; i++){
        float res = scythe::pbeta(psi_border[i+1], a, b) ;
        out_array[i] = log((res - prev) + PSEUDO) ;
        prev = res ;
    }
}


//void prob_data_sample_given_psi(float* out_array, float sample, float all_sample, psi_distr_t & psi_border,
//                                int nbins, float alpha_prior, float beta_prior){
//
//    const float a = sample + alpha_prior ;
//    const float b = (all_sample - sample) + beta_prior ;
//
////cout << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " betap="<< beta_prior <<"\n" ;
//    float prev = scythe::pbeta(psi_border[0], a, b) ;
//    for (int i=0; i<nbins; i++){
//        float res = scythe::pbeta(psi_border[i+1], a, b) ;
//        out_array[i] = log((res - prev) + PSEUDO);
//        prev = res ;
//    }
//}

void psi_posterior(psiLSV*lsvObj, psi_distr_t& psi_border, int nbins){

//cout << "psi_posterior 01 msamples: "<< msamples << " njunc: "<< njunc << " ir?: "<< is_ir << "\n" ;

    const int njunc = lsvObj->get_num_ways() ;
    const int msamples = lsvObj->samps[0].size() ;

    vector<psi_distr_t> alpha_beta_prior(njunc, psi_distr_t(2)) ;
    get_prior_params(alpha_beta_prior, njunc, lsvObj->is_ir()) ;

    psi_distr_t all_m(msamples, 0.0) ;
    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m[m] += lsvObj->samps[j][m] ;
        }
    }

    for (int j=0; j<njunc; j++){
        const float alpha = alpha_beta_prior[j][0] ;
        const float beta = alpha_beta_prior[j][1] ;
        psi_distr_t temp_mupsi(msamples) ;
        for (int m=0; m<msamples; m++){

            const float jnc_val = lsvObj->samps[j][m]  ;
            const float all_val = all_m[m] ;
            psi_distr_t psi_lkh (nbins, 0.0) ;
            temp_mupsi[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh, jnc_val, all_val, psi_border, nbins, alpha, beta) ;
            const float Z = logsumexp(psi_lkh, nbins) ;
            for (int i=0; i< nbins; i++){
                psi_lkh[i] -= Z ;
                lsvObj->post_psi[j][i] += exp(psi_lkh[i]) ;
            }
        }
        sort (temp_mupsi.begin(), temp_mupsi.end()) ;
        lsvObj->mu_psi[j] = median(temp_mupsi) ;
        for (int i=0; i<nbins; i++){
            lsvObj->post_psi[j][i] /= msamples ;
        }
    }
    all_m.clear() ;
    lsvObj->clear_samps() ;
//cout << "OUT LOOP\n" ;
}

void deltapsi_posterior(dpsiLSV*lsvObj, vector<psi_distr_t>& prior_matrix, psi_distr_t& psi_border, int nbins){

    const int njunc = lsvObj->get_num_ways() ;
    const int msamples = lsvObj->cond_sample1[0].size() ;

    vector<psi_distr_t> alpha_beta_prior(njunc, psi_distr_t(2)) ;
    get_prior_params(alpha_beta_prior, njunc, lsvObj->is_ir()) ;

    psi_distr_t all_m1(msamples, 0.0) ;
    psi_distr_t all_m2(msamples, 0.0) ;

    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m1[m] += lsvObj->cond_sample1[j][m] ;
            all_m2[m] += lsvObj->cond_sample2[j][m] ;
        }
    }

    for (int j=0; j<njunc; j++){
        const float alpha = alpha_beta_prior[j][0] ;
        const float beta = alpha_beta_prior[j][1] ;

        vector<psi_distr_t> dpsi_matrix(nbins, psi_distr_t(nbins, 0)) ;
        psi_distr_t temp_mupsi1(msamples) ;
        psi_distr_t temp_mupsi2(msamples) ;
        for (int m=0; m<msamples; m++){
            float jnc_val = lsvObj->cond_sample1[j][m] ;
            float all_val = all_m1[m] ;
            psi_distr_t psi_lkh1 (nbins, 0.0) ;
            temp_mupsi1[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh1, jnc_val, all_val, psi_border, nbins, alpha, beta) ;

            jnc_val = lsvObj->cond_sample2[j][m] ;
            all_val = all_m2[m] ;
            psi_distr_t psi_lkh2 (nbins, 0.0) ;
            temp_mupsi2[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh2, jnc_val, all_val, psi_border, nbins, alpha, beta) ;

            const float Z1 = logsumexp(psi_lkh1, nbins) ;
            const float Z2 = logsumexp(psi_lkh2, nbins) ;
            for (int i=0; i< nbins; i++){
                psi_lkh1[i] -= Z1 ;
                lsvObj->post_psi1[j][i] += exp(psi_lkh1[i]) ;
                psi_lkh2[i] -= Z2 ;
                lsvObj->post_psi2[j][i] += exp(psi_lkh2[i]) ;
            }

            vector<psi_distr_t> A(nbins, psi_distr_t(nbins, 0)) ;
            for (int x=0; x< nbins; x++){
                for(int y=0; y<nbins; y++){
                    A[x][y] = (psi_lkh1[x] + psi_lkh2[y]) + prior_matrix[x][y] ;
                }
            }

            const float Z = logsumexp_2D(A, nbins) ;
            for (int x=0; x<nbins; x++) {
                for (int y=0; y<nbins; y++) {
                    A[x][y] -= Z ;
                    dpsi_matrix[x][y] += exp(A[x][y]) ;
                }
            }
        }
        sort (temp_mupsi1.begin(), temp_mupsi1.end()) ;
        sort (temp_mupsi2.begin(), temp_mupsi2.end()) ;
        lsvObj->mu_psi1[j] = median(temp_mupsi1) ;
        lsvObj->mu_psi2[j] = median(temp_mupsi2) ;

        for (int i=0; i<nbins; i++){
            lsvObj->post_psi1[j][i] /= msamples ;
            lsvObj->post_psi2[j][i] /= msamples ;
            for (int i2=0; i2<nbins; i2++){
                dpsi_matrix[i][i2] /= msamples ;
            }
        }

        collapse_matrix(lsvObj->post_dpsi[j], dpsi_matrix, nbins) ;
        dpsi_matrix.clear() ;
    }
    all_m1.clear() ;
    all_m2.clear() ;
    lsvObj->clear() ;
}


void get_samples_from_psi(float* osamps, hetLSV* lsvObj, int psi_samples, psi_distr_t& psi_border,
                          int nbins, int cidx, int fidx){

//    cout<< "MM1\n" ;
    const int njunc = lsvObj->get_num_ways() ;
    const int j_offset = lsvObj->get_junction_index() ;
    const int msamples = lsvObj->samps[0].size() ;

    if (!lsvObj->is_enabled()){
        for (int j=0; j<njunc; j++){
            lsvObj->mu_psi[cidx][fidx][j] = -1 ;
        }
    }

//cout<< "MM2\n" ;
//
    vector<psi_distr_t> alpha_beta_prior(njunc, psi_distr_t(2)) ;
    get_prior_params(alpha_beta_prior, njunc, lsvObj->is_ir()) ;

    psi_distr_t psi_space(nbins) ;
    for(int i=0; i<nbins; i++){
        psi_space[i] = (psi_border[i] + psi_border[i+1]) / 2 ;
    }

    psi_distr_t all_m(msamples, 0.0) ;
    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m[m] += lsvObj->samps[j][m] ;
        }
    }

   for (int j=0; j<njunc; j++){
        const float alpha = alpha_beta_prior[j][0] ;
        const float beta = alpha_beta_prior[j][1] ;
        psi_distr_t temp_mupsi(msamples) ;
        psi_distr_t temp_postpsi (nbins, 0.0) ;
        for (int m=0; m<msamples; m++){
            const float jnc_val = lsvObj->samps[j][m] ;
            const float all_val = all_m[m] ;
            psi_distr_t psi_lkh (nbins, 0.0) ;

            temp_mupsi[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh, jnc_val, all_val, psi_border, nbins, alpha, beta) ;
            const float Z = logsumexp(psi_lkh, nbins) ;
            for (int i=0; i< nbins; i++){
                psi_lkh[i] -= Z ;
                temp_postpsi[i] += exp(psi_lkh[i]) ;
//                lsvObj->post_psi[cidx][j][i] += exp(psi_lkh[i]) ;
            }
            psi_lkh.clear() ;
        }
        sort (temp_mupsi.begin(), temp_mupsi.end()) ;
        lsvObj->mu_psi[cidx][fidx][j] = median(temp_mupsi) ;
        for (int i=0; i<nbins; i++){
            temp_postpsi[i] /= msamples ;
            lsvObj->post_psi[cidx][j][i] += temp_postpsi[i] ;
//            lsvObj->post_psi[cidx][j][i] /= msamples ;
        }
        if (psi_samples == 1){
            osamps[j+j_offset] = lsvObj->mu_psi[cidx][fidx][j] ;
        }else{
            default_random_engine generator ;
            discrete_distribution<int> psi_distribution (lsvObj->post_psi[cidx][j].begin(),
                                                         lsvObj->post_psi[cidx][j].end());
            for(int i=0; i<psi_samples; i++){
                float p = psi_distribution(generator) ;
                const int idx_2d = ((j+j_offset)*psi_samples) + i ;
                osamps[idx_2d] = psi_space[p] ;
            }
        }
    }
    all_m.clear() ;
    return ;

}

void test_calc(vector<psi_distr_t>& oPvals, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples, float quant){
//void test_calc(float* oPvals, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples, float quant){

    const int nstats = (HetStatsObj->statistics).size() ;
    const int n1 = lsvObj->cond_sample1.size() ;
    const int n2 = lsvObj->cond_sample2.size() ;
    const int njunc = lsvObj->get_num_ways() ;

    for (int j=0; j<njunc; j++){

        vector<vector<float>> pval_vect (nstats, vector<float>(psamples)) ;
        for(int s=0; s<psamples; s++){

            default_random_engine generator;
            normal_distribution<double> dist(0.002, 0.001);


            vector<float> csamps ;
            vector<int> labels ;

            for (int i=0; i<n1; i++){
                csamps.push_back(lsvObj->cond_sample1[i][j][s] + dist(generator)) ;
                labels.push_back(0) ;
            }

            for (int i=0; i<n2; i++){

                csamps.push_back(lsvObj->cond_sample2[i][j][s] + dist(generator)) ;
                labels.push_back(1) ;
            }

            auto p = sort_permutation <float>(csamps, less<float>() ) ;
            csamps = apply_permutation(csamps, p);
            labels = apply_permutation(labels, p);




            for(int i=0; i<nstats; i++){
                pval_vect[i][s] = (float)(HetStatsObj->statistics)[i]->Calc_pval(csamps, labels) ;
            }
        }
        for(int i=0; i<nstats; i++){

            float mmm = quantile(pval_vect[i], quant) ;
            oPvals[j][i] = mmm ;

        }
    }
}
