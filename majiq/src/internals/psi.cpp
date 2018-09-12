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

void prob_data_sample_given_psi2(psi_distr_t& out_array, float sample, float all_sample, psi_distr_t & psi_border,
                                int nbins, float alpha_prior, float beta_prior){

    const float a = sample + alpha_prior ;
    const float b = (all_sample - sample) + beta_prior ;

//cout << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " betap="<< beta_prior <<"\n" ;
    float prev = scythe::pbeta(psi_border[0], a, b) ;
    for (int i=0; i<nbins; i++){
        float res = scythe::pbeta(psi_border[i+1], a, b) ;
        out_array[i] = log((res - prev) + 1e-300) ;
        prev = res ;
    }
}


void prob_data_sample_given_psi(float* out_array, float sample, float all_sample, psi_distr_t & psi_border,
                                int nbins, float alpha_prior, float beta_prior){

    const float a = sample + alpha_prior ;
    const float b = (all_sample - sample) + beta_prior ;

//cout << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " betap="<< beta_prior <<"\n" ;
    float prev = scythe::pbeta(psi_border[0], a, b) ;
    for (int i=0; i<nbins; i++){
        float res = scythe::pbeta(psi_border[i+1], a, b) ;
        out_array[i] = log((res - prev) + 1e-300);
        prev = res ;
    }
}

void psi_posterior(vector<psi_distr_t>& i_psi, float* o_mupsi, float* o_postpsi,
                   int msamples, int njunc, int nbins, bool is_ir){

//cout << "psi_posterior 01 msamples: "<< msamples << " njunc: "<< njunc << " ir?: "<< is_ir << "\n" ;

    vector<psi_distr_t> alpha_beta_prior(njunc, psi_distr_t(2)) ;
    get_prior_params(alpha_beta_prior, njunc, is_ir) ;

//    TODO: we can move this up to avoid recalculation
    const float bsize = 1.0 / nbins ;
    psi_distr_t psi_border (nbins+1) ;
    for(int i=0; i<=nbins; i++){
        psi_border[i] = i*bsize ;
    }

    float * all_m = (float*) calloc(msamples, sizeof(float)) ;
    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m[m] += i_psi[j][m] ;
        }
    }

    for (int j=0; j<njunc; j++){
        const float alpha = alpha_beta_prior[j][0] ;
        const float beta = alpha_beta_prior[j][1] ;
        psi_distr_t temp_mupsi(msamples) ;
        for (int m=0; m<msamples; m++){

            const float jnc_val = i_psi[j][m] ;
            const float all_val = all_m[m] ;
            float * psi_lkh =  (float*) calloc(nbins, sizeof(float)) ;
            temp_mupsi[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh, jnc_val, all_val, psi_border, nbins, alpha, beta) ;
            const float Z = logsumexp(psi_lkh, nbins) ;
            for (int i=0; i< nbins; i++){
                psi_lkh[i] -= Z ;
                const int idx_2d = (j*nbins) + i ;
                o_postpsi[idx_2d] += exp(psi_lkh[i]) ;
            }
            free(psi_lkh) ;
        }
        sort (temp_mupsi.begin(), temp_mupsi.end()) ;
        o_mupsi[j] = median(temp_mupsi) ;
        for (int i=0; i<nbins; i++){
            const int idx_2d = (j*nbins) + i ;
            o_postpsi[idx_2d] /= msamples ;
        }
    }
    free(all_m) ;
//cout << "OUT LOOP\n" ;
}


void deltapsi_posterior(vector<psi_distr_t>& i_psi1, vector<psi_distr_t>& i_psi2, float* prior_matrix,
                        float* o_mupsi1, float* o_mupsi2, float* o_postpsi1, float* o_postpsi2,
                        float* o_posterior_dpsi, int msamples, int njunc, int nbins, bool is_ir){

    vector<psi_distr_t> alpha_beta_prior(njunc, psi_distr_t(2)) ;
    get_prior_params(alpha_beta_prior, njunc, is_ir) ;

    const float bsize = 1.0 / nbins ;
    psi_distr_t psi_border (nbins+1) ;
    for(int i=0; i<=nbins; i++){
        psi_border[i] = i*bsize ;
    }

    float * all_m1 = (float*)calloc(msamples, sizeof(float)) ;
    float * all_m2 = (float*)calloc(msamples, sizeof(float)) ;
    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m1[m] += i_psi1[j][m] ;
            all_m2[m] += i_psi2[j][m] ;
        }
    }

    const int nbins_dpsi = (nbins*2) - 1 ;
    for (int j=0; j<njunc; j++){
        const float alpha = alpha_beta_prior[j][0] ;
        const float beta = alpha_beta_prior[j][1] ;
        float * dpsi_matrix = (float*)calloc(nbins*nbins, sizeof(float)) ;
        psi_distr_t temp_mupsi1(msamples) ;
        psi_distr_t temp_mupsi2(msamples) ;
        for (int m=0; m<msamples; m++){

            float jnc_val = i_psi1[j][m] ;
            float all_val = all_m1[m] ;
            float * psi_lkh1 =  (float*) calloc(nbins, sizeof(float)) ;
            temp_mupsi1[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh1, jnc_val, all_val, psi_border, nbins, alpha, beta) ;

            jnc_val = i_psi2[j][m] ;
            all_val = all_m2[m] ;
            float * psi_lkh2 =  (float*) calloc(nbins, sizeof(float)) ;
            temp_mupsi2[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh2, jnc_val, all_val, psi_border, nbins, alpha, beta) ;

            const float Z1 = logsumexp(psi_lkh1, nbins) ;
            const float Z2 = logsumexp(psi_lkh2, nbins) ;
            for (int i=0; i< nbins; i++){
                const int idx_2d = (j*nbins) + i ;
                psi_lkh1[i] -= Z1 ;
                o_postpsi1[idx_2d] += exp(psi_lkh1[i]) ;
                psi_lkh2[i] -= Z2 ;
                o_postpsi2[idx_2d] += exp(psi_lkh2[i]) ;
            }

            float * A = (float*)calloc(nbins*nbins, sizeof(float)) ;
            for (int x=0; x< nbins; x++){
                for(int y=0; y<nbins; y++){
                    const int idx_2d = (x*nbins) + y ;
                    A[idx_2d] = (psi_lkh1[x] + psi_lkh2[y]) + prior_matrix[idx_2d] ;
                }
            }

            const float Z = logsumexp(A, nbins*nbins) ;
            for (int i=0; i<(nbins*nbins); i++) {
                A[i] -= Z ;
                dpsi_matrix[i] += exp(A[i]) ;
            }
            free(psi_lkh1) ;
            free(psi_lkh2) ;
            free(A) ;

        }
        sort (temp_mupsi1.begin(), temp_mupsi1.end()) ;
        sort (temp_mupsi2.begin(), temp_mupsi2.end()) ;
        o_mupsi1[j] = median(temp_mupsi1) ;
        o_mupsi2[j] = median(temp_mupsi2) ;
        for (int i=0; i<nbins; i++){
            int idx_2d = (j*nbins) + i ;
            o_postpsi1[idx_2d] /= msamples ;
            o_postpsi2[idx_2d] /= msamples ;
            for (int i2=0; i2<nbins; i2++){
                idx_2d = (i*nbins) + i2 ;
                dpsi_matrix[idx_2d] /= msamples ;
            }
        }
        collapse_matrix(&o_posterior_dpsi[j*nbins_dpsi], dpsi_matrix, nbins) ;
        free(dpsi_matrix) ;
    }
    free(all_m1) ;
    free(all_m2) ;
}


void get_samples_from_psi(float* osamps, hetLSV* lsvObj, int psi_samples, psi_distr_t& psi_border,
                          int nbins, int cidx, int fidx){

//    cout<< "MM1\n" ;

    const int njunc = lsvObj->get_num_ways() ;
    const int j_offset = lsvObj->get_junction_index() ;
    const int msamples = lsvObj->samps[0].size() ;

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

        for (int m=0; m<msamples; m++){
            const float jnc_val = lsvObj->samps[j][m] ;
            const float all_val = all_m[m] ;
            psi_distr_t psi_lkh (nbins, 0.0) ;

            temp_mupsi[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi2(psi_lkh, jnc_val, all_val, psi_border, nbins, alpha, beta) ;
            const float Z = logsumexp(psi_lkh, nbins) ;
            for (int i=0; i< nbins; i++){
                psi_lkh[i] -= Z ;
                lsvObj->post_psi[cidx][j][i] += exp(psi_lkh[i]) ;
            }
            psi_lkh.clear() ;
        }
        sort (temp_mupsi.begin(), temp_mupsi.end()) ;
        lsvObj->mu_psi[cidx][fidx][j] = median(temp_mupsi) ;
        for (int i=0; i<nbins; i++){
            lsvObj->post_psi[cidx][j][i] /= msamples ;
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
