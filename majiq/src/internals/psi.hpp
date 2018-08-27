#ifndef PSI_H
#define PSI_H
#include <iostream>
#include <random>
#include <algorithm>
#include <string>
#include <math.h>
#include <omp.h>
#include "stats/stats.hpp"

using namespace std ;

typedef vector<float> psi_distr_t ;
typedef pair<int, int> pair_int_t ;


inline float median(psi_distr_t a){
    const int n = a.size() ;
    sort(a.begin(), a.end()) ;
    if (n % 2 != 0) return (float)a[n/2] ;
    return (float)(a[(n-1)/2] + a[n/2])/2.0 ;
}

inline void get_prior_params( vector<psi_distr_t>& o_priors, int njunc, bool ir){

    int nways = ir ? njunc -1 : njunc ;
    float alpha = 1.0/ nways ;
    float fnjunc = (float)njunc ;

    if (ir){
        alpha *= (fnjunc / (fnjunc+1)) ;
        for (int i=0; i<fnjunc-1; i++){
            o_priors[i][0] = alpha ;
            o_priors[i][1] = 1 - alpha ;
        }
        o_priors[njunc-1][0] = 1 / (fnjunc + 1) ;
        o_priors[njunc-1][1] = 1 - alpha ;
    }else{
        const float beta = (fnjunc-1) / fnjunc ;
        for (int i=0; i<fnjunc; i++){
            o_priors[i][0] = alpha ;
            o_priors[i][1] = beta ;
        }
    }
    return ;
}

inline float quantile(vector<float> set, float quant){
    const int n = set.size() ;
    float idx = n * quant ;
    float c_idx = min(ceilf(idx), (float)n) ;
    if (c_idx == idx){
        int idxp = min(idx+1, (float)n) ;
        return (set[idx] + set[idxp]) / 2 ;
    }else{
        return set[c_idx] ;
    }
}


inline float logsumexp(float nums[], size_t ct){
    float max_exp = nums[0], sum = 0.0 ;
    size_t i ;

    for(i = 1; i < ct; i++){
        max_exp = (nums[i] > max_exp) ? nums[i] : max_exp ;
    }

    for(i = 0; i < ct; i++){
        sum += exp(nums[i] - max_exp) ;
    }
    return log(sum) + max_exp ;
}

template <typename T, typename Compare>
vector<size_t> sort_permutation(const std::vector<T>& vec, const Compare& compare)
{
    vector<size_t> p(vec.size()) ;
    iota(p.begin(), p.end(), 0) ;
    sort(p.begin(), p.end(), [&](size_t i, size_t j){ return compare(vec[i], vec[j]); }) ;
    return p ;
}

template <typename T>
vector<T> apply_permutation(const vector<T>& vec, const vector<size_t>& p)
{
    vector<T> sorted_vec(vec.size()) ;
    transform(p.begin(), p.end(), sorted_vec.begin(), [&](size_t i){ return vec[i]; }) ;
    return sorted_vec ;
}

inline float calc_mupsi(const float sample, const float all_sample, float alpha, float beta){
    return (sample + alpha) / (all_sample + alpha + beta) ;
}

inline void collapse_matrix(float* o_dpsi, float* matrix, int nbins){
    for (int i=0; i<nbins; i++){
       for (int j=0; j<nbins; j++){
            o_dpsi[j-i+(nbins-1)] += matrix[i*nbins + j] ;
       }
    }
}


void prob_data_sample_given_psi(float out_array[], float sample, float all_sample, psi_distr_t & psi_border,
                                int nbins, float alpha_prior, float beta_prior) ;
void psi_posterior(vector<psi_distr_t> & i_psi, float* o_mupsi, float* o_postpsi,
                   int msamples, int njunc, int nbins, bool is_ir) ;

void deltapsi_posterior(vector<psi_distr_t>& i_psi1, vector<psi_distr_t>& i_psi2, float* prior_matrix,
                        float* o_mu_psi1, float* o_mu_psi2, float* o_post_psi1, float* o_post_psi2,
                        float* o_posterior_dpsi, int msamples, int njunc, int nbins, bool is_ir) ;

void get_samples_from_psi(vector<psi_distr_t>& i_psi, float* osamps, float* o_mu_psi, float* o_postpsi1,
                          int psi_samples, int j_offset, psi_distr_t psi_border, int njunc, int msamples, int nbins,
                          bool is_ir) ;

void test_calc(float* oPvals, vector<float*> samples1, vector<float*> samples2, HetStats* HetStatsObj,
               int njunc, int psamples, float quant) ;

psi_distr_t get_psi_border(int nbins) ;

#endif