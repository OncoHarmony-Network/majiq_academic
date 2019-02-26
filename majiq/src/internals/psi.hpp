#ifndef PSI_H
#define PSI_H
#include <iostream>
#include <random>
#include <algorithm>
#include <string>
#include <math.h>
#include <omp.h>
#include "qLSV.hpp"
#include "stats/stats.hpp"

using namespace std ;

typedef vector<float> psi_distr_t ;
typedef pair<int, int> pair_int_t ;
#define PSEUDO 1e-20

inline float median(psi_distr_t a){
    const int n = a.size() ;
    sort(a.begin(), a.end()) ;
    if (n % 2 != 0) return (float)a[n/2] ;

//    const float m = (float)(a[(n-1)/2] + a[n/2])/2.0 ;
//cerr << "MEDIAN: " << m << " :: " ;
//    for (auto const &c: a)
//        cerr << c << ", " ;
//    cerr << "\n" ;
//    return m ;

    return (float)(a[(n-1)/2] + a[n/2])/2.0 ;
}

inline void get_prior_params( vector<psi_distr_t>& o_priors, int njunc, bool ir){

    int nways = ir ? njunc -1 : njunc ;
    float alpha = 1.0/ nways ;
    float fnjunc = (float)njunc ;

//    if (ir){
    if( 0){
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

inline float logsumexp(psi_distr_t& nums, size_t ct){
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

inline float logsumexp_2D(vector<psi_distr_t>& nums, size_t ct) {
    float max_exp = nums[0][0], sum = 0.0 ;
    size_t i, j ;

    for(i = 0; i < ct; i++){
        for(j = 1; j < ct; j++){
            max_exp = (nums[i][j] > max_exp) ? nums[i][j] : max_exp ;
        }
    }

    for(i = 0; i < ct; i++){
        for(j = 0; j < ct; j++){
            sum += exp(nums[i][j] - max_exp) ;
        }
    }
    return log(sum) + max_exp ;
}


inline float my_mean(const std::vector<float>& numbers) {
    if (numbers.empty())
        return std::numeric_limits<float>::quiet_NaN() ;

    return std::accumulate(numbers.begin(), numbers.end(), 0.0) / numbers.size() ;
}

inline float my_variance(const float mean, const std::vector<float>& numbers) {
    if (numbers.size() <= 1u)
        return std::numeric_limits<float>::quiet_NaN() ;

    auto const add_square = [mean](float sum, float i) {
        auto d = i - mean ;

//cerr << "sum: " << sum << " i: " << i << " mean: " << mean << "\n" ;
        return sum + d*d ;
    };
    double total = std::accumulate(numbers.begin(), numbers.end(), 0.0, add_square) ;
    return total / (numbers.size() - 1) ;
}



template <typename T, typename Compare>
vector<size_t> sort_permutation(const std::vector<T>& vec, const Compare& compare) {
    vector<size_t> p(vec.size()) ;
    iota(p.begin(), p.end(), 0) ;
    sort(p.begin(), p.end(), [&](size_t i, size_t j){ return compare(vec[i], vec[j]); }) ;
    return p ;
}

template <typename T>
vector<T> apply_permutation(const vector<T>& vec, const vector<size_t>& p) {
    vector<T> sorted_vec(vec.size()) ;
    transform(p.begin(), p.end(), sorted_vec.begin(), [&](size_t i){ return vec[i]; }) ;
    return sorted_vec ;
}

inline float calc_mupsi(const float sample, const float all_sample, float alpha, float beta){
//cerr << "sample: " << sample << " allsample: " << all_sample << " alpha: " << alpha << " beta: " << beta << "\n" ;
    return (sample + alpha) / (all_sample + alpha + beta) ;
}

inline void collapse_matrix(psi_distr_t& o_dpsi, vector<psi_distr_t>& matrix, int nbins){

    for (int i=0; i<nbins; i++){
        for (int j=0; j<nbins; j++){
            o_dpsi[j-i+(nbins-1)] += matrix[i][j];
        }
    }
}

void prob_data_sample_given_psi(float out_array[], float sample, float all_sample, psi_distr_t & psi_border,
                                int nbins, float alpha_prior, float beta_prior) ;

void psi_posterior(psiLSV* lsvObj, psi_distr_t& psi_border, int nbins) ;

void deltapsi_posterior(dpsiLSV* lsvObj, vector<psi_distr_t>& prior_matrix, psi_distr_t& psi_border, int nbins) ;

void get_samples_from_psi(float* osamps, hetLSV* lsvObj, int psi_samples, psi_distr_t& psi_border,
                          int nbins, int cidx, int fidx) ;
void get_samples_from_psi2(vector<psi_distr_t>& i_psi, float* osamps, float* o_mu_psi, float* o_postpsi,
                          int psi_samples, int j_offset, psi_distr_t& psi_border2, int njunc, int msamples, int nbins,
                          bool is_ir) ;

void test_calc(vector<psi_distr_t>& oPvals, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples, float quant) ;
void get_psi_border(psi_distr_t& psi_border, int nbins) ;

void adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins) ;
pair<float, float> calculate_beta_params(float mean, float vari) ;
void calc_mixture_pdf(psi_distr_t& o_mixpdf, vector<pair<float, float>>& beta_param, psi_distr_t& pmix,
                      psi_distr_t& psi_border, int nbins) ;

#endif