#ifndef PSI_H
#define PSI_H
#include <iostream>
#include <random>
#include <algorithm>
#include <string>
#include <math.h>
#include <cmath>
#include <omp.h>
#include "qLSV.hpp"
#include "stats/stats.hpp"

using namespace std ;

typedef vector<float> psi_distr_t ;
typedef pair<int, int> pair_int_t ;
#define PSEUDO 1e-20


/**
 * Beta-distribution prior (1 / njunc, (njunc - 1) / njunc)
 */
inline std::pair<float, float> get_prior_params(int njunc) {
    float alpha = 1.0 / njunc;
    float beta = 1.0 - alpha;
    return std::make_pair(alpha, beta);
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

void get_samples_from_psi(
    float* osamps, hetLSV* lsvObj, int psi_samples, int visualization_samples,
    psi_distr_t& psi_border, int nbins, int cidx, int fidx, std::mt19937 &generator
);

void test_calc(vector<psi_distr_t>& mean_pvalues, vector<psi_distr_t>& sample_pvalues, psi_distr_t& oScore, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples, float quant) ;
void get_psi_border(psi_distr_t& psi_border, int nbins) ;

int adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins) ;
pair<float, float> calculate_beta_params(float mean, float vari) ;
void calc_mixture_pdf(psi_distr_t& o_mixpdf, vector<pair<float, float>>& beta_param, psi_distr_t& pmix,
                      psi_distr_t& psi_border, int nbins) ;

#endif
