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
 * Obtains median of values between first and last
 *
 * @param first, last random-access iterators. Note: values will be reordered
 */
template<class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type median(It first, It last) {
    // should be random iterator
    static_assert(std::is_same<
            std::random_access_iterator_tag, typename std::iterator_traits<It>::iterator_category
            >::value,
            "median() only accepts random-access iterators as input\n");
    // how many elements?
    const auto n = std::distance(first, last);
    // 3 cases
    if (n < 1) {
        return std::numeric_limits<typename std::iterator_traits<It>::value_type>::quiet_NaN();
    } else if (n == 1) {
        return *first;
    } else if (n == 2) {
        // no need to sort anything
        return (first[0] + first[1]) / 2;
    } else {
        // sort to get the value that would be in the middle if sorted
        It it_below = first + ((n - 1) / 2);
        std::nth_element(first, it_below, last);
        if (n % 2 == 1) {
            // odd means that this is the median
            return *it_below;
        } else {
            // we need the value that follows and take the average
            return (*it_below + *std::min_element(it_below + 1, last)) / 2;
        }
    }
}

/**
 * Beta-distribution prior (1 / njunc, (njunc - 1) / njunc)
 */
inline std::pair<float, float> get_prior_params(int njunc) {
    float alpha = 1.0 / njunc;
    float beta = 1.0 - alpha;
    return std::make_pair(alpha, beta);
}

/**
 * Obtains the specified quantile of values between first and last
 *
 * @param first, last random-access iterators. Note: values will be reordered
 * @param quantile the quantile to evaluate
 *
 * Obtains quantile between observations using linear interpolation (i.e. NumPy
 * default)
 */
template<class It>
inline
#if __cplusplus >= 201703L
constexpr  // required for c++17 and later
#endif
typename std::iterator_traits<It>::value_type quantile(It first, It last, float quantile) {
    // should be random iterator
    static_assert(std::is_same<
            std::random_access_iterator_tag, typename std::iterator_traits<It>::iterator_category
            >::value,
            "quantile() only accepts random-access iterators as input\n");
    // how many elements?
    const auto n = std::distance(first, last);
    if (n < 1) {
        return std::numeric_limits<typename std::iterator_traits<It>::value_type>::quiet_NaN();
    } else if (n == 1) {
        return *first;
    } else {
        // otherwise, we need to determine where we belong between 0 and n-1
        if (quantile <= 0) {
            return *std::min_element(first, last);
        } else if (quantile >= 1) {
            return *std::max_element(first, last);
        } else {
            // between which indexes in sorted array is our quantile?
            const float idx_float = (n - 1) * quantile;
            const float idx_floor_float = std::floor(idx_float);
            It it_floor = first + static_cast<int>(idx_floor_float);
            // partial sort of values to put correct value at it_floor
            std::nth_element(first, it_floor, last);
            // get values around idx_float
            const auto value_below = *it_floor;
            const auto value_above = *std::min_element(it_floor + 1, last);
            // linear interpolation between the two:
            return value_below + (idx_float - idx_floor_float) * (value_above - value_below);
        }
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
