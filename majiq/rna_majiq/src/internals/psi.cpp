#include <random>
#include <algorithm>
#include <string>
#include "boost/math/distributions/beta.hpp"
#include "boost/random/beta_distribution.hpp"
#include <math.h>
#include "psi.hpp"
#include "majiq_utils.hpp"
#include "qLSV.hpp"
#include "stats/stats.hpp"

#include <iostream>
#include <fstream>

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
    boost::math::beta_distribution<float> beta_dist(a, b);
//cerr << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " betap="<< beta_prior <<"\n" ;
    float prev = boost::math::cdf(beta_dist, psi_border[0]);
    for (int i=0; i<nbins; i++){
        float res = boost::math::cdf(beta_dist, psi_border[i + 1]);
        out_array[i] = log((res - prev) + PSEUDO) ;
        prev = res ;
    }
}


void psi_posterior(psiLSV*lsvObj, psi_distr_t& psi_border, int nbins){

//cout << "psi_posterior 01 msamples: "<< msamples << " njunc: "<< njunc << " ir?: "<< is_ir << "\n" ;
    const int njunc = lsvObj->get_num_ways() ;
    const int msamples = lsvObj->samps[0].size() ;

    // prior for LSV with njunc connections
    const std::pair<float, float> alpha_beta_prior = get_prior_params(njunc);
    const float alpha = alpha_beta_prior.first;
    const float beta = alpha_beta_prior.second;

    psi_distr_t all_m(msamples, 0.0) ;
    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m[m] += lsvObj->samps[j][m] ;
        }
    }

    for (int j=0; j<njunc; j++){
        psi_distr_t temp_mupsi(msamples) ;
        for (int m=0; m<msamples; m++){

            const float jnc_val = lsvObj->samps[j][m]  ;
            const float all_val = all_m[m] ;
            psi_distr_t psi_lkh (nbins, 0.0) ;
            temp_mupsi[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            prob_data_sample_given_psi(psi_lkh, jnc_val, all_val, psi_border, nbins, alpha, beta) ;
            using majiq::logsumexp;
            const float Z = logsumexp(psi_lkh.begin(), psi_lkh.end());
            for (int i=0; i< nbins; i++){
                psi_lkh[i] -= Z ;
                lsvObj->post_psi[j][i] += exp(psi_lkh[i]) ;
            }
        }
        using majiq::median;
        lsvObj->mu_psi[j] = median(temp_mupsi.begin(), temp_mupsi.end());
        for (int i=0; i<nbins; i++){
            lsvObj->post_psi[j][i] /= msamples ;
        }
    }
    all_m.clear() ;
    lsvObj->clear_samps() ;
}

void deltapsi_posterior(dpsiLSV*lsvObj, vector<psi_distr_t>& prior_matrix, psi_distr_t& psi_border, int nbins){

    const int njunc    = lsvObj->get_num_ways() ;
    const int msamples = lsvObj->cond_sample1[0].size() ;

    // prior for LSV with njunc connections
    const std::pair<float, float> alpha_beta_prior = get_prior_params(njunc);
    const float alpha = alpha_beta_prior.first;
    const float beta = alpha_beta_prior.second;

    psi_distr_t all_m1(msamples, 0.0) ;
    psi_distr_t all_m2(msamples, 0.0) ;

    for (int j=0; j<njunc; j++){
        for (int m=0; m<msamples; m++) {
            all_m1[m] += lsvObj->cond_sample1[j][m] ;
            all_m2[m] += lsvObj->cond_sample2[j][m] ;
        }
    }

    for (int j=0; j<njunc; j++){
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

            using majiq::logsumexp;
            const float Z1 = logsumexp(psi_lkh1.begin(), psi_lkh1.end());
            const float Z2 = logsumexp(psi_lkh2.begin(), psi_lkh2.end());
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

            using majiq::logsumexp_2D;
            const float Z = logsumexp_2D(A);
            for (int x=0; x<nbins; x++) {
                for (int y=0; y<nbins; y++) {
                    A[x][y] -= Z ;
                    dpsi_matrix[x][y] += exp(A[x][y]) ;
                }
            }
        }
        using majiq::median;
        lsvObj->mu_psi1[j] = median(temp_mupsi1.begin(), temp_mupsi1.end());
        lsvObj->mu_psi2[j] = median(temp_mupsi2.begin(), temp_mupsi2.end());

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


/**
 * Given input LSV reads, update statistics in lsvObj and osamps
 *
 * @param osamps array of output psi samples
 * @param lsvObj pointer to LSV with readrate information that will have
 * @param psi_samples number of samples to generate per junction for stats
 * @param visualization_samples number of samples to generate per junction for visualization
 * @param psi_border ends of bins for visualization/discretizing posterior
 * @param nbins number of bins for discretizing posterior
 * @param cidx, fidx comparison/file index
 * @param generator source of randomness for sampling distributions
 */
void get_samples_from_psi(
    float* osamps, hetLSV* lsvObj, int psi_samples, int visualization_samples,
    psi_distr_t& psi_border, int nbins, int cidx, int fidx, std::mt19937 &generator
) {
    const int osamps_shape1 = psi_samples + 1;  // shape of axis 1 of osamps
    const int njunc = lsvObj->get_num_ways() ;
    const int j_offset = lsvObj->get_junction_index() ;
    const int msamples = lsvObj->samps[0].size() ;
    const int total_samples = std::max(psi_samples, visualization_samples);
    // uniform distribution over the bootstrap replicates
    uniform_int_distribution<int> source_distribution(0, msamples - 1);

    if (!lsvObj->is_enabled()){
        // this experiment should not be quantified
        // Use -1 (generally, negative value) to indicate missing value
        for (int j=0; j<njunc; j++){
            // E(PSI) for junction is -1
            lsvObj->mu_psi[cidx][fidx][j] = -1.;
            // output psisamples should be -1
            for (int i = 0; i < osamps_shape1; ++i) {
                // j_offset is first junction in LSV, 2d index --> 1d index
                const int idx_2d = ((j + j_offset) * osamps_shape1) + i;
                // set osamps to -1
                osamps[idx_2d] = -1.;
            }
        }
        // we do not want to do any further quantification in this case
        return;
    }

    // prior for LSV with njunc connections
    const std::pair<float, float> alpha_beta_prior = get_prior_params(njunc);
    const float alpha = alpha_beta_prior.first;
    const float beta = alpha_beta_prior.second;

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
        // mean/samping distribution per bootstrap replicate
        psi_distr_t temp_mupsi(msamples) ;
        vector<boost::random::beta_distribution<float>> distributions;
        distributions.reserve(msamples);
        for (int m = 0; m < msamples; ++m) {
            // per junction vs total read counts
            const float jnc_val = lsvObj->samps[j][m] ;
            const float all_val = all_m[m] ;
            // obtain posterior mean for current replicate
            temp_mupsi[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
            // obtain parameters for bootstrap replicate distribution
            const float a = jnc_val + alpha;
            const float b = (all_val - jnc_val) + beta;
            distributions.push_back(boost::random::beta_distribution<float>(a, b));
        }
        // get median of posterior means for this junction and save it
        using majiq::median;
        const float psi_mean = median(temp_mupsi.begin(), temp_mupsi.end());
        lsvObj->mu_psi[cidx][fidx][j] = psi_mean;
        osamps[(j + j_offset) * osamps_shape1] = psi_mean;  // first value is psi mean

        // perform sampling for testing/visualization ({psi,visualization}_samples)
        psi_distr_t &visualization_distribution = lsvObj->post_psi[cidx][j];
        float visualization_impact = 1. / total_samples;  // to distribution
        for (int i = 0; i < total_samples; ++i) {
            // sample psi from uniform mixture of betas distribution
            const float sampled_psi = distributions[source_distribution(generator)](generator);
            // track sampled psi for statistical tests?
            if (i < psi_samples) {
                // we are sampling for test statistic and this is sample to use
                // get index for output sample (off by 1 since first value is psi_mean
                const int idx_2d = ((j + j_offset) * osamps_shape1) + i + 1;
                // generate output: sample from uniformly selected bootstrap distribution
                osamps[idx_2d] = sampled_psi;
            }
            // track sampled psi for visualization always
            // get position of sampled_psi among bins
            int idx_bin = static_cast<int>(
                std::upper_bound(psi_border.begin(), psi_border.end(), sampled_psi)
                - psi_border.begin()
            ) - 1;
            if (idx_bin < 0) idx_bin = 0;
            if (idx_bin >= nbins) idx_bin = nbins - 1;
            // update unnormalized visualization distribution
            visualization_distribution[idx_bin] += visualization_impact;
        }
    }
    all_m.clear() ;
    return ;
}

void test_calc(vector<psi_distr_t>& mean_pvalues, vector<psi_distr_t>& sample_pvalues, psi_distr_t& oScore, HetStats* HetStatsObj, hetLSV* lsvObj,
                                                                        int psamples, float quant){

    const int nstats = (HetStatsObj->statistics).size() ;
    const int n1 = lsvObj->cond_sample1.size() ;
    const int n2 = lsvObj->cond_sample2.size() ;
    const int njunc = lsvObj->get_num_ways() ;

    for (int j=0; j<njunc; j++){

        vector<vector<float>> pval_vect (nstats, vector<float>(psamples + 1)) ;
        psi_distr_t score_vect(psamples + 1) ;
        for(int s=0; s<psamples + 1; s++){

            vector<float> csamps ;
            vector<int> labels ;
            // reserve space in advance to avoid expanding vector multiple times
            csamps.reserve(n1 + n2);
            labels.reserve(n1 + n2);

            for (int i=0; i<n1; i++){
                if (lsvObj->cond_sample1[i][j][s] >= 0) {
                    // add only nonmissing quantifications
                    csamps.push_back(lsvObj->cond_sample1[i][j][s]);
                    labels.push_back(0);
                }
            }

            for (int i=0; i<n2; i++){
                if (lsvObj->cond_sample2[i][j][s] >= 0) {
                    // add only nonmissing quantifications
                    csamps.push_back(lsvObj->cond_sample2[i][j][s]);
                    labels.push_back(1);
                }
            }

            auto p = sort_permutation <float>(csamps, less<float>() ) ;
            csamps = apply_permutation(csamps, p);
            labels = apply_permutation(labels, p);


            for(int i=0; i<nstats; i++){
                pval_vect[i][s] = (float)(HetStatsObj->statistics)[i]->Calc_pval(csamps, labels, &score_vect[s]) ;
            }
        }
        for(int i=0; i<nstats; i++){
            mean_pvalues[j][i] = pval_vect[i][0];  // pvalue from mean
            using majiq::quantile;
            sample_pvalues[j][i] = quantile(pval_vect[i].begin() + 1, pval_vect[i].end(), quant);
            if ((HetStatsObj->names)[i] == "TNOM"){
                using majiq::median;
                float ss = median(score_vect.begin() + 1, score_vect.end());
                oScore[j] = ss ;
            }
        }
    }
}


void calc_mixture_pdf(psi_distr_t& o_mixpdf, vector<pair<float, float>>& beta_param, psi_distr_t& pmix,
                      psi_distr_t& psi_border, int nbins){

    int bpara_idx = 0 ;
    float sum = 0.0 ;
    for (const auto &bparam: beta_param){
        psi_distr_t bincdf(nbins, 0.0) ;
        prob_data_sample_given_psi(bincdf, 0.0, 0.0, psi_border, nbins, bparam.first, bparam.second) ;
        for (int i=0; i<nbins; i++){
            const float k = exp(bincdf[i]) ; //* pmix[bpara_idx] ;
            o_mixpdf[i] += k ;
            sum += k ;
        }
        bpara_idx ++ ;
    }

    for (int i=0; i<nbins; i++){
        o_mixpdf[i] /= sum ;
    }
}


pair<float, float> calculate_beta_params(float mean, float vari){
    float p, a , b ;

    p = ((mean*(1 - mean)) / vari) - 1 ;
    a = mean * p ;
    b = (1 - mean) * p ;
cerr << "BETA PARAMS: " << "a: " << a << " b: " << b << " p: " << p  << " vari: " << vari  << " mean: " << mean << "\n" ;
    return pair<float, float>(a, b) ;
}

int adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins){

    psi_distr_t psi_border(nbins+1) ;
    psi_distr_t dpsi_border(nbins+1) ;
    psi_distr_t hist(nbins, 0.0) ;
    psi_distr_t center_dst, spike_dst ;
    psi_distr_t p_mixture (3) ;

    vector<pair<float, float>> beta_params (3) ;

    const float bsize = 2.0 / nbins ;

    for(int i=0; i<=nbins; i++){
        dpsi_border[i] = i*bsize - 0.975;
        psi_border[i] = (dpsi_border[i] +1) /2 ;
    }

    psi_border[0] = 0 ;
    psi_border[nbins] = 1 ;
    dpsi_border[0] = -1 ;
    dpsi_border[nbins] = 1 ;

    sort(emp_dpsi.begin(), emp_dpsi.end()) ;
    int i = 0 ;
    int ub = 26, lb = 21 ;




    for(const auto &v: emp_dpsi){
        if (v>= dpsi_border[i] && v< dpsi_border[i+1]){
            hist[i] ++ ;
        }
        else if (v< dpsi_border[i]){
            cerr<< "THIS SHOULD NOT HAPPEN\n" ;
        }
        else{
            while(i < nbins && v>= dpsi_border[i+1]){
                i++ ;
            }
        }
        if(abs(v)<= dpsi_border[ub] && abs(v)> dpsi_border[lb]){
            center_dst.push_back((v+1)/2) ;
        }
        else if(abs(v)<= dpsi_border[lb]){
            spike_dst.push_back((v+1)/2) ;
        }
    }

    float total_cnt =  emp_dpsi.size() ;
    float spike_cnt =  spike_dst.size() ;
    float center_cnt = center_dst.size() ;


    if (spike_cnt == 0 || center_cnt == 0){
        return -1 ;
    }


    p_mixture[0] = (total_cnt - (spike_cnt + center_cnt)) / total_cnt;
    p_mixture[1] = center_cnt / total_cnt ;
    p_mixture[2] = spike_cnt / total_cnt ;

    beta_params[0] =  std::make_pair(1, 1) ;

    beta_params[1] = calculate_beta_params(0.5,
            majiq::variance<1>(center_dst.begin(), center_dst.end()));
    beta_params[2] = calculate_beta_params(0.5,
            majiq::variance<1>(spike_dst.begin(), spike_dst.end()));

    calc_mixture_pdf(o_mixtpdf, beta_params, p_mixture, psi_border, nbins) ;

    return 0 ;
}


