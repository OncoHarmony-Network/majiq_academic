#include <random>
#include <algorithm>
#include <string>
#include "scythestat/distributions.h"
#include <math.h>
#include "psi.hpp"
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
//cerr << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " betap="<< beta_prior <<"\n" ;
    float prev = scythe::pbeta(psi_border[0], a, b) ;
    for (int i=0; i<nbins; i++){
        float res = scythe::pbeta(psi_border[i+1], a, b) ;
        out_array[i] = log((res - prev) + PSEUDO) ;
        prev = res ;
    }

//    cerr<< "PSI: " ;
//    for (int i=0; i<nbins; i++)
//        cerr << "[" << psi_border[i] <<"] " << out_array[i] ;
//    cerr << "[" << psi_border[nbins]  << "]\n" ;

}


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

    const int njunc    = lsvObj->get_num_ways() ;
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
        const float beta  = alpha_beta_prior[j][1] ;

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

//cerr << "GET PRIOR ";
//    for (const auto &c: alpha_beta_prior){
//cerr << "alpha: " << c[0] << " beta: " << c[1] << "\n" ;
//
//    }

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
//cerr << "BETA PARAMS: " << "a: " << a << " b: " << b << " p: " << p  << " temp: " << t <<
//        " vari: " << vari  << " mean: " << mean << "\n" ;
    return pair<float, float>(a, b) ;
}

void print_mixture( psi_distr_t& dpsi_mean, psi_distr_t& hst, psi_distr_t& mixt){

    for (long unsigned int i=0; i < hst.size(); i++){
        const float  dpsi = (dpsi_mean[i] + dpsi_mean[i+1] ) /2 ;
        cerr << dpsi << ", " << hst[i] << ", " << mixt[i]<< "\n" ;
    }
}

void print_vals( psi_distr_t& hst, string fname ){

    ofstream myfile ;
    myfile.open (fname) ;

    for (long unsigned int i=0; i < hst.size(); i++){
        myfile << hst[i]<< "\n" ;
    }
    myfile.close() ;
}


void adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins){

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

//    cerr<< "PSI\n" ;
//    for (int i=0; i<=nbins; i++)
//   cerr << i << ": " << dpsi_border[i] << "\n" ;


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
    print_vals(center_dst, "./center_dst.csv") ;
    print_vals(spike_dst, "./spike_dst.csv") ;

    float total_cnt =  emp_dpsi.size() ;
    float spike_cnt =  spike_dst.size() ;
    float center_cnt = center_dst.size() ;

    p_mixture[0] = (total_cnt - (spike_cnt + center_cnt)) / total_cnt;
    p_mixture[1] = center_cnt / total_cnt ;
    p_mixture[2] = spike_cnt / total_cnt ;

    beta_params[0] =  std::make_pair(1, 1) ;

    float cnt_mean = my_mean(center_dst) ;
    float spk_mean = my_mean(spike_dst) ;
    beta_params[1] = calculate_beta_params(0.5, my_variance(cnt_mean, center_dst)) ;
    beta_params[2] = calculate_beta_params(0.5, my_variance(spk_mean, spike_dst)) ;


//     _em_beta_mix(D, p_mixture, beta_params, num_iter, min_ratio=1e-5, logger=logger, plotpath=plotpath, nj=njunc,
//                 labels=labels)

    calc_mixture_pdf(o_mixtpdf, beta_params, p_mixture, psi_border, nbins) ;

//    print_mixture(dpsi_border,  hist, o_mixtpdf) ;


}


