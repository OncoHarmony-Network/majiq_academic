#include <random>
#include <algorithm>
#include <string>
#include "scythestat/distributions.h"
#include <math.h>
#include "psi.hpp"

using namespace std ;

void prob_data_sample_given_psi(float* out_array, float sample, float all_sample, psi_distr_t & psi_border,
                                int nbins, float alpha_prior, float beta_prior){

    const float a = sample + alpha_prior ;
    const float b = (all_sample - sample) + beta_prior ;

    float prev = scythe::pbeta(psi_border[0], a, b) ;
//cout << "a=" << sample << "+"<< alpha_prior <<" b=" << b << " nbins=" << nbins << " prev=" << prev << " betap="<< beta_prior <<"\n" ;
    for (int i=0; i<nbins; i++){
        float res = scythe::pbeta(psi_border[i+1], a, b) ;
        out_array[i] = log((res - prev) + 1e-300);
        prev = res ;
    }
}

void psi_posterior(vector<psi_distr_t>& i_psi, float* o_mupsi, float* o_postpsi,
                   int msamples, int njunc, int nbins, bool is_ir){

cout << "psi_posterior 01 msamples: "<< msamples << " njunc: "<< njunc << " ir?: "<< is_ir << "\n" ;

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
cout << "IN THE LOOP\n" ;
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
        o_mupsi[j] = median(temp_mupsi) ;
        for (int i=0; i<nbins; i++){
            const int idx_2d = (j*nbins) + i ;
            o_postpsi[idx_2d] /= msamples ;
        }
    }
    free(all_m) ;
cout << "OUT LOOP\n" ;
}

//
////TODO: 1 prior_matrix should be in log space
////TODO: 2 fix mu_psi1 that should be a set of elements and median
//// All arrays are precreated
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
            temp_mupsi2[m] = calc_mupsi(jnc_val, all_val, alpha, beta) ;
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
                    A[idx_2d] = (psi_lkh1[x] * psi_lkh2[y]) + prior_matrix[idx_2d] ;
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
//
//
//
//vector<vector<vector<float>>> heterogen_posterior(float psi[], float mu_psi[]){
//
//    vector<vector<vector<float>>> samps(num_ways,  vector<vector<float>>(psi_samples, vector<float>())) ;
//   alpha_beta_prior = get_prior_params(num_ways, ir) ;
//
//    const float bsize = 1.0 / nbins ;
//    vector<float> psi_border vector(nbins) ;
//    for(int i=0; i<=nbins; i++){
//        psi_border[i] = i*bsize ;
//    }
//    vector<float> all_m vector(m_samples) ;
//    for (int j=0; j<num_ways; j++){
//        for (int i=0; i<m_samples; i++) {
//            all_m[m] += psi[j, m] ;
//        }
//    }
//    for (int j=0; j<njunc; j++){
//        const float alpha = alpha_beta_prior[j][0] ;
//        const float beta = alpha_beta_prior[j][1] ;
//        for (int m=0; m<msamples; m++){
//
//            const float jnc_val = psi[p_id, m] ;
//            const float all_val = all_m[m] ;
//            vector<float> psi_lkh vector(nbins) ;
//            mu_psi[j] += calc_mupsi(jnc_val, all_val, alpha, beta) ;
//            prob_data_sample_given_psi(psi_lkh, jnc_val, all_val, psi_border, nbins, alpha, beta) ;
//            psi_lkh = log(psi_lkh) ;
//            const float Z = logsumexp(psi_lkh, nbins) ;
//            for (i=0; i< nbins; i++){
//                psi_lkh[i] -= Z ;
//                post_psi[j][i] += exp(psi_lkh[i]) ;
//            }
//            free(psi_lkh) ;
//        }
//        for (i=0; i<nbins; i++){
//            post_psi[j][i] /= msamples ;
//        }
//
//        if (psi_samples == 1){
//
//        }else{
//            default_random_engine generator ;
//            discrete_distribution<float> psi_distribution (post_psi.begin(), post_psi.end());
//            for(i=0; i<psi_samples; i++){
//                float p = psi_distribution(generator) ;
//                samps[j][i].push_back(p) ;
//            }
//        }
////        _samples_from_psi(post_psi, mu_psi[exp, p_idx], vwindow, psi_samples, nbins)
//
//    }
//    free(all_m) ;
//    return samps ;
//}
//
//
//void test_calc(vector<vector<vector<float>>> samples1, vector<vector<vector<float>>> samples2, list<string> stats,
//               int njunc, int psamples){
//
//    const int nstats = stats.size() ;
//
//    for (int j=0; j<njunc; j++){
//        for(int s=0; s<psamples; s++){
//            vector<float> csamps ;
//            vector<int> labels ;
//
//            for( const auto& v: samples1[j][s]){
//                csamps.push_back(v) ;
//                labels.push_back(0) ;
//            }
//            for( const auto& v: samples2[j][s]){
//                csamps.push_back(v) ;
//                labels.push_back(1) ;
//            }
//
//            auto p = sort_permutation(csamps, [](T const& a, T const& b){ return(a<=b) });
//
//            vectorA = apply_permutation(vectorA, p);
//            vectorB = apply_permutation(vectorB, p);
//
//            for(int i=0; i<nstats; i++){
//                stats[i].operator(csamps, clabels) ;
//
//            }
//        }
//
//    }
//
//}
//
//
//
//
//
//
//
//
//
//    cts = [insamps[0].shape[0], insamps[1].shape[0]]
//    # lsv_id, lsv_type = lsv_attrs
//    # samps = []
//    # cts = []
//    # for grp in samples.values():
//    #     cts.append(len(grp))
//    #     samps.append(grp)
//    samps = np.concatenate(insamps, axis=0)
//    labels = np.concatenate((np.zeros(cts[0]), np.ones(cts[1])), axis=0).astype(int)
//    nfiles, njuncs, nsamples = samps.shape
//    outstats = np.ones(shape=(njuncs, len(stats)))
//    if all([nsamps >= minsamps for nsamps in cts]):
//        for jn_idx in range(njuncs):
//            cur_output = np.zeros(shape=(len(stats), nsamples))
//            for samp_idx in range(nsamples):
//                csamps = samps[:, jn_idx, samp_idx]
//                clabels = labels[csamps != -1]
//                npos = clabels.sum()
//                nneg = clabels.size - npos
//                if npos < minsamps or nneg < minsamps:
//                    continue
//                csamps = csamps[csamps != -1]
//                asort = csamps.argsort()
//                csamps = csamps[asort]
//                clabels = clabels[asort]
//                for stat_idx, stat_name in enumerate(stats):
//                    cur_output[stat_idx, samp_idx] = operator[stat_name].operator(csamps, clabels)
