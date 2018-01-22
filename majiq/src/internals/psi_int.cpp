#include <boost/math/special_functions/beta.hpp>

float * _prob_data_sample_given_psi(float sample, float all_sample, float alpha_0, float beta_0){
    float a = sample + alpha_0;
    float b = (all_sample - sample) + beta_0;
    float* bincdf

    iterator_range< range_detail::integer_iterator_with_step<Integer, StepSize> > x = irange(0, 1.01, bsize);
    bincdf = ibetac(a, b, x)
    float bin_test = project(bincdf, 0, bsize-1) - project(bincdf, 1, bsize);
    return bin_test
}





void delta_posterior(float ** data1, float ** data2, float ** prior, int p_idx, int msamples, int nexp,
                     int nbins, float alpha_0, float beta_0) {
//    deltapsi implementation

    for (int i; i<m_samples; i ++)
    {

    }




}



void main(){

    float * kk = _prob_data_sample_given_psi(10, 15, 0.5, 0.5);
    printf("RESULTS");
    printf(kk);

}