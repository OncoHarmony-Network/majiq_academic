#ifndef MIXTURE_MOMENTS_H
#define MIXTURE_MOMENTS_H

#include "beta_mixture_container.h"
#include "numpy/ndarraytypes.h"

namespace BetaMixture {
    /**
     * outer loop of numpy gufunc computing moments(a, b) with signature (M)(M)->()
     */
    template<typename _npy_type>
    static void mixture_moments_outer(char **args, npy_intp *dimensions, npy_intp* steps, void* data) {
        // outer loop dimensions and index
        npy_intp dim_broadcast = *dimensions++;
        npy_intp idx_broadcast;
        // strides on each variable for outer loop
        npy_intp str_a = *steps++;
        npy_intp str_b = *steps++;
        npy_intp str_mean = *steps++;
        npy_intp str_variance = *steps++;
        // core dimensions
        npy_intp dim_ab = dimensions[0];
        // what are the inner strides?
        npy_intp inner_str_a = steps[0];
        npy_intp inner_str_b = steps[1];

        // outer loop on broadcasted variables
        for (idx_broadcast = 0; idx_broadcast < dim_broadcast;
                ++idx_broadcast,
                args[0] += str_a, args[1] += str_b,
                args[2] += str_mean, args[3] += str_variance) {
            // get iterator over beta distributions
            detail::BetaDistributionMixture<_npy_type> dist(
                    dim_ab, args[0], inner_str_a, args[1], inner_str_b);
            // compute moments
            const auto moments = dist.moments();
            // save mean/variance to appropriate locations
            detail::_access<_npy_type>(args[2]) = moments.mean;
            detail::_access<_npy_type>(args[3]) = moments.variance;
        }
        return;
    }
}

#endif
