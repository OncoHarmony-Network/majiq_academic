#ifndef MIXTURE_QUANTILE_H
#define MIXTURE_QUANTILE_H

#include "beta_mixture_container.h"
#include "numpy/ndarraytypes.h"

namespace BetaMixture {
    template<typename _npy_type>
    static void mixture_quantile_outer(char **args, npy_intp *dimensions, npy_intp* steps, void* data) {
        // outer loop dimensions and index
        npy_intp dim_broadcast = *dimensions++;
        npy_intp idx_broadcast;
        // strides on each variable for outer loop
        npy_intp str_q = *steps++;
        npy_intp str_a = *steps++;
        npy_intp str_b = *steps++;
        npy_intp str_out = *steps++;
        // core dimensions
        npy_intp dim_ab = dimensions[0];
        // what are the inner strides?
        npy_intp inner_str_a = steps[0];
        npy_intp inner_str_b = steps[1];

        // outer loop on broadcasted variables
        for (idx_broadcast = 0; idx_broadcast < dim_broadcast;
                ++idx_broadcast,
                args[0] += str_q, args[1] += str_a,
                args[2] += str_b, args[3] += str_out) {
            // get quantile we are evaluating
            const _npy_type q = detail::_access<_npy_type>(args[0]);
            // get iterator over beta distributions
            detail::BetaDistributionMixture<_npy_type> dist(dim_ab,
                    args[1], inner_str_a, args[2], inner_str_b);
            // compute quantile and save
            detail::_access<_npy_type>(args[3]) = dist.quantile(q);
        }
        return;
    }
}

#endif
