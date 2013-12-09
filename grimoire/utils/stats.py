import numpy as np
import globals

def junction_expression( junc_pos_array, gc, gc_factor_per_pos, k ):

    num_pos = junc_pos_array.shape[0]
    print "VAL",num_pos
    prob = np.zeros(num_pos,dtype = np.dtype('float'))
    for yy,xx in enumerate(junc_pos_array):
            frac = globals.weigh_factor[xx]
            if xx == 0:
                prob[yy] = frac
            elif xx> 0 :
                prob[yy] = xx
    sm = prob.sum()
    prob = prob/sm
    multi_res = np.random.multinomial(k*num_pos,prob)
    
    junc_cov = 0
    print "MULTI",multi_res
    print "GC", gc
    for jj,mm in enumerate(multi_res) :
        if junc_pos_array[jj] <= 0 : continue
        junc_cov += junc_pos_array[jj] * mm * gc_factor_per_pos(gc[jj])
    
    return junc_cov




