#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <algorithm>
#include <string>
#include <math.h>
#include <omp.h>
#include "psi.hpp"

void get_aggr_coverage(map<string, vector<psi_distr_t>>& output, string lsv_id, float * coverage,
                                                                        int njunc, int msamples){

//                if weight_fname != "":
//                    cov = cov * weights[lsv_id][fidx]
    if (output.count(lsv_id) > 0){
        int cidx = 0 ;
        for(int xx=0; xx< njunc; xx++){
            for(int yy=0; yy< njunc; yy++){

                output[lsv_id][xx][yy] = output[lsv_id][xx][yy] + coverage[cidx] ;
                cidx ++ ;
            }
        }
    }else{
        #pragma omp critical
        {
            output[lsv_id] = vector<psi_distr_t>(njunc, psi_distr_t(msamples)) ;
        }
        int cidx = 0 ;
        for(int xx=0; xx< njunc; xx++){
            for(int yy=0; yy< njunc; yy++){
                output[lsv_id][xx][yy] = coverage[cidx] ;
                cidx ++ ;
            }
        }
    }
}
#endif