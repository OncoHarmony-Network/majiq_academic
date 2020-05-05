#ifndef TTEST_H
#define TTEST_H
#include "MathFunctions.h"
#include "testStats.hpp"
#include "scythestat/distributions.h"
#include <random>
#include <algorithm>
#include <string>
#include <assert.h>
#include <unordered_map>
#include <iostream>

#define MAXCLASS 2

using namespace std ;
namespace MajiqStats{

    class tTest: public MajiqStats::TestStat{

        private:
            unordered_map<tTRecord, double, hash_fn> _pval_cache ;

        public:
            double Calc_pval(vector<float>& data, vector<int>& labels, float* score){

                // calculate mean and var for the expression values and num of samples in each class.
                double Count[2] ;
                double Mean[2] ;
                double Var[2] ;
                for(int i = 0; i < MAXCLASS; i++ ) {
                    Count[i] = 0 ;
                    Mean[i] = 0 ;
                    Var[i] = 0 ;
                }
                for(int i = 0; i < (int)data.size(); i++ ){
                    if( labels[i] >= 0 ){
                        double x = data[i] ;
                        Count[labels[i]]++ ;
                        Mean[labels[i]] += x ;
                        Var[labels[i]] += x*x ;
                    }
                }
                for(int i = 0; i < 2; i++ ) {
                    if( Count[i] > 0 ) {
                        Mean[i] /= Count[i] ;
                        Var[i] /= Count[i] ;
                        Var[i] -= Mean[i]*Mean[i] ;
                        if( Var[i] < 0.0000000001 )
                            Var[i] = 0.0000000001 ;
                    }
                }
                if (Count[0] < 2 || Count[1] < 2) {
                    // unable to calculate sample variance, set p value to 1
                    return 1.;
                }
                double _Stat = sqrt(Count[0]+Count[1] - 2) * ( Mean[0] - Mean[1] ) ;
                _Stat /= sqrt((1.0/Count[0] + 1.0/Count[1])*(Count[0]*Var[0] + Count[1]*Var[1])) ;

                double _Deg = Count[0]+Count[1] - 2 ;

                // if the two distributions are equal, then U ~ T(Mean[1]+Mean[2] -2)

                double U = fabs(_Stat) ;
//                double P, Q ;
//                cumt( &U, &_Deg, &P, &Q ) ;
                double Q = scythe::pt(U, _Deg) ;
//                double Q = 0 ;
                return (1-Q) ;
            }
    };
}
#endif
