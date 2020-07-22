#ifndef TNOM_H
#define TNOM_H
#include "MathFunctions.h"
#include <random>
#include <algorithm>
#include <string>
#include <unordered_map>
#include "testStats.hpp"
#include <iostream>
#define MAXCLASS 2

using namespace std ;
namespace MajiqStats{

    class TNOM: public MajiqStats::TestStat{
        private:

            unordered_map<tTRecord, double, hash_fn> _pval_cache ;

            /**
            * compute for the classes distribution: number of "mistakes" in this side of the vector:
            * number of samples - number of samples in biggest class in this side.
            * for example: the vector ++---|-+++ will get 2 for left side, 1 for right side.
            * @param N - int vector hold num of samples for each class
            * @return 'Sum - Max'
            **/
            double Loss(int *N) const {
                double Sum = 0 ;
                double Max = 0.0 ;
                for( int i = 0; i < MAXCLASS ; i++ ){
                    Sum += abs(N[i]) ;
                    if( N[i] > Max )
                        Max = N[i] ;
                }
                return Sum - Max;
            }

            // Count number of paths from (0,0) to (Len,Offset)
            double LogPathNum(int Len, int Offset) {
                if( Offset < 0 )
                    Offset = -Offset;
                if( Offset > Len )
                    return -HUGE_VAL;

                #ifdef DEBUG
                    cerr << "   LogPathNum(" << Len << ", " << Offset << ") = "
                         << exp(lchoose( (Len + Offset)/2, Len )) << endl;
                #endif

                return lchoose( (Len + Offset)/2, Len );
            }


            double ComputePValue( int Neg, int Pos, int Score ){
                tTRecord R = {Neg, Pos, (double)Score} ;

                #ifdef DEBUG
                  cerr << "TNOM:ComputePValue(" << Neg << ", " << Pos << ", " ;
                  cerr << Score << ")" <<endl ;
                #endif

                int tc ;
                #pragma omp critical
                    tc = _pval_cache.count(R) ;

                if( tc > 0 )
                    return _pval_cache[R] ;

                // Need to reach A to get Score with positives on the left
                int A = Pos - Score ;
                // Need to reach -B to get Score with positives on the right
                int B = Neg - Score ;
                int L = Neg+Pos ;
                int O = Pos - Neg ;

                // if Neg or Pos = Score, every partition will give score Score or less
                // (things can only get better).
                if( A == 0 || B == 0 )
                    // Trivial partition does the job
                    return 1 ;

                // Use different accumulants to hold positive & negative counts
                double PNum = -HUGE_VAL ;
                double NNum = -HUGE_VAL ;

                int Ni = 2*B;
                int Pi = 2*A;
                int sign = 1;

                while( abs(Ni+O) <= L || abs(Pi-O) <= L ){
                    if( sign == 1 ){
                        PNum = AddLog(PNum, LogPathNum( L, Pi - O ));
                        PNum = AddLog(PNum, LogPathNum( L, -Ni - O ));
                    } else {
                        NNum = AddLog(NNum, LogPathNum( L, Pi - O ));
                        NNum = AddLog(NNum, LogPathNum( L, -Ni - O ));
                    }
                    sign *= -1 ;
                    int OldPi = Pi ;
                    int OldNi = Ni ;

                    // Recursion rule based on reflection principle
                    Pi = 2*A + OldNi ;
                    Ni = 2*B + OldPi ;

                    #ifdef DEBUG
                        cerr << " -> " << OldNi << " " << OldPi << " "
                         << exp(NNum) << " " << exp(PNum) << endl ;
                    #endif
                }

                // Normalize by total number of paths
                double NPath = LogPathNum(L, O) ;
                PNum -= NPath ;
                NNum -= NPath ;

                double PValue = exp(PNum) ;
                if( NNum > -HUGE_VAL )
                    PValue -= exp(NNum ) ;

                #ifdef DEBUG
                    cerr << "PValue( " << Neg << ", " << Pos << ", "
                         << Score << " ) -> " << exp(PNum) << " " << exp(NNum) << " "
                         << PValue << "\n" ;
                #endif
                #pragma omp critical
                    _pval_cache[R] = PValue ;
                return PValue ;
            }

        public:
            double Calc_pval(vector<float>& data, vector<int>& labels, float* score){
                // n = number of columms in the data set
                // BestLoss - minimum num of mistakes

                int n = data.size() ;
                double BestLoss = n ;

                int RightClass[MAXCLASS] ;
                int LeftClass[MAXCLASS] ;

                /**
                * first step - Right class is all the samples in the vector, Left class is empty.
                * count how many samples are in every class.
                * for example: there are 13 "1" and 4 "0" so RightClass[1]=13, RightClass[0]=4.
                **/
                for( int i = 0; i < MAXCLASS; i++ )
                    RightClass[i] = LeftClass[i] = 0;

                for( int i = 0; i < n; i++ ){
                    if(labels[i] >= 0)
                        RightClass[labels[i]]++;
                }

                /** mainloop
                * calculate number of "mistakes"
                * L - num of mistakes in right + left classes
                * X - expression level in i's sample (ordered), or hug_val if it's end of the vector.
                **/
                double LastValue = -HUGE_VAL;
                BestLoss = Loss(LeftClass) + Loss(RightClass) ;

                for( int i = 0; i <= n; i++ ){
                    if( i == n || labels[i] >= 0 ){


                        double L = Loss(LeftClass) + Loss(RightClass);
                        double X = (i < n) ? data[i] : HUGE_VAL;
                        if(L < BestLoss && X != LastValue) {
                            BestLoss = L;
                        }
                        // Move class
                        // Move one step to thr right, and change Right and Left class
                        if( i < n ) {
                            RightClass[labels[i]]--;
                            LeftClass[labels[i]]++;
                        }
                        LastValue = X;
                    }
                }
                *score = (float)BestLoss;
                return ComputePValue(LeftClass[0], LeftClass[1], (int)BestLoss) ;
            }
    };
}
#endif