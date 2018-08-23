#ifndef INFOSCORE_H
#define INFOSCORE_H
#include "MathFunctions.h"
#include "testStats.hpp"
#include <random>
#include <algorithm>
#include <string>
#include <assert.h>
#include <map>

#define MAXCLASS 2

using namespace std ;
namespace MajiqStats{

    class InfoScore: public MajiqStats::TestStat{

        private:

            struct tInfoRecord{
                int neg ;
                int pos ;
                double score ;

                bool operator<(const tInfoRecord& rhs) const{
                        return (neg < rhs.neg || (neg >= rhs.neg && pos >= rhs.pos) ||
                                (neg >= rhs.neg && pos >= rhs.pos && score < rhs.score)) ;
                }

            };

            map<tInfoRecord, double> _pval_cache ;

            /**
             * for each part of the vector compute Entropy = sum{ numOfSamplesInClass * log(SumOfSamplesInVector/numOfSamplesInClass) }
             * @param N - pointer to int array which hold for every calss the number of samples in it.
             * @return the Entropy calculated.
             **/
            double Entropy(int *N, int M) const {
                double Sum = 0 ;
                double E = 0.0 ;

                for( int i = 0; i < M; i++ ){
                    Sum += N[i] ;
                }
                for( int i = 0; i < M; i++ ){
                    if( N[i] > 0 ){
                        E -= N[i] * log( N[i]/Sum ) ;
                    }
                }
                return E;
            }

            double ComputePathScore(int Neg, int Pos, int k, int O) const {
                int Left[3] ;
                int Right[3] ;

                Left[1] = (k - O)/2 ;
                Left[2] = (k + O)/2 ;
                Right[1] = Neg - Left[1] ;
                Right[2] = Pos - Left[2] ;

                return Entropy(Left, 2) + Entropy(Right, 2) ;
            }


           /**
            * main method for computing PVal
             * @param Neg - number of negative labels
             * @param Pos - number of positive labels
             * @param Score - threshold used for computing Pval. (PVal = Prob(g,U(n,p)) =< Score)
             * @return Pvalue with the given parameters.
             **/
            double ComputePValue( int Neg, int Pos, double Score ){

                tInfoRecord R = {Neg, Pos, Score} ;
                if (_pval_cache.count(R)){
                    return _pval_cache[R] ;
                }


                int L = Pos+Neg ;
                int O = Pos-Neg ;
                vector<double> M1(L+1) ;
                vector<double> M2(L+1) ;

                //  cerr << "InfoScore::ComputePValue(" << Neg << ", " << Pos << ", " << Score << ")" << endl;

                vector<double>* OldCounts = &M1 ;
                vector<double>* NewCounts = &M2 ;
//                int k = 1 ;
                int LastMax = 0 ;
                int LastMin = 0 ;
                for( int i = 0; i < L+1; i++ ){
                    (*OldCounts)[i] = -HUGE_VAL;
                }
                (*OldCounts)[Neg] = 0.0;

                // An alternative count of "Bad" paths that achieve the score
                double BadPaths = -HUGE_VAL;

                for( int k = 1; k <= L; k++ ) {
                    int M = LastMax+1 ;
                    int m = LastMin-1 ;
                    if( M > Pos )
                        M = Pos ;
                    if( M > O + L-k )
                        M = O + L-k ;
                    if( m < -Neg )
                        m = -Neg ;
                    if( m < O - (L-k) )
                        m = O - (L-k) ;

                    int NewMin = Pos ;
                    int NewMax = -Neg ;

                //    cerr << " k = " << k << " " << " BadPaths = " << exp(BadPaths) << "\n";
                    double S = 0 ;
                    for( int i = m; i <= M; i++ ) {
                        (*NewCounts)[i+Neg] = -HUGE_VAL ;

                        double C = -HUGE_VAL ;
                        if( (i - 1 >= LastMin) && (i -1 <= LastMax) )
                            C = AddLog(C, (*OldCounts)[i-1+Neg]) ;
                        if( (i + 1 >= LastMin) && (i + 1 <= LastMax) )
                            C = AddLog(C, (*OldCounts)[i+1+Neg]) ;

                        if( C > -HUGE_VAL ) {
                            S = ComputePathScore(Neg, Pos, k, i) ;

                            // this assert is crucial to aviod optimization bug!!!
                            assert( S >= 0.0 ) ;

                            //	cerr << "   C[" << i << "] = ";
                            //	cerr << exp(C) << " ";
                            //	cerr << " S = " <<  S ;

                            if( S > Score ) {
                            //	  cerr << " good";
                                (*NewCounts)[i+Neg] = C ;

                                if( i < NewMin )
                                    NewMin = i ;
                                if( i > NewMax )
                                    NewMax = i ;
                            } else {
                                // We hit the boundry of allowable region.
                                // Every path from here to the end is bad
                                int dx = L - k ;
                                int dy = O - i ;
                                int p = (dx + dy)/2;
                                assert( (dy >= 0 ) ? (dy <= dx) : (-dy <= dx) ) ;
                                assert( p >= 0 && p <= dx ) ;

                                //cerr << "   C[" << i << "] = ";
                                //cerr << exp(C) << " ";
                                //cerr << " S = " <<  S ;

                                //
                                // p = number of positive steps we need to make
                                C += lchoose( p, dx ) ;
                                //cerr << " bad " << exp(lchoose( p, dx )) << " Total = " << exp(C);
                                //cerr << endl;
                                BadPaths = AddLog( BadPaths, C ) ;
                            }
                            //	cerr << endl;
                        }
                    }
                    if( NewMin > NewMax )
                      return 1.0 ;

                    LastMin = NewMin ;
                    LastMax = NewMax ;
                    vector<double>* Temp = OldCounts ;
                    OldCounts = NewCounts ;
                    NewCounts = Temp ;
                }

                #ifdef DEBUG
                  cerr << "Paths = " << exp( (*OldCounts)[O+Neg]) << " of "
                       << exp(lchoose(Neg,Neg+Pos)) << " "
                       << 1 - exp( (*OldCounts)[O+Neg] - lchoose(Neg,Neg+Pos)) << " "
                       << (exp(lchoose(Neg,Neg+Pos)) - exp( (*OldCounts)[O+Neg]))/exp(lchoose(Neg,Neg+Pos)) << "\n" ;
                  cerr << "BadPaths = " << exp( BadPaths ) << " of "
                       << exp(lchoose(Neg,Neg+Pos)) << " "
                       << exp( BadPaths - lchoose(Neg,Neg+Pos)) << "\n" ;
                #endif

                double PVal = exp( BadPaths - lchoose(Neg,Neg+Pos)) ;
                _pval_cache[R] = PVal ;
                return PVal ;
            }


        public:
            double Calc_pval(vector<float> data, vector<int> labels){

                int n = data.size() ;
                double BestLoss = n ;

                int RightClass[MAXCLASS] ;
                int LeftClass[MAXCLASS] ;
                float _Threshold = 0 ;

                /**
                * first step - Right class is all the samples in the vector, Left class is empty.
                * count how many samples are in every class.
                * for example: there are 13 "1" and 4 "0" so RightClass[1]=13, RightClass[0]=4.
                **/
                for( int i = 0; i < MAXCLASS; i++ )
                    RightClass[i] = LeftClass[i] = 0 ;

                for( int i = 0; i < n; i++ ){
                    if(labels[i] <= 0)
                        RightClass[labels[i]]++ ;
                }


                /** main loop
                * calculate minimal Entropy
                * L - Entropy in right + left classes
                * X - expression level in i's sample (ordered), or hug_val if it's end of the vector.
                *
                * _Threshold hold the avarage expression level in the place we are standing
                **/
                double LastValue = -HUGE_VAL ;
                for(int i = 0; i <= n; i++ ) {

                    if( i == n || labels[i] <= 0 ){

                        double L = Entropy(LeftClass, MAXCLASS) + Entropy(RightClass, MAXCLASS) ;
                        double X = (i < n) ? data[i] : HUGE_VAL ;

                        if( i == 0 || (L < BestLoss && X != LastValue )  ) {
                            BestLoss = L ;
                            if( i == 0 ) {
                                _Threshold = data[i] - 1 ;
                            } else {
                                if( i < n ) {
                                    _Threshold = (data[i] + LastValue )/2 ;
                                } else {
                                    _Threshold = data[n-1] + 1 ;
                                }
                            }
                        }

                        // Move class
                        if( i < n ) {
                            RightClass[labels[i]] -- ;
                            LeftClass[labels[i]] ++ ;
                        }
                        LastValue = X ;
                    }
                }
                return ComputePValue( LeftClass[0], LeftClass[1], BestLoss) ;
            }

    };
}
#endif
