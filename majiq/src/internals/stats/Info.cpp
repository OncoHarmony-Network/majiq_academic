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

            /**
             * for each part of the vector compute Entropy = sum{ numOfSamplesInClass * log(SumOfSamplesInVector/numOfSamplesInClass) }
             * @param N - pointer to int array which hold for every calss the number of samples in it.
             * @return the Entropy calculated.
             **/
            double Entropy(int *N) const {
                double Sum = 0 ;
                double E = 0.0 ;

                for( int i = 1; i < ClassNum(); i++ ){
                    Sum += N[i] ;
                }
                for( int i = 1; i < ClassNum(); i++ ){
                    if( N[i] > 0 ){
                        E -= N[i] * log( N[i]/Sum ) ;
                    }
                }
                return E;
            }

           /**
            * main method for computing PVal
             * @param Neg - number of negative labels
             * @param Pos - number of positive labels
             * @param Score - threshold used for computing Pval. (PVal = Prob(g,U(n,p)) =< Score)
             * @return Pvalue with the given parameters.
             **/
            double ComputePValue( int Neg, int Pos, double Score ){
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
                return PVal ;
            }


        public:
            double Calc_pval(vector<float> data, vector<int> labels){

                int n = data.size() ;
                double BestLoss = n ;

                int RightClass[MAXCLASS] ;
                int LeftClass[MAXCLASS] ;
                float _Threshold = 0

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

                        double L = Entropy(LeftClass) + Entropy(RightClass) ;
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
                return PValue( LeftClass[1], LeftClass[2], BestLoss) ;
            }
}
#endif






















/************************
File: Info.cpp
Description: class tInfoScore calculate for a given gene g, Score and P-Value, using the
             Mutual Information algorithm. the method is similiar to TNOM, but instead
	     "counting mistakes" it measures the entropy for every  decision rule.
**************************/
#include <MainLib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#include <ExpressionData.h>
#include <Info.h>
#include <MathFunctions.h>

tPValueCache tInfoScore::_PValueCache;

/**
 * Ctor
 * @param Data - reference to the data base
 * @param g - number of gene to be scored
 * @param Labels - reference to the vector that hold tha samples labels
 **/
tInfoScore::tInfoScore(tData const& Data, int g, tLabels const& Labels)
  : tRankScore( Data, g, Labels )
{
}
// Dtor 
tInfoScore::~tInfoScore()
{}

double
tInfoScore::ComputePathScore(int Neg, int Pos, int k, int O) const
{
  int Left[3];
  int Right[3];

  if ( ClassNum() != 2+1 )
    Data().ExitWithError(" Number of classes in classification file is illegal. Check there are no extra \'Enters\' at the end of Labels file. ");
  
  Left[1] = (k - O)/2;
  Left[2] = (k + O)/2;
  Right[1] = Neg - Left[1];
  Right[2] = Pos - Left[2];

  return Entropy(Left) + Entropy(Right);
}

/**
 * Pval: first look at the cache, if Pval was not updated there, compute it, and insert to cache.
 * @param Neg - number of negative labels
 * @param Pos - number of positive labels
 * @param Score - threshold used for computing Pval. (PVal = Prob(g,U(n,p)) =< Score)
 * @return Pvalue with the given parameters.
 **/
double
tInfoScore::PValue( int Neg, int Pos, double Score )
{
  pair<bool,double> i = _PValueCache.Retrieve(Neg, Pos, Score);
  if( i.first )
  {
    //    double p = ComputePValue(Neg, Pos, Score );
    //    if( p != i.second )
    //    cerr << "PValue("<< Neg << ", " << Pos << ", " << Score << ") = " << p << " " << i.second << "\n";
    return i.second;
  }

  double p = ComputePValue(Neg, Pos, Score );
  _PValueCache.Insert(Neg, Pos, Score, p);

  return p;
}

/**
 * main method for computing PVal
 * @param Neg - number of negative labels
 * @param Pos - number of positive labels
 * @param Score - threshold used for computing Pval. (PVal = Prob(g,U(n,p)) =< Score)
 * @return Pvalue with the given parameters.
 **/
double
tInfoScore::ComputePValue( int Neg, int Pos, double Score )
{
  int L = Pos+Neg;
  int O = Pos-Neg;
  vector<double> M1(L+1);
  vector<double> M2(L+1);

  //  cerr << "InfoScore::ComputePValue(" << Neg << ", " << Pos << ", " << Score << ")" << endl;
  
  vector<double>* OldCounts = &M1;
  vector<double>* NewCounts = &M2;
  int k = 1;
  int LastMax = 0;
  int LastMin = 0;
  int i;
  for( i = 0; i < L+1; i++ )
    (*OldCounts)[i] = -HUGE_VAL;
  (*OldCounts)[Neg] = 0.0;

  // An alternative count of "Bad" paths that achieve the score
  double BadPaths = -HUGE_VAL;
  
  for( k = 1; k <= L; k++ )
  {
    int M = LastMax+1;
    int m = LastMin-1;
    if( M > Pos )
      M = Pos;
    if( M > O + L-k )
      M = O + L-k;
    if( m < -Neg )
      m = -Neg;
    if( m < O - (L-k) )
      m = O - (L-k);

    int NewMin = Pos;
    int NewMax = -Neg;

    //    cerr << " k = " << k << " " << " BadPaths = " << exp(BadPaths) << "\n";
    double S = 0;
    for( i = m; i <= M; i++ )
    {
      (*NewCounts)[i+Neg] = -HUGE_VAL;

      double C = -HUGE_VAL;
      if( (i - 1 >= LastMin) && (i -1 <= LastMax) )
	C = AddLog(C, (*OldCounts)[i-1+Neg]);
      if( (i + 1 >= LastMin) && (i + 1 <= LastMax) )
	C = AddLog(C, (*OldCounts)[i+1+Neg]);

      if( C > -HUGE_VAL )
      {
	S = ComputePathScore(Neg, Pos, k, i);

	// this assert is crucial to aviod optimization bug!!!
	assert( S >= 0.0 );
	
	//	cerr << "   C[" << i << "] = ";
	//	cerr << exp(C) << " ";
	//	cerr << " S = " <<  S ;

	if( S > Score )
	{
	  //	  cerr << " good";
	  (*NewCounts)[i+Neg] = C;

	  if( i < NewMin )
	    NewMin = i;
	  if( i > NewMax )
	    NewMax = i;
	} else
	{
	  // We hit the boundry of allowable region.
	  // Every path from here to the end is bad
	  int dx = L - k;
	  int dy = O - i;
	  int p = (dx + dy)/2;
	  assert( (dy >= 0 ) ? (dy <= dx) : (-dy <= dx) );
	  assert( p >= 0 && p <= dx );

	  //cerr << "   C[" << i << "] = ";
	  //cerr << exp(C) << " ";
	  //cerr << " S = " <<  S ;

	  //
	  // p = number of positive steps we need to make
	  C += lchoose( p, dx );
	  //cerr << " bad " << exp(lchoose( p, dx )) << " Total = " << exp(C);
	  //cerr << endl;
	  BadPaths = AddLog( BadPaths, C );
	}
	//	cerr << endl;
      }
    }
    if( NewMin > NewMax )
      return 1.0;
    
    LastMin = NewMin;
    LastMax = NewMax;
    vector<double>* Temp = OldCounts;
    OldCounts = NewCounts;
    NewCounts = Temp;
  }

#ifdef DEBUG
  cerr << "Paths = " << exp( (*OldCounts)[O+Neg]) << " of "
       << exp(lchoose(Neg,Neg+Pos)) << " "
       << 1 - exp( (*OldCounts)[O+Neg] - lchoose(Neg,Neg+Pos)) << " "
       << (exp(lchoose(Neg,Neg+Pos)) - exp( (*OldCounts)[O+Neg]))/exp(lchoose(Neg,Neg+Pos)) << "\n";
  cerr << "BadPaths = " << exp( BadPaths ) << " of "
       << exp(lchoose(Neg,Neg+Pos)) << " "
       << exp( BadPaths - lchoose(Neg,Neg+Pos)) << "\n";
#endif
  
  double PVal = exp( BadPaths - lchoose(Neg,Neg+Pos));
  return PVal;
}

/**
 * for each part of the vector compute Entropy = sum{ numOfSamplesInClass * log(SumOfSamplesInVector/numOfSamplesInClass) }
 * @param N - pointer to int array which hold for every calss the number of samples in it.
 * @return the Entropy calculated.
 **/
double
tInfoScore::Entropy(int *N) const
{
  double Sum = 0;
  double E = 0.0;
  int i;

  for( i = 1; i < ClassNum(); i++ )
    Sum += N[i];
  for( i = 1; i < ClassNum(); i++ )
    if( N[i] > 0 )
      E -= N[i] * log( N[i]/Sum );
  
  return E;
}

/**
 * compute Info score for the associated gene.
 * this is an iterative method that check for every possible distribution the Entropy,
 * and choose the distribution with the minimal Entropy.
 * Right and Left class hold the number of samples from each class.
 * @return PVal of the score that was found
 **/
double
tInfoScore::operator()()
{
   // n = number of columms in the data set
   // BestLoss - minimum num of mistakes  
  int i;
  int n = Data().NumberCond();
  double BestLoss = n;

  int RightClass[MAXCLASS];
  int LeftClass[MAXCLASS];
    
  for( i = 0; i < ClassNum(); i++ )
    RightClass[i] = LeftClass[i] = 0;

  /**
   * first step - Right class is all the samples in the vector, Left class is empty.
   * count how many samples are in every class.
   **/
  for( i = 0; i < n; i++ )
    if(  Data().Obs(Gene(), i) )
      RightClass[Label(i)]++;

  /** main loop
   * calculate minimal Entropy
   * L - Entropy in right + left classes
   * X - expression level in i's sample (ordered), or hug_val if it's end of the vector.
   *
   * _Threshold hold the avarage expression level in the place we are standing
   **/
  double LastValue = -HUGE_VAL;
  for(i = 0; i <= n; i++ )
    if( i == n || Data().Obs( Gene(), Order()[i] ) )
    {
      double L = Entropy(LeftClass) + Entropy(RightClass);
      double X = (i < n) ? Data().Level(Gene(),Order()[i]) : HUGE_VAL;
      if( i == 0 || (L < BestLoss && X != LastValue )  )
      {
	BestLoss = L;
	if( i == 0 )
	  _Threshold = Data().Level(Gene(),Order()[i]) - 1;
	else
	  if( i < n )
	    _Threshold = (Data().Level(Gene(),Order()[i]) + LastValue )/2;
	  else
	    _Threshold = Data().Level(Gene(),Order()[n-1]) + 1;
      }
	
      // Move class
      if( i < n )
      {
	RightClass[Label(Order()[i])]--;
	LeftClass[Label(Order()[i])]++;
      }
      LastValue = X;
    }
  //  return BestLoss/(LeftClass[1]+LeftClass[2]);
  if ( ClassNum() !=  2+1 )
    Data().ExitWithError(" Number of classes in classification file is illegal. Check there are no extra \'Enters\' at the end of Labels file. ");
  //  assert( ClassNum() == 2+1 );
  return PValue( LeftClass[1], LeftClass[2], BestLoss);
}

/**
 * insert the threshold value to buffer writer.
 * @param o - reference to the buffer writer
 **/
void
tInfoScore::Print(ostream& o) const
{
  o << _Threshold;
}

// @return true if the object can copute PVal.
bool
tInfoScore::ComputesPVal() const
{
  return true;
}
