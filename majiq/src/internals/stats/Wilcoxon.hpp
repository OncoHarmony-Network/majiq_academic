#ifndef WILCOXON_H
#define WILCOXON_H
#include "MathFunctions.h"
#include "testStats.hpp"
#include <random>
#include <algorithm>
#include <string>
#include <assert.h>
#include <map>

#define MAXCLASS 2

using namespace std ;
namespace stats{

    class Wilcoxon: public TestStat{

        private:

            struct tWilcoxonCountRecord{
                int n ;
                int k ;
                int s ;

                bool operator<(const tWilcoxonCountRecord& rhs) const{
                        return (n < rhs.n || (n >= rhs.n && k >= rhs.k) ||
                                (n >= rhs.n && k >= rhs.k && s < rhs.s)) ;
                }


            };

             map< tWilcoxonCountRecord, double, less<tWilcoxonCountRecord >> _pval_cache ;

            double Count(int n, int k, int s ){
                if( s < 0 )
                    return -HUGE_VAL ;
                if( k == 0 ){
                    if( s > 0 )
                        return -HUGE_VAL ;
                    else
                        return lgamma(n+1) ;
                }

                if( k == n ){
                    if( s == (k*(k+1))/2 )
                        return  lgamma(n+1) ;
                    else
                        return -HUGE_VAL ;
                }

                tWilcoxonCountRecord R = {n,k,s} ;
                map<tWilcoxonCountRecord, double, less<tWilcoxonCountRecord> >::iterator i = _pval_cache.find(R) ;

                if( i != _pval_cache.end() )
                    return (*i).second ;

                assert( k > 0 ) ;
                double c1 = Count( n -1 , k, s ) ;
                //  cout << "Count(" <<  n-1 << "," << k << "," << s << ") = " << c1 << " " << exp(c1)  << "\n";

                double c2 = Count( n -1 , k-1, s - n ) ;
                //  cout << "Count(" <<  n-1 << "," << k-1 << "," << s-n << ") = " << c2 << " " << exp(c2)  << "\n";

                c1 += log( 1.0*(n - k) ) ;
                c2 += log( 1.0*k ) ;
                double c = AddLog(c1, c2) ;

                map<tWilcoxonCountRecord, double, less<tWilcoxonCountRecord> >::value_type v(R,c) ;
                _pval_cache.insert( v ) ;

                return c ;
            }

        public:
            double Calc_pval(vector<float> data, vector<int> labels){
                int n1 = 0 ;
                int n2 = 0 ;
                double s = 0 ;

                int nsampls = data.size() ;

                int j = 0 ;
                int i = 0 ;
                while( i < nsampls ){
                    if( labels[i] >= 0 ) {
                        double x = data[i] ;
                        int m1 = 0 ;
                        int m2 = 0 ;
                        // as long as expression level doesn't change, count how many samples labeled "1" (m1)
                        // and how many labeled "2" (m2).
                        for( ; i < nsampls && data[i]; i++ )
                            if( labels[i] >= 0 ){
                                if( labels[i] == 0 )
                                    m1++ ;
                                if( labels[i] == 1 )
                                    m2++ ;
                            }
                            if( m1+m2 > 0 ){
                                if( m1 > 0 )
                                    s += m1*(j + 0.5*(m1+m2+1)) ;
                                j += m1 + m2 ;
                                n1 += m1 ;
                                n2 += m2 ;
                            }
                    }
                }
                double _ZScore = (s - 0.5*n1*(n1+n2+1))/sqrt(n1*n2*(n1+n2+1)/12.0);

                double PValue = -HUGE_VAL ;
                int s1;
                int s0 = (int)floor(s+0.5) ; // added int cast
                if( _ZScore < 0 ){
                    for( s1 = (n1*(n1+1))/2; s1 <= s0; s1++ ){
                        PValue = AddLog(PValue, Count(n1+n2,n1,s1)) ;
                    }
                }else{
                    for( s1 = ((n1+n2)*(n1+n2+1) - n2*(n2-1))/2; s1 >= s0; s1-- ){
                        PValue = AddLog(PValue, Count(n1+n2,n1,s1)) ;
                    }
                }
                PValue -= lgamma(n1+n2 + 1) ;

                //  cerr << "Wilcoxon PValue (" << n1+n2 << ", " << n1 << "," << s <<") = " << exp(PValue)*2 << " compare to " << 2*GaussCDF(-fabs(_ZScore), 0, 1) << "\n";

                return exp(PValue)*2 ;
            }
    } ;
}
#endif