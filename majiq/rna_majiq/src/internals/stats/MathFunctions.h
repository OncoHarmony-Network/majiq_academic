#ifndef __MathFunctions_h
#define __MathFunctions_h

#include <math.h>

inline
double lchoose( int k, int N )
{
  return lgamma(N+1) - lgamma(k+1) - lgamma(N-k+1);
}

inline double l2gamma(double x )
{
  return lgamma(x) / log(2.0);
}

inline double digamma(double x)
{
  return ( log(x) - 1/(2*x) - 1/(12*pow(x,2)) +
	   1/(120*pow(x,4)) - 1/(252*pow(x,6)));
}

inline double digamma1(double x)
{
  return (1/x + 1/(2*pow(x,2) + 1/(6*pow(x,3)) -
	   1/(30*pow(x,5)) + 1/(42*pow(x,7)) ));
}

inline
double
AddLog(double x, double y )
{
  if( x == -HUGE_VAL ) return y;
  if( y == -HUGE_VAL ) return x;
  if( x >= y )
    y -= x;
  else
  {
    double t = x;
    x = y;
    y = t - x;
  }
  return x + log(1 + exp(y));
}


inline double
GaussCDF(double x, double mu, double sigma2 )
{
  double z = (x -mu)/sqrt(sigma2);
  return 0.5*( 1 + erf(z / M_SQRT2 ));
}

inline double
GaussPDF(double x, double mu, double sigma2 )
{
  double z = (x -mu);
  return 1/sqrt(2*M_PI*sigma2)*exp(-0.5*z*z/sigma2);
}

inline double
LogGaussPDF(double x, double mu, double sigma2 )
{
  double z = (x -mu);
  return -0.5*(z*z/sigma2 + log(2*M_PI*sigma2));
}

inline double
Log2(double x)
{
  return log(x)/log(2.0);
}
inline double
Exp2(double x)
{
  return exp(x * log(2.0));
}

double LogBinomialTail(double p, int k, int N);
void cumbin(double*,double*,double*,double*,double*,double*);

#endif
