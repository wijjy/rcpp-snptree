#ifndef TESTSTATS_H__
#define TESTSTATS_H__

template <typename T>
double GStat(T x, T n, double p)
{
  double E0 = static_cast<double>(n)*p;
  double E1 = static_cast<double>(n)-E0;
  return 2.0*(static_cast<double>(x)*log(static_cast<double>(x)/E0) 
              + static_cast<double>(n-x)*log(static_cast<double>(n-x)/E1));
}  

template <typename T>
double AbsSevonStat(T x, T n, double p) {
  return fabs(static_cast<double>(x-static_cast<double>(n)*p)/sqrt(static_cast<double>(n)*p*(1.-p)));    
}


template <typename T>
double SqSevonStat(T x, T n, double p) {
  double z =  fabs(static_cast<double>(x-static_cast<double>(n)*p)/sqrt(static_cast<double>(n)*p*(1.-p)));    
  return z*z;
}

#ifdef USE_R
#include "Rmath.h"

template <typename T>
double LogP(T x, T n, double p)
{
  if (static_cast<double>(x)/static_cast<double>(n) < p)
    return -pbinom(static_cast<double>(x),p,static_cast<unsigned int>(n),1,1);
  else
    return -pbinom(static_cast<double>(x-1),p,static_cast<double>(n),1,1);
}  


template <typename T>
double LogPNorm(T x, T n, double p)
{
  double phat=static_cast<double>(x)/static_cast<double>(n) ;
  if (phat < p)
    return -pbinom(static_cast<double>(x),p,static_cast<unsigned int>(n),1,1)
    +pbinom(static_cast<double>(x),phat,static_cast<unsigned int>(n),1,1);
  else
    return -pbinom(static_cast<double>(x-1),p,static_cast<unsigned int>(n),0,1)
    +pbinom(static_cast<double>(x),phat,static_cast<unsigned int>(n),1,1);
}  

#else
#include "gsl/gsl_cdf.h"


template <typename T>
double LogP(T x, T n, double p)
{
  if (static_cast<double>(x)/static_cast<double>(n) < p)
    return -log(gsl_cdf_binomial_P(static_cast<unsigned int>(x),p,static_cast<unsigned int>(n)));
  else
    return -log(gsl_cdf_binomial_Q(static_cast<unsigned int>(x-1),p,static_cast<unsigned int>(n)));
}  


template <typename T>
double LogPNorm(T x, T n, double p)
{
  double phat=static_cast<double>(x)/static_cast<double>(n) ;
  if (phat < p)
    return -log(gsl_cdf_binomial_P(static_cast<unsigned int>(x),p,static_cast<unsigned int>(n)))
    +log(gsl_cdf_binomial_P(static_cast<unsigned int>(x),phat,static_cast<unsigned int>(n)));
  else
    return -log(gsl_cdf_binomial_Q(static_cast<unsigned int>(x-1),p,static_cast<unsigned int>(n)))
    +log(gsl_cdf_binomial_P(static_cast<unsigned int>(x),phat,static_cast<unsigned int>(n)));
}  



#endif

#endif
