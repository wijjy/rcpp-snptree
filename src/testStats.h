#ifndef TESTSTATS_H__
#define TESTSTATS_H__


#include <Rcpp.h>
#include <vector>
#include "Rmath.h"
#include "R.h"

//using namespace Rcpp;

//two pass mean and variance calculation
template <typename T>
void avevar(const std::vector<T> &xx, double &mean, double &v) {
  mean=std::accumulate(xx.begin(),xx.end(),0.0)/static_cast<double>(xx.size());
  double ep=0.0,var=0.0;
  for (size_t ii=0;ii<xx.size();ii++) {
    double s=static_cast<double>(xx[ii])-mean;
    ep+=s;
    var+=s*s;
  }
  v = (var-ep*ep/static_cast<double>(xx.size()))/static_cast<double>(xx.size()-1);
  
}


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



template <typename T>
double LogP(T x, T n, double p)
{
  if (static_cast<double>(x)/static_cast<double>(n) < p)
    return -R::pbinom(static_cast<double>(x),p,static_cast<unsigned int>(n),1,1);
  else
    return -R::pbinom(static_cast<double>(x-1),p,static_cast<double>(n),1,1);
}  


template <typename T>
double LogPNorm(T x, T n, double p)
{
  double phat=static_cast<double>(x)/static_cast<double>(n) ;
  if (phat < p)
    return -R::pbinom(static_cast<double>(x),p,static_cast<unsigned int>(n),1,1)
    +R::pbinom(static_cast<double>(x),phat,static_cast<unsigned int>(n),1,1);
  else
    return -R::pbinom(static_cast<double>(x-1),p,static_cast<unsigned int>(n),0,1)
    +R::pbinom(static_cast<double>(x),phat,static_cast<unsigned int>(n),1,1);
}  

#endif
