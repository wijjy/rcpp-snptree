#ifndef UTIL_H_ijw
#define UTIL_H_ijw
#include <Rcpp.h>
#include <vector>

bool matches_hap(const Rcpp::IntegerVector &test,
                 const Rcpp::IntegerVector &target, 
                 const Rcpp::IntegerVector &target_positions );


bool SNPmatch( const Rcpp::IntegerMatrix &hap, const std::vector<int> &rows, int col);

template<typename T>
bool mismatch(std::vector<T> &labels, T *cases, int nc);

template<typename T>
bool mismatch2(std::vector<T> &labels, T *cases, int nc);

int count_intersection(std::vector<int> &a, const Rcpp::IntegerVector &b);
/** Note that both of these should be sorted                               */
int count_intersection(std::vector<int> &a, int *cases, int nc);

#endif