#include <Rcpp.h>


bool matches_hap(const Rcpp::IntegerVector &test,
                 const Rcpp::IntegerVector &target, 
                 const Rcpp::IntegerVector &target_positions ) {
  for (int ii=0; ii<target_positions.size(); ii++)  {
    if (target[target_positions[ii]] != -1)
      if (test(target_positions[ii]) != target(target_positions[ii]))
        return false;
  }
  return true;
}


/** Do all the SNPs at position col for rows rows match?                         */
bool SNPmatch(const Rcpp::IntegerMatrix &hap, const std::vector<int> &rows, int col) 
{
  int SNP = hap(rows[0], col);
  for (size_t j=0; j<rows.size(); j++) {
    if (hap.at(rows.at(j), col) != SNP) return false;      // BOUNDS CHECK
  }
  return true;
} 
/** Are any of the labels not cases?                                             */
template<typename T>
bool mismatch(std::vector<T> &labels, T *cases, int nc)
{
  typename std::vector<T>::iterator curr=labels.begin();
  typename std::vector<T>::iterator end=labels.end();
  typename std::vector<T>::iterator fo=std::find_first_of(curr,end,cases,cases+nc);
  if (fo==end) return false; // no labels in cases - all controls
  if (fo!=curr) return true; // first label not a case
  curr++;
  
  while (curr != end) {
    // find the first label (from curr) that is in the cases
    std::vector<int>::iterator fo=std::find_first_of(curr,end,cases,cases+nc);
    if (fo!=curr)  //the first label is not  a case
      return true;
    curr++;
  }
  //curr==end so all labels must be in cases!
  return false;
}

/** Note that both of these should be sorted                               */
int count_intersection(std::vector<int> &a, const Rcpp::IntegerVector &b)
{
  std::vector<int>::iterator al=a.end(),ii=a.begin();
  Rcpp::IntegerVector::const_iterator bl=b.end(),jj=b.begin();
  int count=0;
  while (ii!=al && jj!=bl) {
    if (*ii<*jj) ++ii;
    else if (*jj<*ii) ++jj;
    else {
      count++;
      ii++;
      jj++;
    }
  }
  return count;
} 
/** Note that both of these should be sorted                               */
int count_intersection(std::vector<int> &a, int *cases, int nc)
{
  std::vector<int>::iterator al=a.end(),ii=a.begin();
  int *jj=cases,*bl=cases+nc;
  int count=0;
  while (ii!=al && jj!=bl) {
    if (*ii<*jj) ++ii;
    else if (*jj<*ii) ++jj;
    else {
      count++;
      ii++;
      jj++;
    }
  }
  return count;
}