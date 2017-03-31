#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List testlist() 
{
   List L = List::create(
    1,
    2,
    seq_len(10)
  );
  L.push_back(7);
  for (int i=0;i<500;i++) L.push_back(seq_len(i));
  return L;
}

<<<<<<< HEAD
// [[Rcpp::export]]
NumericVector testvec() 
{
  stop("WTF");
  return NumericVector(10, 1.0);

}
=======
>>>>>>> 9f795b88ea059100451d9e9df8b3953bd71c0c88

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
a <- testlist()
<<<<<<< HEAD
testvec()
=======
>>>>>>> 9f795b88ea059100451d9e9df8b3953bd71c0c88
*/




