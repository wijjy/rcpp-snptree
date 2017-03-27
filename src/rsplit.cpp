#include <Rcpp.h>
using namespace Rcpp;

#include "rcpp_binode.h"
#include "splitter.h"


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
SEXP  simple_split(IntegerMatrix d, IntegerVector positions) {
 
    splitter *s = new splitter(d, d.nrow(), d.ncol());            // define the splitter object s
    for (int i=0; i< positions.size(); i++) s->split(positions[i]);   // split at positions
  
 
  Rcpp::XPtr< splitter > pt(s, true);
  return pt;
}



