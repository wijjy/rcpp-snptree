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
    splitter *s = new splitter(d);                                     // define the splitter object s
    for (int i=0; i<positions.size(); i++) s->split(positions[i]-1);   // split at positions
    Rcpp::XPtr< splitter > pt(s, true);                                // get pointer as SEXP
    return pt;
}


// [[Rcpp::export]]
int  nleaves(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);
  return s->nleaves();
}

// [[Rcpp::export]]
IntegerMatrix  case_control_leaves(SEXP ptr, IntegerVector cases) {
  Rcpp::XPtr< splitter > s(ptr);
  IntegerMatrix xxx(s->nleaves(), 2);
  s->getCaseControlLeaves(xxx, cases);
  return xxx;
}

// [[Rcpp::export]]
IntegerVector  leaf_count(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);
  IntegerVector counts(s->nleaves());
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  int index=0;
  while (!ii.isend()) {
    counts[index++] = (*ii)->labels.size();
    ii.nextLeaf();
  }
  return counts;
}

// [[Rcpp::export]]
NumericMatrix  leaf_positions(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  NumericMatrix ypositions(leaves);
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  ypositions(0, 0) = 0;
  ypositions(0, 1) = (*ii)->labels.size();
  ii.nextLeaf();
  int index=1;
  while (!ii.isend()) {
    ypositions(index, 0) = gap + ypositions(index-1, 1);
    ypositions(index, 1) = ypositions(index, 0) + (*ii)->labels.size();
    index++;
    ii.nextLeaf();
  }
  return ypositions;
}

// [[Rcpp::export]]
NumericMatrix  node_positions(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  NumericMatrix ypositions(leaves);
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  ypositions(0, 0) = 0;
  ypositions(0, 1) = (*ii)->labels.size();
  ii.nextLeaf();
  int index=1;
  while (!ii.isend()) {
    ypositions(index, 0) = gap + ypositions(index-1, 1);
    ypositions(index, 1) = ypositions(index, 0) + (*ii)->labels.size();
    index++;
    ii.nextLeaf();
  }
  return ypositions;
}


// [[Rcpp::export]]
NumericMatrix  stumps(SEXP ptr, int pos1, int pos2, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  NumericMatrix ypositions(leaves);
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  ypositions(0, 0) = 0;
  ypositions(0, 1) = (*ii)->labels.size();
  ii.nextLeaf();
  int index=1;
  while (!ii.isend()) {
    ypositions(index, 0) = gap + ypositions(index-1, 1);
    ypositions(index, 1) = ypositions(index, 0) + (*ii)->labels.size();
    index++;
    ii.nextLeaf();
  }
  return ypositions;
}

// [[Rcpp::export]]
NumericMatrix  get_coordinates(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  s->calculate_top_bottom(gap);   // gets the tops and bottoms for the and the rest of the tree
  
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  while (!ii.isend()) {
    Rprintf("range for leaf = (%g, %g)\n", (*ii)->range.first, (*ii)->range.second);
    ii.nextLeaf();
  }
  // have the tops and bottoms for all nobes, now just have to take the coordinates in order.
  NumericMatrix coords(2*leaves+3*(leaves-1), 3);
  
  s->get_coordinates(coords);
  
  return coords;
}




