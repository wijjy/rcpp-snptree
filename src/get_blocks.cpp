#include <Rcpp.h>
#include "rcpp_binode.h"
#include "splitter.h"
#include "ijw_rand.h" 

/** Return a bifurcation diagram as a set of blocks, each block being a node 
 * in the bifurcation tree.  This doesn't work for log yet, but I can try that later.  Far 
 * more imprtant to work on the movement and selection
 */

// [[Rcpp::export]]
Rcpp::NumericMatrix get_blocks(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);
  if (s->root()->range.first==0) 
    Rcpp::stop("You have not calculated the positions of the leaves.");
  
  int leaves = s->nleaves();
  Rcpp::NumericMatrix boxes(2*leaves-2, 5);  // last one left for the root if needed
  
  int index=0;
  NLRIterator<binode> ii(s->root());
  while (!ii.isend()) {
    if ((*ii)->isleaf())
      Rcpp::stop("should never get to a leaf in this function");
    // left
    boxes(index, 0) = (*ii)->position+1;
    boxes(index, 1) = (*ii)->left->position+1;
    boxes(index, 2) = (*ii)->range.first;
    boxes(index, 3) = (*ii)->left->range.first;
    boxes(index, 4) = (*ii)->left->height();
    // right
    index++;
    boxes(index, 0) = (*ii)->position+1;
    boxes(index, 1) = (*ii)->right->position+1;
    boxes(index, 2) = (*ii)->range.first+(*ii)->left->height();
    boxes(index, 3) = (*ii)->right->range.first;
    boxes(index, 4) = (*ii)->right->height();
    ii.nextInternal();
    index++;
  }
  return boxes;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix get_id_blocks(SEXP ptr, Rcpp::IntegerVector id) {
  Rcpp::XPtr< splitter > s(ptr);
 // if (s->root()->range.first==0)
 //   s->calculate_id_top_bottom(id);   // gets the tops and bottoms 
  
  int leaves = s->nleaves();
  Rcpp::NumericMatrix boxes(2*leaves-2, 5);  // last one left for the root if needed
  
  int index=0;
  NLRIterator<binode> ii(s->root());
  while (!ii.isend()) {
    if ((*ii)->isleaf())
      Rcpp::stop("should never get to a leaf in this function");
    // left
    boxes(index, 0) = (*ii)->position+1;
    boxes(index, 1) = (*ii)->left->position+1;
    boxes(index, 2) = (*ii)->range.first;
    boxes(index, 3) = (*ii)->left->range.first;
    boxes(index, 4) = (*ii)->left->height();
    // right
    index++;
    boxes(index, 0) = (*ii)->position+1;
    boxes(index, 1) = (*ii)->right->position+1;
    boxes(index, 2) = (*ii)->range.first+(*ii)->left->height();
    boxes(index, 3) = (*ii)->right->range.first;
    boxes(index, 4) = (*ii)->right->height();
    ii.nextInternal();
    index++;
  }
  return boxes;
}


/*** R
library(rcppsnptree)
data(snptreeExample)
id <- sort(sample(nrow(haps), 1000))
sp <- simple_split(haps, 1:24)
set_leaf_position(sp, 24)
calc_node_ranges(sp,  gap=40)

blocks<- get_blocks(sp)
id_blocks <- get_id_blocks(sp, id = id)
## little plot
plot(-1,-1, xlim=range(blocks[, 1:2]), 
     ylim=range(c(blocks[, 3:4], blocks[, 3:4]+blocks[, 5])), 
     type="n", axes=FALSE, xlab="", ylab="")
apply(blocks, 1, plot_block)
apply(id_blocks, 1, plot_block, col="lightblue")
*/
