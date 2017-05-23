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
    boxes(index, 2) = (*ii)->range.first + (*ii)->left->height();
    boxes(index, 3) = (*ii)->right->range.first;
    boxes(index, 4) = (*ii)->right->height();
    ii.nextInternal();
    index++;
  }
  return boxes;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix get_blocks3(SEXP ptr) {
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
    boxes(index, 2) = (*ii)->range.first + (*ii)->left->height();
    boxes(index, 3) = (*ii)->right->range.first;
    boxes(index, 4) = (*ii)->right->height();
    ii.nextInternal();
    index++;
  }
  return boxes;
}
/** Return a matrix of blocks which have the same order as the nodes
  * from APE           */
// [[Rcpp::export]]
Rcpp::NumericMatrix get_id_blocks(SEXP ptr, Rcpp::IntegerVector id) {
  Rcpp::XPtr< splitter > s(ptr);
  
  s->calculate_id_top_bottom(id);   // gets the tops and bottoms 
  
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
    boxes(index, 2) = (*ii)->id_range.first;
    boxes(index, 3) = (*ii)->left->id_range.first;
    boxes(index, 4) = (*ii)->left->id_height();
    // right
    index++;
    boxes(index, 0) = (*ii)->position+1;
    boxes(index, 1) = (*ii)->right->position+1;
    boxes(index, 2) = (*ii)->id_range.first+(*ii)->left->id_height();
    boxes(index, 3) = (*ii)->right->id_range.first;
    boxes(index, 4) = (*ii)->right->id_height();
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
apply(id_blocks, 1, plot_block, col="red")


library(ARG)
library(rcppsnptree)
opar <- par(mfrow=c(2,2), mar=c(1,1,1,1))
test_plot <- function(ss, rec, var, gap=10,...) {

  b <- simARG(ss, sites=10*var, rec,...)
  b <- mutate(b, var=var)
  
  sp <- simple_split(b$haplotype, 1:var)
  set_leaf_position(sp, b$var+1)
  calc_node_ranges(sp,  gap=gap)
  
  blocks<- get_blocks(sp)
  plot(-1,-1, xlim=range(blocks[, 1:2]), 
       ylim=range(c(blocks[, 3:4], blocks[, 3:4]+blocks[, 5])), 
       type="n", axes=FALSE, xlab="", ylab="")
  apply(blocks, 1, plot_block)
  }
  
test_plot(1000, 0.0001, 50, gap=20)
test_plot(1000, 0.001, 50, gap=20)
test_plot(1000, 0.01, 50, gap=20)
test_plot(1000, 0.1, 50, gap=20)
par(opar)
 
## Does this allow use to look at the correlation between a pair of sites.  
  

# this can also be done with population subdivision
  




*/
