#include <Rcpp.h>

#include "rcpp_binode.h"
#include "splitter.h"
#include "ijw_rand.h" 

/** Get the coordinates of the outside of a bifurcation diagram as 
 * a series of points.  The points are x, y and the "curve" varible from xspline
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix  get_coordinates(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  s->calculate_top_bottom(gap);   // gets the tops and bottoms 
  // for the leaves and the rest of the tree
  int leaves = s->nleaves();
  Rcpp::NumericMatrix coords(2*leaves+3*(leaves-1), 3);
  
  s->get_coordinates(coords);
  for (int i=0; i<coords.nrow(); i++) {
    coords(i, 0) +=1; 
  }
  
  return coords;
}
/** Return a bifurcation diagram as a set of blocks, each block being a node 
 * in the bifurcation tree
 */

// [[Rcpp::export]]
Rcpp::NumericMatrix  get_blocks(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  s->calculate_top_bottom(gap);   // gets the tops and bottoms 
                                  // for the leaves and the rest of the tree
  int leaves = s->nleaves();
  Rcpp::NumericMatrix boxes(2*leaves-1, 6);  // last one left for the root if needed
  
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
    boxes(index, 5) = 0;     // This is an indicator for up, down etc.
    // right
    index++;
    boxes(index, 0) = (*ii)->position+1;
    boxes(index, 1) = (*ii)->right->position+1;
    boxes(index, 2) = (*ii)->range.first+(*ii)->left->height();
    boxes(index, 3) = (*ii)->right->range.first;
    boxes(index, 4) = (*ii)->right->height();
    boxes(index, 5) = 1.0;
    ii.nextInternal();
    index++;
  }
  
  return boxes;
}

/*** R
library(rcppsnptree)
data(snptreeExample)
split_right <- simple_split(haps, 1:24)
a <- get_coordinates(split_right, gap=300)
plot_block <- function(v,  col="lightgrey", ...) {
  x <- c(v[1], v[2], v[2], v[1])
  y <- c(v[3], v[4], v[4]+v[5], v[3]+v[5])
  polygon(x,y, col=col, border=col, ...)
}
blocks <- get_blocks(split_right, gap=100)
blocks <- blocks[-nrow(blocks), ]
blocks[blocks[,2]==0,2 ] <- 24
plot(range(blocks[,1:2]), range(c(blocks[,3:4]),c(blocks[,3:4])+blocks[,5]), axes=FALSE, xlab="", ylab="")

plot_block(blocks[1,], col=2)
plot_block(blocks[2,], col=3)
plot_block(blocks[3,], col=4)
plot_block(blocks[4,], col=4)
plot_block(blocks[5,], col=4)
plot_block(blocks[6,], col=6)


apply(blocks, 1, plot_block)
plot_block(blocks[26,])
polygon(boxes[,1], boxes[,2])
*/


 


