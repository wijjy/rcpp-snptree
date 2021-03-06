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
  if (s->root()->range.first==0)
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
/** Set the location of the leaf nodes.  Since these are generally past the last 
 * split there is no information on the location of the split passed to the function
 */
// [[Rcpp::export]]
void set_leaf_position(SEXP ptr, double pos) {
  Rcpp::XPtr< splitter > s(ptr);
  std::list<binode *>::iterator bi=s->begin_leaf();
  while (bi!=s->end_leaf()) {
    (*bi)->position=pos;
    bi++;
  }
}

/** Return a representation of a node.  I use the same 
 * indexing as ape, so that we have 1..nleaves for the leaves and
 * then node nleaves+1 is the deepest node on the left brnach than 
 * is not a leaf.  So inorder" on the nodes
 */

// [[Rcpp::export]]
Rcpp::List leaf(SEXP ptr, int index) {
  Rcpp::XPtr< splitter > s(ptr);
  LRNIterator<binode> ii(s->root());

  ii.nextLeaf();   // the root can never be a leaf
  while (!ii.isend()) {
    if (index==1) {
       Rcpp::NumericVector range(2);
      range(0) = (*ii)->range.first;
      range(1) = (*ii)->range.second;
      
      return Rcpp::List::create(Rcpp::Named("position")=(*ii)->position,
                         Rcpp::Named("range")=range);
    }
    
    ii.nextLeaf();
    index--;
  }
}

// [[Rcpp::export]]
Rcpp::List internal_node(SEXP ptr, int index) {
  Rcpp::XPtr< splitter > s(ptr);
  LRNIterator<binode> ii(s->root());
  while (!ii.isend()) {
    if (index==1) {
      return (*ii)->list_node();
    }
    ii.nextInternal();
    index--;
  }
  return 0;
}
/** Return a representation of a node that should be in the same order as blocks  */
// [[Rcpp::export]]
Rcpp::List nodeb(SEXP ptr, int index) {
  Rcpp::XPtr< splitter > s(ptr);
  NLRIterator<binode> ii(s->root());
  while (!ii.isend()) {
    if (index==1) {
      return (*ii)->left->list_node_long();
    } else if (index==2) {
      return (*ii)->right->list_node_long();
    }
    ii.nextInternal();
    index -= 2;
  }
  return 0;
}

  
  
// [[Rcpp::export]]
void calc_node_ranges(SEXP ptr, double gap, bool log=false) {
  Rcpp::XPtr< splitter > s(ptr);    
  if (log==FALSE)
    s->calculate_top_bottom(gap); 
  else
    s->calculate_top_bottom_log2(gap); 
  
}




// [[Rcpp::export]]
Rcpp::NumericMatrix get_blocks2(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  if (s->root()->range.first==0)
    s->calculate_top_bottom(gap);   // gets the tops and bottoms 
  
  int leaves = s->nleaves();
  Rcpp::NumericMatrix boxes(2*leaves-2, 5);  // last one left for the root if needed
  
  int index=0;
  std::list<binode *>::const_iterator ii=s->begin_internal(); 
                                     
  while (ii != s->end_internal()) {
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
    index++;
    ii++;
  }
  return boxes;
}
 
 
// [[Rcpp::export]]
Rcpp::IntegerVector individuals_matching_haplotype(SEXP ptr, const Rcpp::IntegerVector &haplotype, 
                                                   const Rcpp::IntegerVector &position) {
  // haplotype can have 0, 1 and -1 (to act as anything)
  // Get the original tree whihc is to act as a scaffold
  Rcpp::XPtr< splitter > s(ptr);
  
  // now get another tree just with the required nodes.  special split just one way
  splitter tmp(*s);                           // define the temporary spliiter object
  for (int i=0;i< position.size();i++) {
    if (haplotype[position[i]] != -1) 
      tmp.split(position[i]-1);   // split at positions
  }
  binode *bl = tmp.first_leaf_matching(position, haplotype);
  return Rcpp::IntegerVector(bl->labels.begin(), bl->labels.end())-1;
  
}

 
/*** R
plot_block <- function(v,  col="lightgrey", ...) {
  x <- c(v[1], v[2], v[2], v[1])
  y <- c(v[3], v[4], v[4]+v[5], v[3]+v[5])
  xspline(x,y, col=col, border=col, shape=0, ...)
}
plot_block <- function(v,  col="lightgrey", ...) {
    x <- c(v[1], (3*v[1]+v[2])/4, v[2], v[2]     ,    (3*v[1]+v[2])/4, v[1])
    y <- c(v[3],   (v[3]+v[4])/2, v[4], v[4]+v[5], (v[3]+v[4])/2+v[5], v[3]+v[5])
    s <- c(   0,              -1,    0,         0,                 -1,    0)

  xspline(x,y, col=col, shape=s, border="darkgrey", open=FALSE, ...)
}
library(rcppsnptree)
data(snptreeExample)

split_right <- simple_split(haps, 13:24)
set_leaf_position(split_right, 24)
calc_node_ranges(split_right, 100)
leaf(split_right, 3)
blocks_left <- get_blocks(split_right, gap=100)

split_left <- simple_split(haps,12:1)
set_leaf_position(split_left, -1)
calc_node_ranges(split_left, 100)
blocks_left <- get_blocks(split_left, gap=100)

node(split_left, 3)




blocks_right <- get_blocks(split_right, gap=100)
blocks_right2 <- get_blocks2(split_right, gap=100)
blocks_left <- get_blocks(split_left, gap=100)
blocks_right <- blocks_right[-nrow(blocks_right), ]
blocks_left <-  blocks_left[-nrow(blocks_left), ]

plot(range(blocks_left[,1:2]), range(c(blocks_left[,3:4]),
                                     c(blocks_left[,3:4])+blocks_left[,5]), 
     axes=FALSE, xlab="", ylab="", type="n")

if (FALSE) {

plot_block(blocks_left[1,], col=2)
plot_block(blocks_left[2,], col=3)
plot_block(blocks_left[3,], col=4)
plot_block(blocks_left[4,], col=4)
plot_block(blocks_left[5,], col=4)
plot_block(blocks_left[6,], col=6)


apply(blocks_left, 1, plot_block)
axis(1)

plot(range(blocks_right[,1:2]), range(c(blocks_right[,3:4]),
                                     c(blocks_right[,3:4])+blocks_right[,5]), 
     axes=FALSE, xlab="", ylab="", type="n")

apply(blocks_right2, 1, plot_block)
axis(1)
plot_block(blocks[26,])
polygon(boxes[,1], boxes[,2])
*/


 


