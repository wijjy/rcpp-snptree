plot_block <- function(v,  col="lightgrey", ...) {
  x <- c(v[1], (3*v[1]+v[2])/4, v[2], v[2]     ,    (3*v[1]+v[2])/4, v[1])
  y <- c(v[3],   (v[3]+v[4])/2, v[4], v[4]+v[5], (v[3]+v[4])/2+v[5], v[3]+v[5])
  s <- c(   0,              -1,    0,         0,                 -1,    0)
  
  xspline(x,y, col=col, shape=s,  open=FALSE, ...)
}

locate_block <- function(blocks) {
  click <- locator(n=1)
  ux <- (click$x > blocks[,1] & click$x < blocks[,2]) | (click$x > blocks[,2] & click$x < blocks[,1])
  uy <- click$y > pmin(blocks[,3], blocks[,4])  & click$y <  pmax(blocks[,3], blocks[,4]+ blocks[,5])
  which(ux & uy)
}


mark_individuals <- function(b) {
  ## click once on each side of the central block and draw those individual 
  ## in both in a dfferent colour
  click_right <- locate_block(b$blocks_right)
  click_left <- locate_block(b$blocks_left)
  
  
  
}


mark_individual <- function(b) {
  ## click once on each side of the central block and draw those individual 
  ## in both in a dfferent colour
  
  click_right <- locate_block(b$blocks_right)
  click_left <- locate_block(b$blocks_left)
  
  
  
}


add_block <- function(blocks, num, ...) {
  plot_block(blocks[num,], ...)
}

recolour_block <- function(blocks, col="green") {
  w <- locate_block(blocks)
  add_block(blocks, w, col=col, border=col)
}



bifurc <- function(haplotypes, centre_position, nleft, nright, gap=100, 
                   realpositions=FALSE, mark_variants=FALSE, col="lightgrey") {
  nsnps <- ncol(haplotypes)  
  if (centre_position-nleft < 1)
    stop("The centre is between centre-1 and centre.  There are not enough SNPs to the left.  nleft is too large")
  if (centre_position + nright -1 > nsnps)
    stop("Too many SNPs to the right of the centre")
  if (centre_position < nsnps || centre_position > 1) {
    split_right <- simple_split(haplotypes, (centre_position):(centre_position+nright-1))
    set_leaf_position(split_right, centre_position+nright-1)
    blocks_right <- get_blocks(split_right, gap=gap)
    split_left <- simple_split(haps, (centre_position-1):(centre_position-nleft))
    set_leaf_position(split_left, centre_position-nleft-2)
    blocks_left <- get_blocks(split_left, gap=gap)
  } else {
    stop("code currently only written for a SNP in the centre")
  }
  ## Line up the left and right plots
  height_diff <- blocks_left[1, 3] - blocks_right[1, 3]
  blocks_right[,3:4] <- blocks_right[,3:4]+height_diff
  ## Get a joining piece for the middle
  x <- c(centre_position-1, centre_position-1, centre_position, centre_position)
  y0 <- blocks_left[1, 3]
  y1 <- y0 + blocks_left[1,5] + blocks_left[2, 5]
  y <- c(y0, y1, y1, y0)

  if (realpositions) {
    blocks_left[, 1] <- position[blocks_left[, 2]]
    blocks_left[, 2] <- position[blocks_left[, 2]]
    blocks_right[, 1] <- position[blocks_right[, 2]]
    blocks_right[, 2] <- position[blocks_right[, 2]]   
    x <- position[x]
  }
  
  opar=par(mar=rep(1,4))
  range_x <- range(c(blocks_left[,1], blocks_left[,2], blocks_right[,1] , blocks_right[,2])) 
  range_y <- range(c(
    blocks_left[, 3:4], 
    blocks_left[, 3:4]+blocks_left[,5], 
    blocks_right[, 3:4], blocks_right[, 3:4]+blocks_right[, 5]))
  plot(range_x, range_y, type="n", xlab="", ylab="", axes=FALSE)
  apply(blocks_left, 1, plot_block, col=col)
  apply(blocks_right, 1, plot_block, col=col)
  xspline(x, y, shape=0, open=FALSE, col=col)
  
  if (mark_variants) {
    pos <- unique(c(range_left[,1], range_right[,1]))
    #rug(pos)
    lapply(pos, function(x) segments(x, 0, x, ry[2], lwd=0.5, col="white"))
  }
  
  par(opar)
  list(blocks_right=blocks_right, blocks_left=blocks_left, 
       left=split_left, right=split_right, centre=centre_position, 
       nleft=nleft, nright=nright)
}


# This bifurcation plot takes the genotype of the centre and uses it to colour the 
# rest of the plot.
# So we can have centre-1 to the left and nsnps-centre to the right in total
# this is another way of exploring the haplotypes.

bifurcb <- function(haplotypes, centre_position, nleft, nright, gap=100, 
                   realpositions=FALSE, mark_variants=FALSE, col="lightgrey") {
  nsnps <- ncol(haplotypes)  
  if (centre_position-nleft < 1)
    stop("The centre is between centre-1 and centre.  There are not enough SNPs to the left.  nleft is too large")
  if (centre_position + nright  > nsnps)
    stop("Too many SNPs to the right of the centre")
  if (centre_position < nsnps || centre_position > 1) {
    split_right <- simple_split(haplotypes, (centre_position+1):(centre_position+nright))
    set_leaf_position(split_right, centre_position+nright+1)
    blocks_right <- get_blocks(split_right, gap=gap)
    split_left <- simple_split(haps, (centre_position-1):(centre_position-nleft))
    set_leaf_position(split_left, centre_position-nleft-2)
    blocks_left <- get_blocks(split_left, gap=gap)
  } else {
    stop("code currently only written for a SNP in the centre")
  }
  ## Line up the left and right plots
  height_diff <- blocks_left[1, 3] - blocks_right[1, 3]
  blocks_right[,3:4] <- blocks_right[,3:4]+height_diff
  ## Get a joining piece for the middle
  x <- c(centre_position-1, centre_position-1, centre_position+1, centre_position+1)
  y0 <- blocks_left[1, 3]
  y1 <- y0 + blocks_left[1,5] + blocks_left[2, 5]
  y <- c(y0, y1, y1, y0)
  
  if (realpositions) {
    blocks_left[, 1] <- position[blocks_left[, 2]]
    blocks_left[, 2] <- position[blocks_left[, 2]]
    blocks_right[, 1] <- position[blocks_right[, 2]]
    blocks_right[, 2] <- position[blocks_right[, 2]]   
    x <- position[x]
  }
  
  opar=par(mar=rep(1,4))
  range_x <- range(c(blocks_left[,1], blocks_left[,2], blocks_right[,1] , blocks_right[,2])) 
  range_y <- range(c(
    blocks_left[, 3:4], 
    blocks_left[, 3:4]+blocks_left[,5], 
    blocks_right[, 3:4], blocks_right[, 3:4]+blocks_right[, 5]))
  plot(range_x, range_y, type="n", xlab="", ylab="", axes=FALSE)
  apply(blocks_left, 1, plot_block, col=col)
  apply(blocks_right, 1, plot_block, col=col)
  xspline(x, y, shape=0, open=FALSE, col=col)
  
  if (mark_variants) {
    pos <- unique(c(range_left[,1], range_right[,1]))
    #rug(pos)
    lapply(pos, function(x) segments(x, 0, x, ry[2], lwd=0.5, col="white"))
  }
  
  par(opar)
  list(blocks_right=blocks_right, blocks_left=blocks_left, 
       left=split_left, right=split_right, centre=centre_position, 
       nleft=nleft, nright=nright)
}

get_bifurc <- function(haplotypes, positions_left, positions_right,  gap=100) {
  nsnps <- ncol(haplotypes)
  r <- range(c(positions_left, positions_right))
  
  if (r[1] <1 || r[2] > nsnps)
    stop("positions should be between 0 and nsnps")
  
  if(max(split_left) >= min(split_right))
    stop("Left splits should all be less than right splits")

  split_right <- simple_split(haplotypes, positions_right)
  split_left <- simple_split(haplotypes, positions_left)
    
  set_leaf_position(split_right, r[2] + 1)
  set_leaf_position(split_left, r[1] - 1)
  
  blocks_right <- get_blocks(split_right, gap=gap)
  blocks_left <- get_blocks(split_left, gap=gap)
 
  ## Line up the left and right plots
  height_diff <- blocks_left[1, 3] - blocks_right[1, 3]
  blocks_right[,3:4] <- blocks_right[,3:4]+height_diff
  centre_left <- max(positions_left)
  centre_right <- min(positions_right)
  ## Get a joining piece for the middle
  y0 <- blocks_left[1, 3]
  y1 <- y0 + blocks_left[1,5] + blocks_left[2, 5]
  centre <- data.frame(x=c(centre_left, centre_left, centre_right, centre_right), 
                       y=c(y0, y1, y1, y0))
  
  if (realpositions) {
    blocks_left[, 1] <- position[blocks_left[, 2]]
    blocks_left[, 2] <- position[blocks_left[, 2]]
    blocks_right[, 1] <- position[blocks_right[, 2]]
    blocks_right[, 2] <- position[blocks_right[, 2]]   
    x <- position[x]
  }
  
  range_x <- range(c(blocks_left[,1], blocks_left[,2], blocks_right[,1] , blocks_right[,2])) 
  range_y <- range(c(
    blocks_left[, 3:4], 
    blocks_left[, 3:4]+blocks_left[,5], 
    blocks_right[, 3:4], blocks_right[, 3:4]+blocks_right[, 5]))
 
  res <- list(blocks_right=blocks_right, 
              blocks_left=blocks_left, 
              left=split_left, 
              right=split_right, 
              centre=centre,
              range_x=range_x, 
              range_y <- range_y,
              positions_left=positions_left,
              positions_right=positions_right,
              haplotypes=haplotypes)
  
  class(res) <- "bifurc"
  return(res)
  
}


plot_bifurc <- function(bb) {
  
}

if (FALSE) {
  library(rcppsnptree)
  data(snptreeExample)
  bifurcb(haps, 13, 12, 11, col="lightblue1")
}


  
