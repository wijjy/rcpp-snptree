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

add_block <- function(blocks, num, ...) {
  plot_block(blocks[num,], ...)
}

recolour_block <- function(blocks, col="green") {
  w <- locate_block(blocks)
  add_block(blocks, w, col=col, border=col)
}

bifurc <- function(haplotypes, centre_position, nleft, nright, gap=100) {
  nsnps <- ncol(haplotypes)  
  if (centre_position-nleft < 1)
    stop("The centre is between centre-1 and centre.  There are not enough SNPs to the left.  nleft is too large")
  if (centre_position + nright -1 > len)
    stop("Too many SNPs to the right of the centre")
  if (centre_position < len || centre_position > 1) {
    split_right <- simple_split(haplotypes, (centre_position):(centre_position+nright-1))
    set_leaf_position(split_right, centre_position+nright)
    blocks_right <- get_blocks(split_right, gap=gap)
    split_left <- simple_split(haps, (centre_position-1):(centre_position-nleft))
    set_leaf_position(split_left, centre_position-nleft-1)
    blocks_left <- get_blocks(split_left, gap=gap)
  } else {
    stop("code currently only written for a SNP in the centre")
  }
  ## Line up the left and right plots
  height_diff <- blocks_left[1, 3] - blocks_right[1, 3]
  blocks_right[,3:4] <- blocks_right[,3:4]+height_diff
  ## Get a joining piece for the middle
  x <- c(centre-1, centre-1, centre, centre)
  y0 <- blocks_left[1,3]
  y1 <- y0 + blocks_left[1,5] + blocks_left[2,5]
  y <- c(y0, y1, y1, y0)

  if (realpositions) {
    blocks_left[, 1] <- position[blocks_left[, 2]]
    blocks_left[, 2] <- position[blocks_left[, 2]]
    blocks_right[, 1] <- position[blocks_right[, 2]]
    blocks_right[, 2] <- position[blocks_right[, 2]]   
    x <- position[x]
  }
  
  opar=par(mar=rep(1,4))
  plot_range <- range(c(blocks_left[,1], blocks_left[,2], blocks_right[,1] , blocks_right[,2])) 
  plot(range(c(coordinates_left[,1], coordinates_right[,1])), ry, 
       type="n", xlab="", ylab="", axes=FALSE)
  xspline(coordinates_left[,1], coordinates_left[,2], shape=coordinates_left[,3], open=FALSE, col="lightgrey", border=FALSE)
  xspline(coordinates_right[,1 ], coordinates_right[,2], shape=coordinates_right[,3], open=FALSE, col="lightgrey", border=FALSE)
  xspline(x, y, shape=0, open=FALSE, col="lightgrey", border="lightgrey")
  
  if (mark_variants) {
    pos <- unique(c(coordinates_left[,1], coordinates_right[,1]))
    #rug(pos)
    lapply(pos, function(x) segments(x, 0, x, ry[2], lwd=0.5, col="white"))
  }
  
  par(opar)
}
  
