
bifurcation_plot <- function(haps, position, centre, nleft, nright, 
                      gap=10, realpositions=TRUE, mark_variants=FALSE) {
  #ss=100;rec=0.01;var=100;centre=50;nleft=10;nright=10;gap=10

  len <- length(position)
  if (centre-nleft < 1)
    stop("The centre is between centre-1 and centre.  nleft is too large")
  if (centre+nright>=len)
    stop("Too many SNPs to the right of the centre")
  if (centre<len || centre > 1) {
    split_right <- simple_split(haps, (centre):(centre+nright-1))
    coordinates_right <- get_coordinates(split_right, gap=gap) 
    coordinates_right[coordinates_right[,1]==0, 1] <- centre+nright-1 
    split_left <- simple_split(haps, (centre-1):(centre-nleft))
    coordinates_left <- get_coordinates(split_left, gap=gap)
    coordinates_left[coordinates_left[,1]==0, 1] <- centre-nleft
  } else {
    stop("code currently only written for a SNP in the centre")
  }
  ## Line up the left and right plots
  height_diff <- coordinates_left[1, 2] - coordinates_right[1, 2]
  coordinates_right[,2] <- coordinates_right[,2]+height_diff
  ## Get a joining peice for the middle
  ry <- range(c(coordinates_left[,2], coordinates_right[,2]))
  x <- c(centre-1, centre-1, centre, centre)
  y <- c(coordinates_left[1, 2], coordinates_right[nrow(coordinates_right), 2])
  y <- c(y, rev(y))
  
  if (realpositions) {
    coordinates_right[,1] <- position[coordinates_right[,1]]
    coordinates_left[,1] <- position[coordinates_left[,1]]
    x <- position[x]
  }
  
  opar=par(mar=rep(1,4))
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