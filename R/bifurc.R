
bifurcation <- function(haplotypes, positions_left, positions_right,  gap=100, log=FALSE) {
  nsnps <- ncol(haplotypes)
  r <- range(c(positions_left, positions_right))
  
  if (r[1] <1 || r[2] > nsnps)
    stop("positions should be between 0 and nsnps")
  
  if(max(positions_left) >= min(positions_right))
    stop("Left splits should all be less than right splits")
  
  split_right <- simple_split(haplotypes, positions_right)
  split_left <- simple_split(haplotypes, positions_left)
  
  set_leaf_position(split_right, r[2])
  set_leaf_position(split_left, r[1] - 2)
  
  calc_node_ranges(split_right, gap, log)
  calc_node_ranges(split_left, gap, log)
  
  blocks_right <- get_blocks(split_right)
  blocks_left  <- get_blocks(split_left)
  
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
  
  
  
  range_x <- range(c(blocks_left[,1], blocks_left[,2], blocks_right[,1] , blocks_right[,2])) 
  range_y <- range(c(
    blocks_left[, 3:4], 
    blocks_left[, 3:4]+blocks_left[,5], 
    blocks_right[, 3:4], blocks_right[, 3:4]+blocks_right[, 5])
  )
  
  res <- list(blocks_right=blocks_right, 
              blocks_left=blocks_left, 
              left=split_left, 
              right=split_right, 
              centre=centre,
              range_x=range_x, 
              range_y = range_y,
              positions_left=positions_left,
              positions_right=positions_right,
              haplotypes=haplotypes)
  
  class(res) <- "bifurc"
  return(res)
}


summary.bifurc <- function(bf) {
  cat("A bifurcating tree on", nrow(bf$haplotypes), "individuals\n")
  cat("A total of", ncol(bf$haplotypes),"SNPs are available to split,\n")
  cat(length(b$positions_left), "are split to the left and", length(b$positions_right), "to the right.\n" )
  cat("There are", nleaves(b$left),"leaves to the left and", nleaves(b$right), "to the right\n")
}

plot.bifurc <- function(bb, col="lightgrey", border="black", ...) {
  opar=par(mar=rep(1,4))
  plot(bb$range_x, bb$range_y, type="n", xlab="", ylab="", axes=FALSE)
  apply(bb$blocks_left, 1, plot_block, col=col, border=border)
  apply(bb$blocks_right, 1, plot_block, col=col, border=border)
  xspline(bb$centre$x, bb$centre$y, shape=0, open=FALSE, col=col, border=border)
  par(opar)
}


if (FALSE) {
  library(rcppsnptree)
  data(snptreeExample)
  b <- bifurcation(haps, 12:1, 13:24)
  print(class(b))
  plot(b)
  summary(b)
  b <- bifurcation(haps, 12:1, 13:24, log=TRUE)
  plot(b)
}

