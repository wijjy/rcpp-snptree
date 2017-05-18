
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
              haplotypes=haplotypes,
              gap=gap)
  
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

recentre <- function(bb) {
  click <- locator(n=1)
  
  x <- round(click$x, 0)
  left_pos <- x:1
  right_pos <- (x+1):ncol(bb$haplotypes)
  dd <- bifurcation(bb$haplotypes, left_pos, right_pos, gap=bb$gap )
  plot(dd)
  return(dd)
}

add_id <- function(bb, labels, col="green") {
  blocks_id_left <- get_id_blocks(bb$left, labels)
  blocks_id_right <- get_id_blocks(bb$right, labels)
  apply(blocks_id_left, 1, plot_block, col=col)
  apply(blocks_id_right, 1, plot_block, col=col)
  
}



locate_bifurc_block <- function(bb) {
  centre_block <- matrix(c(bb$centre$x[c(1,3)], bb$centre$y[1], bb$centre$y[1], 
                           bb$centre$y[2]-bb$centre$y[1]), nrow=1)
  blocks <- rbind(bb$blocks_left, bb$blocks_right, centre_block)
  
  click <- locator(n=1)

  ux <- (click$x > blocks[,1] & click$x < blocks[,2]) | (click$x > blocks[,2] & click$x < blocks[,1])
  uy <- click$y > pmin(blocks[,3], blocks[,4])  & click$y <  pmax(blocks[,3], blocks[,4]+ blocks[,5])
  w <- which(ux & uy)
  if (length(w)!=1)
    stop("need a single block, try to click within the block")
  if (w==nrow(blocks))
    stop("you have selected the centre.  no change made to plot")
  
  print(w)
  if (w <= nrow(bb$blocks_left))
    print(nodeb(bb$left, w))
  else 
    print(nodeb(bb$right, w-nrow(bb$blocks_left)))
}




if (FALSE) {
  library(rcppsnptree)
  
#  library(ARG)
#  a <- simARG(c(250,250), 5000, r=0.001, growthmodel = "exponential(10)", migmatrix = "Island(2, 10)")
#  b <- mutate(a, var=50)
#  simhaps <- b$haplotype
#  simlocation <- b$location
#  simdata <- list(haps=simhaps, location=b$location) 
#  save(simdata, file="simdata.rda")
  data(simdata)

  b <- bifurcation(simdata$haps, 25:1, 26:50, gap=10,log=FALSE)
  
  print(class(b))
  plot(b)
  node <- locate_bifurc_block(b)
  add_id(b, node$labels+1)
  node2 <- locate_bifurc_block(b)
  add_id(b, node2$labels, col="red")
  summary(b)
  d <- recentre(b)
add_id(d, locate_bifurc_block(d)$labels, col="red")
}

