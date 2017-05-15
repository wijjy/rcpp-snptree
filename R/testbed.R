## explore the splits

split_plot <- function(haps) {
  npos <- ncol(haps)
  res <- numeric(npos)
  positions <- 1:npos
  for (i in c(1:npos)) {
    sp <- simple_split(haps, positions)
    positions <- c(positions[npos], positions[-npos])
    res[i] <- nleaves(sp)
    print(positions)
  }
  res
  
}



if (FALSE) {
  library(rcppsnptree)
  data(snptreeExample)
  split_plot(haps)
}

