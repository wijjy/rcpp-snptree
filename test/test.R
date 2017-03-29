data("snptreeExample")
a <- simple_split(haps, position)
nleaves(a)
cases <- 1:1000
case_control_leaves(a,cases)
leaves <- leaf_count(a)
sum(leaves)
res <- get_coordinates(a, gap=100) 
res[res[,1]==-1,1] <- max(res[,1])+1 
plot(range(res[,1]), range(res[,2]), type="n", xlab="", ylab="")
xspline(res[,1], res[,2], shape=res[,3], open=FALSE, col="lightgrey", border=FALSE)


















plot_stump <- function(x1, x2, y1, h1, y2, h2, gap) {
  x <- c(x1, (x1+x2)/2,    x2,    x2,    x1,        x2,        x2,    x1)
  y <- c(y1, y1+h2       ,    y2, y2+h2, y1+h2, y2+h2+gap, y2+gap+h1, y1+h1)
  s <- c(0,       0.5,     0,     0,   0.5,         0,         0,     0)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
  xspline(x,y, shape=s, open=FALSE, repEnds=TRUE, col="lightgrey", border=FALSE)
}

plot_stump(0.5, 1, 0.4, 0.2, 0.6, 0.15, 0.05)