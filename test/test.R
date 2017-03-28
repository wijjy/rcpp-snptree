data("snptreeExample")
a <- simple_split(haps, 5:ncol(haps))
nleaves(a)
cases <- 1:1000
case_control_leaves(a,cases)
leaf_count(a)



plot_stump <- function(x1, x2, y1, h1, y2, h2, gap) {
  x <- c(x1, (x1+x2)/2,    x2,    x2,    x1,        x2,        x2,    x1)
  y <- c(y1, y1+h2       ,    y2, y2+h2, y1+h2, y2+h2+gap, y2+gap+h1, y1+h1)
  s <- c(0,       0.5,     0,     0,   0.5,         0,         0,     0)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
  xspline(x,y, shape=s, open=FALSE, repEnds=TRUE, col="lightgrey", border=FALSE)
  
}

plot_stump(0.5, 1, 0.4, 0.2, 0.6, 0.15, 0.05)