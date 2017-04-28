plot_block <- function(v,  col="lightgrey", ...) {
  x <- c(v[1], (3*v[1]+v[2])/4, v[2], v[2]     ,    (3*v[1]+v[2])/4, v[1])
  y <- c(v[3],   (v[3]+v[4])/2, v[4], v[4]+v[5], (v[3]+v[4])/2+v[5], v[3]+v[5])
  s <- c(   0,              -1,    0,         0,                 -1,    0)
  
  xspline(x,y, col=col, shape=s,  open=FALSE, ...)
}

locate_block <- function(blocks) {
  click <- locator(n=1)
  ux <- click$x > blocks[,1] & click$x < blocks[,2]
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

plot_left_right <- function(left_blocks, right_blocks) {
  
}
