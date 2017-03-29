data("snptreeExample")
a <- simple_split(haps, position)
nleaves(a)
cases <- 1:1000
case_control_leaves(a,cases)
leaves <- leaf_count(a)
sum(leaves)
res <- get_coordinates(a, gap=40) 
res[res[,1]==-1,1] <- max(res[,1])+1 
plot(range(res[,1]), range(res[,2]), type="n", xlab="", ylab="", axes=TRUE)
xspline(res[,1], res[,2], shape=res[,3], open=FALSE, col="lightgrey", border=FALSE)




library(ARG)
b <- simARG(200, sites=1000, rec=0.01)
b <- mutate(b, var=100)

split_b <- simple_split(b$haplotype, 1:10)
nleaves(split_b)
cases <- 1:100
case_control_leaves(split_b, cases)
leaves <- leaf_count(split_b)
sum(leaves)


res <- get_coordinates(split_b, gap=10) 
res[res[,1]==-1,1] <- max(res[,1])+1 

plot(range(res[,1]), range(res[,2]), type="n", xlab="", ylab="", axes=TRUE)
xspline(res[,1], res[,2], shape=res[,3], open=FALSE, col="lightgrey", border=FALSE)



split_right <- simple_split(b$haplotype, 51:60)
coordinates_right <- get_coordinates(split_right, gap=10) 
coordinates_right[coordinates_right[,1]==-1,1] <- max(coordinates_right[,1])+1 

split_left <- simple_split(b$haplotype, 50:40)
coordinates_left <- get_coordinates(split_left, gap=10)
coordinates_left[coordinates_left[,1]==-1,1] <- min(coordinates_left[coordinates_left[,1]!=-1,1])-1 

height_diff <- coordinates_left[1, 2] - coordinates_right[1, 2]
coordinates_right[,2] <- coordinates_right[,2]+height_diff

x <- c(max(coordinates_left[,1])-0.02, min(coordinates_right[,1])+0.02)
x <- rep(x, each=2)
y <- c(coordinates_left[1,2], coordinates_right[nrow(coordinates_right),2])
y <- c(y, rev(y))

plot(range(c(coordinates_left[,1], coordinates_right[,1])), 
     range(c(coordinates_left[,2], coordinates_right[,2])), 
           type="n", xlab="", ylab="", axes=TRUE)
xspline(coordinates_left[,1], coordinates_left[,2], shape=coordinates_left[,3], open=FALSE, col="lightgrey", border=FALSE)

xspline(coordinates_right[,1], coordinates_right[,2], shape=coordinates_right[,3], open=FALSE, col="lightgrey", border=FALSE)
xspline(x, y, shape=0, open=FALSE, col="lightgrey", border="lightgrey")




plot_stump <- function(x1, x2, y1, h1, y2, h2, gap) {
  x <- c(x1, (x1+x2)/2,    x2,    x2,    x1,        x2,        x2,    x1)
  y <- c(y1, y1+h2       ,    y2, y2+h2, y1+h2, y2+h2+gap, y2+gap+h1, y1+h1)
  s <- c(0,       0.5,     0,     0,   0.5,         0,         0,     0)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
  xspline(x,y, shape=s, open=FALSE, repEnds=TRUE, col="lightgrey", border=FALSE)
}

plot_stump(0.5, 1, 0.4, 0.2, 0.6, 0.15, 0.05)
