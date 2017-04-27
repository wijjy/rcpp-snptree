library(rcppsnptree)
data(snptreeExample)
bifurcation_plot(haps, position, centre = 13, nleft=12, nright=10, gap=100)
bifurcation_plot(haps, position, centre = 1, nleft=0, nright=24, gap=100)
bifurcation_plot(haps, position, centre = 24, nleft=23, nright=1, gap=100)
bifurcation_plot(haps, position, centre = 13, nleft=7, nright=10, gap=100)




library(rcppsnptree)
data(snptreeExample)
split_right <- simple_split(haps, 13:24)
get_blocks(split_right, 10)



library(ARG)


