library(rcppsnptree)

## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
data(snptreeExample)

## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control


trait1 <- rnorm(nrow(haps))
tst <- spli tQTLTest(haps, trait1)
tst$p.value

## now add a small value equal to 1/5 of a standard deviation nto all those in
## the last haplotype
trait1[s$labels[['45']]] <- trait1[s$labels[['45']]] + 0.2
tst <- splitQTLTest(haps, trait1)
tst$p.value