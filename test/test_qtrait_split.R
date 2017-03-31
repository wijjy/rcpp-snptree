library(rcppsnptree)
library(ape)

## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
data(snptreeExample)
## simulate a trait
trait1 <- rnorm(nrow(haps))
s <- simple_split(haps, position)
sstat <- qtrait_statistics(s, trait1, statPick="A")
nl <- leaf_count(s)
plot(get_phylo(s), show.tip.label = FALSE)
tiplabels(round(sstat, 2), cex = 0.5)

hapstring <- factor(apply(haps,1,paste,collapse=""))
hapstring <- as.integer(hapstring)


boxplot(trait1 ~ hapstring)
anova(lm(trait1~hapstring))
tst <- qtrait_test(haps, trait1)
tst$p.value

s <- simple_split(haps, position)
samples <- node_labels(s)

## now add a small value equal to 1/5 of a standard deviation nto all those in
## the last haplotype
trait1[unlist(samples[1:9])] <- trait1[unlist(samples[1:9])] + 0.2
boxplot(trait1 ~ hapstring)

tst <- qtrait_test(haps, trait1, maxk = 9, pickStat = "A")
tst$p.value

sstat <- qtrait_statistics(s, trait1, statPick="A")
plot(get_phylo(s), show.tip.label = FALSE)
hl <- cut(sstat, c(-1e+100, -2, 2, 1e+100))
cols <- c("green", "white", "red")[as.numeric(hl)]
tiplabels(round(sstat, 2), cex = 0.5, bg=cols)

anova(lm(trait1~hapstring))
