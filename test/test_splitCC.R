library(rcppsnptree)
data(snptreeExample)
## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
haplotype <- apply(haps, 1, paste, collapse="")
table(haplotype, sample)
chisq.test(table(haplotype, sample), simulate=TRUE, B=2000)
## We should be able to repeat this by looking at the tips of
## out haplotype tree only, whihc we can do by setting maxk=45
s <- splitTestCC(haps, which(sample=="Case"), reps=1000, pickStat="G", maxk=45)
print(s$p.value)

s <- splitTestCC(haps, which(sample=="Case"), reps=1000, pickStat="AbsSevon")
print(s$p.value)

s <- splitTestCC(haps, which(sample=="Case"), reps=10000, pickStat="Gtest")
print(s$p.value)

## try with a random sample of cases
random_cases <- sort(sample(nrow(haps), size=nrow(haps)/2, replace=FALSE))
s <- splitTestCC(haps , random_cases, reps=1000, pickStat="G")
print(s$p.value)

library(snptree)
snptree::random_cases <- sort(sample(nrow(haps), size=nrow(haps)/2, replace=FALSE))
