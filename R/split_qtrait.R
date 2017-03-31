

"qtrait_test" <- function(d, qtrait, positions, maxk = 4, reps = 1000, pickStat = "A") {
    if (missing(positions)) {
        positions <- 1:ncol(d)
    }
    if (length(positions) != ncol(d)) {
        stop("length of positions does not match columns")
    }
    n <- nrow(d)
    if (length(qtrait) != n) {
        stop("number of QTLs does not match number of haplotypes")
    }
    
    if (!(pickStat %in% c("A", "Z", "P", "N"))) {
        stop("Error, statPick should be one of A,P,Z,N")
    }
    
    rs <- rcppsplittestqtrait(d, qtrait, positions, reps, maxk, pickStat)
  
    p_ranks <- apply(rs, 2, function(x) (rank(x)[1]))
    p <- 1 - (p_ranks - 0.5)/(reps + 1)
    
    list(testStat = rs[1,], randTestStats = rs[-1,], p.value = p)
    
}


