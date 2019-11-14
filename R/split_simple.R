
"split_simple" <- function(d, SplitPositions, positions, quiet = TRUE) {
  
    if (missing(SplitPositions)) 
        SplitPositions <- 1:ncol(d)
    if (missing(positions)) 
        positions <- 1:ncol(d)
    if (length(positions) != ncol(d)) 
        stop("length of positions does not match columns")
    if (min(SplitPositions) < 1 || max(SplitPositions) > ncol(d)) 
        stop("values in Splitpositions should lie between 1 and ncol(d)")
    

    a <- rcpp_split_simple(d, SplitPositions);

    haplotype_string <- apply(d, 1, paste, collapse = "")
    firstlabel <- sapply(a$labels, function(x) x[1])  ## haplotypes at the 
    tip_haplotypes <- haplotype_string[firstlabel]
    names(tip_haplotypes) <- 1:(a$Nnode+1)
    
    if (!quiet) cat(a$Nnode+1, " leaves on tree\n")
    
    b <- list(tree = a, 
              nodepos = a$nodepos[1:(a$Nnode)] + 1, 
              n = nrow(d), 
              tip.haplotypes = tip_haplotypes)
    
    class(b) <- c("split")
    return(b)
}

