
"split_simple" <- function(d, SplitPositions, positions, quiet = TRUE) {
  
    if (missing(SplitPositions)) 
        SplitPositions <- 1:ncol(d)
    if (missing(positions)) 
        positions <- 1:ncol(d)
    if (length(positions) != ncol(d)) 
        stop("length of positions does not match columns")
    if (min(SplitPositions) < 1 || max(SplitPositions) > ncol(d)) 
        stop("values in Splitpositions should lie between 1 and ncol(d)")
    

    a <- rcpp_split_simple( d, SplitPositions-1);

    a$labels <- tapply(a$labels, rep(1:a$nleaves, a$leafcount), c)

    bb <- list(edge = a$edge, 
               Nnode = a$nleaves - 1, 
               edge.length = rep(1, 2 * (a$nleaves - 1)), 
               tip.label = a$labels)
    
    class(bb) <- c("phylo")
    
    haplotype_string <- apply(d, 1, paste, collapse = "")
    firstlabel <- sapply(a$labels, function(x) x[1])  ## haplotypes at the 
    tip_haplotypes <- haplotype_string[firstlabel]
    
    if (!quiet) 
        cat(leaves, " leaves on tree\n")
    
    b <- list(tree = bb, 
              nodepos = a$nodepos[1:(a$nleaves - 1)] + 1, 
              n = nrow(d), 
              labels = a$labels, 
              tip.haplotypes = tip_haplotypes)
    
    class(b) <- c("split")
    return(b)
}

