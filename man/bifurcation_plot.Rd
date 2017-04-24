\name{bifurcation_plot}
\alias{bifurcation_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a bifurcation diagram.
}
\description{
Takes the centre of the plot and a set of haplotypes and produces a bifurcation 
diagram to the left and the right.
}
\usage{
bifurcation_plot(haps, position, centre, nleft, nright, gap = 10, realpositions = TRUE, mark_variants = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{haps}{
A matrix of 0s and 1's which gives the haplotypes by row
}
  \item{position}{
The positions of the SNPS (if used)
}
  \item{centre}{
The SNP to be used as the centre.  If we take the SNPs to be labelled 1 to n, then the split at position
k means that SNPs to the left of position k, starting with k-1, k-2 are split for the left tree.
}
  \item{nleft}{
The number of SNPs to the left to split
}
  \item{nright}{
The number of SNPs to the right to split.  If we have n SNPs in total, then centre+nright must be less than n.
}
  \item{gap}{
The gap (given relative to the number of individuals) between branches to the right
}
  \item{realpositions}{
%%     ~~Describe \code{realpositions} here~~
}
  \item{mark_variants}{
%%     ~~Describe \code{mark_variants} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
