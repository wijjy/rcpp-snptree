\name{leaf}
\alias{leaf, node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get information about nodes or leaves on a tree.
}
\description{
Returns information about nodes or leaves on a tree.  The index uses the same 
order as the APE tree ordering.
}
\usage{
leaf(ptr, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptr}{
A pointer to the tree}
  \item{index}{
The index to take (o offset).
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
Ian Wilson \emai{ian.wilson@newcastle.ac.uk}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ptr, index) 
{
    .Call("rcppsnptree_leaf", PACKAGE = "rcppsnptree", ptr, index)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
