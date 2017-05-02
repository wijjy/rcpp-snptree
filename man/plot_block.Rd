\name{plot_block}
\alias{plot_block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a single block.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_block(v, col = "lightgrey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (v, col = "lightgrey", ...) 
{
    x <- c(v[1], (3 * v[1] + v[2])/4, v[2], v[2], (3 * v[1] + 
        v[2])/4, v[1])
    y <- c(v[3], (v[3] + v[4])/2, v[4], v[4] + v[5], (v[3] + 
        v[4])/2 + v[5], v[3] + v[5])
    s <- c(0, -1, 0, 0, -1, 0)
    xspline(x, y, col = col, shape = s, open = FALSE, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
