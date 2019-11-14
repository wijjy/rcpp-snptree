\name{get_blocks}
\alias{get_blocks, get_id_blocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the positions of blocks that build up a bifurcation diagram from a snptree.
}
\description{

}
\usage{
get_blocks(ptr, gap = 1)
get_id_blocks(ptr, id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptr}{
%%     ~~Describe \code{ptr} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
}
\item{id}{

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
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
set_leaf_position}
\examples{
data(snptreeExample)

split_right <- simple_split(haps, 1:20)
set_leaf_position(split_right, 25)
calc_node_ranges(split_right, 100)
leaf(split_right, 3)

split_left <- simple_split(haps,12:1)
set_leaf_position(split_left, 0)
calc_node_ranges(split_left, 100)
leaf(split_left, 3)
node(split_right, 6)
id <- sort(sample(nrow(haps), 100))

blocks_right <- get_blocks(split_right, gap=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
