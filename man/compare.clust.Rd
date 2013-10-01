\name{compare.clust}
\alias{compare.clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Make all estimated clusters agree with a pivot allocation
}
\description{
Given a pivot allocation vector, a set of simulated allocations and a set of permutations from different relabelling algorithms, this function relabels the permutations so that all methods maximize their similarity with the pivot. This is helpful when comparing different different label switching algorithms. 
}
\usage{
compare.clust(pivot.clust,perms,z,K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pivot.clust}{
%%     ~~Describe \code{mcmc} here~~
a pivot allocation vector of the \eqn{n} observations among the \eqn{K} clusters.
}
  \item{perms}{
%%     ~~Describe \code{permutations} here~~
a list containing \eqn{f} permutation arrays, as returned by \code{\link{label.switching}} function.
}
  \item{z}{
%%     ~~Describe \code{permutations} here~~
a set of simulated allocation arrays.
}
  \item{K}{
%%     ~~Describe \code{permutations} here~~
number of mixture components
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{similarity }{\eqn{(f+1) K\times (f+1)} matrix containing the similarity coefficient of the resulting clusters.}
\item{clusters }{\eqn{f\times n} array of single best clusterings, relabelled in order to maximize their similarity with \code{pivot.clust}.}
}
\author{
%%  ~~who you are~~
Panagiotis Papastamoulis
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{label.switching}} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
