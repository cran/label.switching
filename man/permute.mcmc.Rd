\name{permute.mcmc}
\alias{permute.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Reorder MCMC samples
}
\description{
This function applies the permutation returned by any relabelling algorithm to a simulated MCMC output.
}
\usage{
permute.mcmc(mcmc, permutations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc}{
%%     ~~Describe \code{mcmc} here~~
\eqn{m\times K\times J} array of simulated MCMC parameters.
}
  \item{permutations}{
%%     ~~Describe \code{permutations} here~~
\eqn{m\times K} dimensional array of permutations.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{output }{\eqn{m\times K\times J} array of reordered MCMC parameters.}
}
\author{
%%  ~~who you are~~
Panagiotis Papastamoulis
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{label.switching}}, \code{\link{ecr}}, \code{\link{ecr.iterative.1}}, \code{\link{ecr.iterative.2}},\code{\link{stephens}},\code{\link{pra}}, \code{\link{sjw}} 
}
\examples{
#load MCMC simulated data
data("mcmc_output")
mcmc.pars<-data_list$"mcmc.pars"
z<-data_list$"z"
K<-data_list$"K"

#apply \code{ecr.iterative.1} algorithm
run<-ecr.iterative.1(z = z, K = 2)
#reorder the MCMC output according to this method:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances of the components
# reordered.mcmc[,,3]: reordered weights of the two components
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
