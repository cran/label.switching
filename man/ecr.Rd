\name{ecr}
\alias{ecr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ECR algorithm (default version)
}
\description{
This function applies the standard version of Equivalence Classes Representatives (ECR) algorithm (Papastamoulis and Iliopoulos, 2010). The set of all allocation variables is partitioned into equivalence classes and exactly one representative is chosen from each class. The practical implementation of this idea is to reorder the output so that all simulated allocation vectors (\code{z}) are as similar as possible with a pivot allocation vector (\code{zpivot}). The useR has to specify this pivot allocation vector as a good allocation of the observations among the mixture components. Some typical choices are the allocations that correspond to the complete or non-complete MAP estimate.  
}
\usage{
ecr(zpivot, z, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zpivot}{\eqn{n}-dimensional integer vector \eqn{(z_1,\ldots,z_n)} with \eqn{z_i\in\{1,\ldots,K\}}, \eqn{i=1,\ldots,n}. 
}
  \item{z}{
%%     ~~Describe \code{z} here~~
\eqn{m\times n} integer array of the latent allocation vectors generated from an MCMC algorithm.
}
  \item{K}{
%%     ~~Describe \code{K} here~~
the number of mixture components (at least equal to 2).
}
}
\details{
\code{zpivot} should be chosen as an allocation vector that corresponds to a high-posterior density area, or in general as an allocation that is considered as a good allocation of the observations among the \eqn{K} components.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{permutations }{\eqn{m\times K} dimensional array of permutations}

}
\references{
Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.
}
\author{
Panagiotis Papastamoulis
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{permute.mcmc}}, \code{\link{label.switching}}, \code{\link{ecr.iterative.1}}, \code{\link{ecr.iterative.2}} 
}
\examples{
#load a toy example: MCMC output consists of the random beta model
#	applied to a normal mixture of \code{K=2} components. The
# 	number of observations is equal to \code{n=5}. The number
#	of MCMC samples is equal to \code{m=300}. The 300 
#	simulated allocations are stored to array \code{z}. The 
#	complete MAP estimate corresponds to iteration \code{mapindex}.
data("mcmc_output")
z<-data_list$"z"
K<-data_list$"K"
mapindex<-data_list$"mapindex"

# mcmc parameters are stored to array \code{mcmc.pars}
mcmc.pars<-data_list$"mcmc.pars"
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances 
# mcmc.pars[,,3]: simulated weights
run<-ecr(zpivot = z[mapindex,],z = z, K = K)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances
# reordered.mcmc[,,3]: reordered weights
}
