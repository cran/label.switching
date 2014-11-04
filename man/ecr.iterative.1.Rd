\name{ecr.iterative.1}
\alias{ecr.iterative.1}
\title{
ECR algorithm (iterative version 1)
}
\description{
This function applies the first iterative version of Equivalence Classes Representatives (ECR) algorithm (Papastamoulis and Iliopoulos, 2010, Rodriguez and Walker, 2012). The set of all allocation variables is partitioned into equivalence classes and exactly one representative is chosen from each class. The difference with the default version of ECR algorithm is that no pivot is required and the method is iterative, until a fixed pivot has been found.
}
\usage{
ecr.iterative.1(z, K, opt_init, threshold, maxiter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{\eqn{m\times n} integer array of the latent allocation vectors generated from an MCMC algorithm.
}
  \item{K}{the number of mixture components (at least equal to 2).
}
  \item{opt_init}{
An (optional) \eqn{m\times K} array of permutations to initialize the algorithm. The identity permutation is used if it is not specified.
}
  \item{threshold}{
An (optional) positive number controlling the convergence criterion. Default value: 1e-6.
}
  \item{maxiter}{
An (optional) integer controlling the max number of iterations. Default value: 100.
}
}
\value{
\item{permutations }{\eqn{m\times K} dimensional array of permutations}
\item{iterations }{integer denoting the number of iterations until convergence}
\item{status }{returns the exit status}
}
\references{
Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.

Rodriguez C.E. and Walker S. (2014). Label Switching in Bayesian Mixture Models: Deterministic relabeling strategies. Journal of Computational and Graphical Statistics. 23:1, 25-45
}
\author{
Panagiotis Papastamoulis
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{permute.mcmc}}, \code{\link{label.switching}}, \code{\link{ecr}}, \code{\link{ecr.iterative.2}}
}
\examples{
#load a toy example: MCMC output consists of the random beta model
# applied to a normal mixture of \code{K=2} components. The number of
# observations is equal to \code{n=5}. The number of MCMC samples is
# equal to \code{m=1000}. The 300 simulated allocations are stored to
# array \code{z}. 
data("mcmc_output")
# mcmc parameters are stored to array \code{mcmc.pars}
mcmc.pars<-data_list$"mcmc.pars"
z<-data_list$"z"
K<-data_list$"K"
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances 
# mcmc.pars[,,3]: simulated weights
# the relabelling algorithm will run with the default initialization
# (no opt_init is specified)
run<-ecr.iterative.1(z = z, K = K)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances
# reordered.mcmc[,,3]: reordered weights
}
