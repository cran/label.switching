\name{dataBased}
\alias{dataBased}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data-based labelling
}
\description{
This function reorders the MCMC output according the data-based relabelling algorithm of Rodriguez and Walker (2014). The idea is to define a loss function which resembles a k-means type diveging measure of cluster centers. After the cluster centers have been estimated, the algorithm finds the optimal permutations that switch every simulated MCMC sample to them.
}
\usage{
dataBased(x, K, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\eqn{n}-dimensional data vector/array. 
}
  \item{K}{
the number of mixture components.
}
  \item{z}{
\eqn{m\times n} integer array of the latent allocation vectors generated from an MCMC algorithm. 
}
}
\value{
\item{permutations }{\eqn{m\times K} dimensional array of permutations}
}
\references{
Rodriguez C.E. and Walker S. (2014). Label Switching in Bayesian Mixture Models: Deterministic relabeling strategies. Journal of Computational and Graphical Statistics. 23:1, 25-45
}
\author{
Panagiotis Papastamoulis
}
\seealso{
\code{\link{permute.mcmc}}, \code{\link{label.switching}}
}
\examples{
#load a toy example: MCMC output consists of the random beta model
# applied to a normal mixture of \code{K=2} components. The number of
# observations is equal to \code{n=5}. The number of MCMC samples is
# equal to \code{m=300}. The 1000 generated MCMC samples are stored 
#to array mcmc.pars. 
data("mcmc_output")
z<-data_list$"z"
K<-data_list$"K"
x<-data_list$"x"
mcmc.pars<-data_list$"mcmc.pars"
# mcmc parameters are stored to array \code{mcmc.pars}
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances of the two components
# mcmc.pars[,,3]: simulated weights of the two components
# Apply dataBased relabelling
run<-dataBased(x = x, K = K, z = z)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances of the components
# reordered.mcmc[,,3]: reordered weights 
}

