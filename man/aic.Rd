\name{aic}
\alias{aic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Artificial Identifiability Constraints
}
\description{
This function relabels the MCMC output by simply ordering a specific parameter. Let \eqn{m}, \eqn{K} and \eqn{J} denote the number of simulated MCMC samples, number of mixture components and different parameter types, respectively.
}
\usage{
aic(mcmc.pars, constraint)
}
\arguments{
  \item{mcmc.pars}{
\eqn{m\times K\times J} array of simulated MCMC parameters.
}
  \item{constraint}{
An integer between 1 and J corresponding to the parameter that will be used to apply the Identifiabiality Constraint. In this case, the MCMC output is reordered according to the constraint \deqn{mcmc.pars[i,1,constraint] < \ldots < mcmc.pars[i,K,constraint],} for all \eqn{i=1,\ldots,m}. If \code{constraint = "ALL"}, all \eqn{J} Identifiability Constraints are applied. 
}
}
\value{
\item{permutations }{an \eqn{m\times K} array of permutations.}
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
mcmc.pars<-data_list$"mcmc.pars"

# mcmc parameters are stored to array \code{mcmc.pars}
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances of the two components
# mcmc.pars[,,3]: simulated weights of the two components
# We will apply AIC by ordering the means
# which corresponds to value \code{constraint=1}
run<-aic(mcmc = mcmc.pars,constraint=1)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances of the components
# reordered.mcmc[,,3]: reordered weights 
}

