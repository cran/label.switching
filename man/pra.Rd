\name{pra}
\alias{pra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PRA algorithm
}
\description{
This function reorders the MCMC output using the geometrically-based Pivotal Reordering Algorithm (PRA) (Marin et al, 2005, Marin and Robert, 2007). The method requires as input the generated MCMC sample and a pivot parameter vector. The user should be careful in order the pivot elements have the same parameters with the generated MCMC output. The simulated MCMC sample should be provided by the useR as a \eqn{m\times K\times J} dimensional array, where \eqn{m} denotes the number of MCMC samples, \eqn{K} denotes the number of mixture components and \eqn{J} corresponds to the number of different parameter types of the model. The pivot should correspond to a high-posterior density point.
}
\usage{
pra(mcmc.pars, pivot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc.pars}{
%%     ~~Describe \code{mcmc.pars} here~~
\eqn{m\times K\times J} array of simulated MCMC parameters.
}
  \item{pivot}{
%%     ~~Describe \code{pivot} here~~
\eqn{K\times J} array containing the parameter that will be used as a pivot.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The positive integer \eqn{J} denotes the number of different parameter types of the model. For example, in a univariate normal mixture model there are \eqn{J = 3} different types: means, variances and weights. In a Poisson mixture there are \eqn{J=2} types: means and weights.
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
%% ~put references to the literature/web site here ~
Marin, J.M., Mengersen, K. and Robert, C.P. (2005). Bayesian modelling and inference on mixtures of distributions. Handbook of Statistics (25), D. Dey and C.R. Rao (eds). Elsevier-Sciences.

Marin, J.M. and Robert, C.P. (2007). Bayesian Core: A Practical Approach to Computational Bayesian Statistics, Springer-Verlag, New York.
}
\author{
Panagiotis Papastamoulis
}
\note{
%%  ~~further notes~~
This algorithm is very fast for small number of components. It is not suggested for large number of components due to the calculation of \eqn{K!} distances for each MCMC sample.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
mapindex<-data_list$"mapindex"

# mcmc parameters are stored to array \code{mcmc.pars}
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances of the two components
# mcmc.pars[,,3]: simulated weights of the two components
# We will apply PRA using as pivot the complete MAP estimate
# which corresponds to \code{mcmc.pars[mapindex,,]}
run<-pra(mcmc = mcmc.pars, pivot = mcmc.pars[mapindex,,])
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances of the components
# reordered.mcmc[,,3]: reordered weights 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
