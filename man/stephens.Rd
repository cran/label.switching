\name{stephens}
\alias{stephens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stephens' algorithm
}
\description{
Stephens (2000) developed a relabelling algorithm that makes the permuted sample points to agree as much as possible on the \eqn{n\times K} matrix of classification probabilities, using the Kullback-Leibler divergence. The algorithm's input is the matrix of allocation probabilities for each MCMC iteration.
}
\usage{
stephens(p, threshold, maxiter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
\eqn{m\times n \times K} dimensional array of allocation probabilities of the \eqn{n} observations among the \eqn{K} mixture components, for each iteration \eqn{t = 1,\ldots,m} of the MCMC algorithm.
}
  \item{threshold}{
An (optional) positive number controlling the convergence criterion. Default value: 1e-6.
}
  \item{maxiter}{
An (optional) integer controlling the max number of iterations. Default value: 100.
}
}
\details{
For a given MCMC iteration \eqn{t=1,\ldots,m}, let \eqn{w_k^{(t)}} and \eqn{\theta_k^{(t)}}, \eqn{k=1,\ldots,K} denote the simulated mixture weights and component specific parameters respectively. Then, the \eqn{(t,i,k)} element of \code{p} corresponds to the conditional probability that observation \eqn{i=1,\ldots,n} belongs to component \eqn{k} and is proportional to \eqn{p_{tik} \propto w_k^{(t)} f(x_i|\theta_k^{(t)}), k=1,\ldots,K}, where \eqn{f(x_i|\theta_k)} denotes the density of component \eqn{k}. This means that:
\deqn{p_{tik} = \frac{w_k^{(t)} f(x_i|\theta_k^{(t)})}{w_1^{(t)} f(x_i|\theta_1^{(t)})+\ldots + w_K^{(t)} f(x_i|\theta_K^{(t)})}.}
In case of hidden Markov models, the probabilities \eqn{w_k} should be replaced with the proper left (normalized) eigenvector of the state-transition matrix.
}
\value{
\item{permutations }{\eqn{m\times K} dimensional array of permutations}
\item{iterations }{integer denoting the number of iterations until convergence}
\item{status }{returns the exit status}
}
\references{
Stephens, M. (2000). Dealing with label Switching in mixture models. Journal of the Royal Statistical Society Series B, 62, 795-809.
}
\author{
Panagiotis Papastamoulis
}
\seealso{
\code{\link{permute.mcmc}}, \code{\link{label.switching}}
}
\examples{
#load a toy example: MCMC output consists of the random beta model
# applied to a normal mixture of \code{K=2} components. The number 
# of observations is equal to \code{n=5}. The number of MCMC samples
# is equal to \code{m=300}. The matrix of allocation probabilities 
# is stored to matrix \code{p}. 
data("mcmc_output")
# mcmc parameters are stored to array \code{mcmc.pars}
mcmc.pars<-data_list$"mcmc.pars"
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances 
# mcmc.pars[,,3]: simulated weights 
# the computed allocation matrix is p
p<-data_list$"p"
run<-stephens(p)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the components
# reordered.mcmc[,,2]: reordered variances
# reordered.mcmc[,,3]: reordered weights
}

