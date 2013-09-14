\name{sjw}
\alias{sjw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Probabilistic relabelling algorithm
}
\description{
Function to apply the probabilistic relabelling strategy of Sperrin et al (2010).
}
\usage{
sjw(mcmc.pars, z, complete, x, init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc.pars}{
%%     ~~Describe \code{mcmc.pars} here~~
\eqn{m\times K\times J} array of simulated MCMC parameters.
}
  \item{z}{
%%     ~~Describe \code{z} here~~
\eqn{m\times n} integer array of the latent allocation vectors generated from an MCMC algorithm.
}
  \item{complete}{
%%     ~~Describe \code{complete} here~~
function that returns the complete log-likelihood of the mixture model.  
}
  \item{x}{
%%     ~~Describe \code{x} here~~
\eqn{n}-dimensional data vector/array
}
  \item{init}{
%%     ~~Describe \code{init} here~~
An (optional) index pointing at the MCMC iteration whose parameters will initialize the algorithm. If it is less or equal to zero, the overall MCMC mean will be used for initialization. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \eqn{x=(x_1,\ldots,x_n)} denote the observed data and \eqn{\boldsymbol{w},\boldsymbol{\theta}} denote the mixture weights and component specific parameters, respectively. Assume that \eqn{K} is the the number of components. Then,
\deqn{L(\boldsymbol{w},\boldsymbol{\theta}|\boldsymbol{x})=\prod_{i=1}^{n}\sum_{k=1}^{K}w_k f_k(x_i|\theta_k),} 
\eqn{i=1,\ldots,n} is the observed likelihood of the mixture model. Given the latent allocation variables \eqn{\boldsymbol{z}=(z_1,\ldots,z_n)}, the complete likelihood of the model is defined as: 
\deqn{L_c(\boldsymbol{w},\boldsymbol{\theta}|\boldsymbol{x},\boldsymbol{z})=\prod_{i=1}^{n}w_{z_{i}}f_{z_{i}}(x_i|\theta_{z_{i}}).}
Then, \code{complete} corresponds to the log of \eqn{L_c} and should take as input the following: a vector of \eqn{n} allocations, the observed data and the parameters of the model as a \eqn{K\times J} array where \eqn{J} corresponds to the different parameter types of the model. See the example for an implementation at a univariate normal mixture.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{permutations }{\eqn{m\times K} dimensional array of permutations}
\item{iterations }{integer denoting the number of iterations until convergence}
}
\references{
%% ~put references to the literature/web site here ~
Sperrin M, Jaki T and Wit E (2010). Probabilistic relabelling strategies for the label switching problem in Bayesian mixture models.  Statistics and Computing, 20(3), 357-366.
}
\author{
%%  ~~who you are~~
Panagiotis Papastamoulis
}
\note{
%%  ~~further notes~~
This algorithm is not suggested for large number of components due to the computational overload: \eqn{K!} permutation probabilities are computed at each MCMC iteration. Moreover, the useR should carefully provide the complete log-likelihood function of the model as input to the algorithm and this makes its use quite complicated.
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
# equal to \code{m=300}.  
data("mcmc_output")
mcmc.pars<-data_list$"mcmc.pars"
z<-data_list$"z"
K<-data_list$"K"
x<-data_list$"x"

# mcmc parameters are stored to array \code{mcmc.pars}
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances
# mcmc.pars[,,3]: simulated weights 
# The number of different parameters for the univariate
# normal mixture is equal to J = 3: means, variances 
# and weights. The generated allocations variables are 
# stored to \code{z}. The observed data is stored to \code{x}.  
# The complete data log-likelihood is defined as follows:
complete.normal.loglikelihood<-function(x,z,pars){
#	x: data (size = n)
#	z: allocation vector (size = n)
#	pars: K\times J vector of normal mixture parameters:
#		pars[k,1] = mean of the k-normal component
#		pars[k,2] = variance of the k-normal component
#		pars[k,3] = weight of the k-normal component
#			k = 1,...,K
	g <- dim(pars)[1] #K (number of mixture components)
	n <- length(x)	#this denotes the sample size
	logl<- rep(0, n)	
 	logpi <- log(pars[,3])
	mean <- pars[,1]
	sigma <- sqrt(pars[,2])
	logl<-logpi[z] + dnorm(x,mean = mean[z],sd = sigma[z],log = TRUE)
	return(sum(logl))
}

#run the algorithm:
run<-sjw(mcmc = mcmc.pars,z = z, 
complete = complete.normal.loglikelihood,x = x, init=0)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two components
# reordered.mcmc[,,2]: reordered variances 
# reordered.mcmc[,,3]: reordered weights
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
