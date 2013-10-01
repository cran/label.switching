\name{label.switching}
\alias{label.switching}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Main calling function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is the main function of the package. It is used to reorder a simulated MCMC sample of the parameters of a mixture (or more general a hidden Markov model) according to six label switching solving methods: ECR algorithm (default version), ECR algorithm (two iterative versions), PRA algorithm, Stephens' algorithm and a probabilistic relabelling algorithm. The input depends on the type of the label switching method. The output contains a list with the permutation returned by each method, the corresponding single best clusterings and the CPU time demanded for each method. In what follows: \eqn{m} denotes the number of MCMC iterations, \eqn{n} denotes the sample size of the observed data, \eqn{K} denotes the number of mixture components and \eqn{J} the number of different types of parameters of the model.
}
\usage{
label.switching(method, zpivot, z, K, prapivot, p, complete, 
			mcmc, sjwinit, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
%%     ~~Describe \code{method} here~~
any non-empty subset of c("ECR","ECR-ITERATIVE-1","PRA","ECR-ITERATIVE-2","STEPHENS","SJW") indicating the desired label-switching solving method.
}
  \item{zpivot}{
%%     ~~Describe \code{zpivot} here~~
\eqn{d\times n}-dimensional array of pivot allocation vectors, where \eqn{d} denotes the number of pivots. This is demanded by the \code{\link{ecr}} method. The method will be applied \eqn{d} times.
}
  \item{z}{
%%     ~~Describe \code{z} here~~
\eqn{m\times n} integer array of the latent allocation vectors generated from an MCMC algorithm. 
}
  \item{K}{
%%     ~~Describe \code{K} here~~
the number of mixture components. This is demanded by the \code{\link{ecr}}, \code{\link{ecr.iterative.1}} and \code{\link{ecr.iterative.2}} methods.
}
  \item{prapivot}{
%%     ~~Describe \code{prapivot} here~~
\eqn{K\times J} array containing the parameter that will be used as a pivot by the \code{\link{pra}} method.
}
  \item{p}{
%%     ~~Describe \code{p} here~~
\eqn{m\times n \times K} dimensional array of allocation probabilities of the \eqn{n} observations among the \eqn{K} mixture components, for each iteration \eqn{t = 1,\ldots,m} of the MCMC algorithm. This is demanded by the \code{\link{ecr.iterative.2}} and \code{\link{stephens}} methods.
}
  \item{complete}{
%%     ~~Describe \code{complete} here~~
function that returns the complete log-likelihood of the mixture model. Demanded by \code{\link{sjw}} method.
}
  \item{mcmc}{
%%     ~~Describe \code{mcmc} here~~
\eqn{m\times K\times J} array of simulated MCMC parameters. Needed by \code{\link{sjw}} and \code{\link{pra}} methods.
}
  \item{sjwinit}{
%%     ~~Describe \code{sjwinit} here~~
An index pointing at the MCMC iteration whose parameters will initialize the \code{\link{sjw}} algorithm (optional). 
}
  \item{data}{
%%     ~~Describe \code{data} here~~
\eqn{n}-dimensional data vector/array. Needed by the \code{\link{sjw}} algorithm.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The methods \code{\link{sjw}} and \code{\link{pra}} are not suggested for large number of components. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{permutations }{an \eqn{m\times K} array of permutations per method.}
\item{clusters }{an \eqn{n} dimensional vector of best clustering of the the observations for each method.}
\item{timings }{CPU time needed for each relabelling method.}
\item{similarity }{correlation matrix between the label switching solving methods in terms of their matching best-clustering allocations.}
}
\author{
%%  ~~who you are~~
Panagiotis Papastamoulis
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ecr}}, \code{\link{ecr.iterative.1}}, \code{\link{ecr.iterative.2}}, \code{\link{stephens}}, \code{\link{pra}}, \code{\link{sjw}} 
}
\examples{
# We will apply four methods:
# ECR, ECR-ITERATIVE-1, PRA and STEPHENS.
# default ECR will use two different pivots.

#load a toy example: MCMC output consists of the random beta model
# applied to a normal mixture of \code{K=2} components. The number of
# observations is equal to \code{n=5}. The number of MCMC samples is
# equal to \code{m=300}. simulated allocations are stored to array \code{z}. 
data("mcmc_output")
mcmc.pars<-data_list$"mcmc.pars"
# mcmc parameters are stored to array \code{mcmc.pars}
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances 
# mcmc.pars[,,3]: simulated weights 
# We will use two pivots for default ECR algorithm:
# the first one corresponds to iteration \code{mapindex} (complete MAP)
# the second one corresponds to iteration \code{mapindex.non} (observed MAP)
# The array \code{p} corresponds to the the allocation probabilities 
z<-data_list$"z"
K<-data_list$"K"
p<-data_list$"p"
mapindex<-data_list$"mapindex"
mapindex.non<-data_list$"mapindex.non"
# The PRA method will use as pivot the iteration that corresponds to
# the non-complete MAP estimate (mapindex). 

#Apply the four methods by typing:

ls<-label.switching(method=c("ECR","ECR-ITERATIVE-1","PRA","STEPHENS"),
zpivot=z[c(mapindex,mapindex.non),],z = z,K = K,
prapivot = mcmc.pars[mapindex,,],p=p,mcmc = mcmc.pars)

#plot the raw and reordered means of the K=2 normal mixture components for each method
par(mfrow = c(2,3))
#raw MCMC output for the means (with label switching)
matplot(mcmc.pars[,,1],type="l",
xlab="iteration",main="Raw MCMC output",ylab = "means")
# Reordered outputs
matplot(permute.mcmc(mcmc.pars,ls$permutations$"ECR-1")$output[,,1],type="l",
xlab="iteration",main="ECR (1st pivot)",ylab = "means")
matplot(permute.mcmc(mcmc.pars,ls$permutations$"ECR-2")$output[,,1],type="l",
xlab="iteration",main="ECR (2nd pivot)",ylab = "means")
matplot(permute.mcmc(mcmc.pars,ls$permutations$"ECR-ITERATIVE-1")$output[,,1],
type="l",xlab="iteration",main="ECR-iterative-1",ylab = "means")
matplot(permute.mcmc(mcmc.pars,ls$permutations$"PRA")$output[,,1],type="l",
xlab="iteration",main="PRA",ylab = "means")
matplot(permute.mcmc(mcmc.pars,ls$permutations$"STEPHENS")$output[,,1],type="l",
xlab="iteration",main="STEPHENS",ylab = "means")

#######################################################
# if the useR wants to apply the SJW algorithm as well:
# The SJW method needs to define the complete log-likelihood of the
# model. For the univariate normal mixture, this is done as follows:

complete.normal.loglikelihood<-function(x,z,pars){
	#x: denotes the n data points
	#z: denotes an allocation vector (size=n)
	#pars: K\times 3 vector of means,variance, weights
	# pars[k,1]: corresponds to the mean of component k
	# pars[k,2]: corresponds to the variance of component k
	# pars[k,3]: corresponds to the weight of component k
	g <- dim(pars)[1]
	n <- length(x)
	logl<- rep(0, n)
 	logpi <- log(pars[,3])
	mean <- pars[,1]
	sigma <- sqrt(pars[,2])
	logl<-logpi[z] + dnorm(x,mean = mean[z],sd = sigma[z],log = T)
	return(sum(logl))
}

# and then run (after removing all #):
#ls<-label.switching(method=c("ECR","ECR-ITERATIVE-1","ECR-ITERATIVE-2",
#"PRA","STEPHENS","SJW"),zpivot=z[c(mapindex,mapindex.non),],z = z,
#K = K,prapivot = mcmc.pars[mapindex,,],p=p,
#complete = complete.normal.loglikelihood,mcmc.pars,
#data = x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
