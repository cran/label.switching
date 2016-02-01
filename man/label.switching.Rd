\name{label.switching}
\alias{label.switching}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main calling function
}
\description{
This is the main function of the package. It is used to reorder a simulated MCMC sample of the parameters of a mixture (or more general a hidden Markov model) according to eight label switching solving methods: ECR algorithm (default version), ECR algorithm (two iterative versions), PRA algorithm, Stephens' algorithm, Artificial Identifiability Constraint (AIC), Data-Based relabelling and a probabilistic relabelling algorithm (SJW). The input depends on the type of the label switching method. The output contains a list with the permutation returned by each method, the corresponding single best clusterings and the CPU time demanded for each method. In what follows: \eqn{m} denotes the number of MCMC iterations, \eqn{n} denotes the sample size of the observed data, \eqn{K} denotes the number of mixture components and \eqn{J} the number of different types of parameters of the model.
}
\usage{
label.switching(method, zpivot, z, K, prapivot, p, complete, 
			mcmc, sjwinit, data, constraint, 
			groundTruth, thrECR, thrSTE, thrSJW, 
			maxECR, maxSTE, maxSJW, userPerm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
%%     ~~Describe \code{method} here~~
any non-empty subset of c("ECR","ECR-ITERATIVE-1","PRA","ECR-ITERATIVE-2","STEPHENS","SJW","AIC","DATA-BASED") indicating the desired label-switching solving method. Also available is the option "USER-PERM" which corresponds to a user-defined set of permutations \code{userPerm}.
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
\eqn{n}-dimensional data vector/array. Needed by the \code{\link{sjw}} and \code{\link{dataBased}} algorithms.
}
  \item{constraint}{
An (optional) integer between 1 and J corresponding to the parameter that will be used to apply the Identifiabiality Constraint. In this casethe mcmc output is reordered according to the constraint \eqn{mcmc[i,1,constraint] < \ldots < mcmc[i,K,constraint]}. If \code{constraint = "ALL"}, all \eqn{J} Identifiability Constraints are applied. Default value: 1.
}
  \item{groundTruth}{
Optional integer vector of \eqn{n} allocations, which are considered as the 'ground truth' allocations of the \eqn{n} observations among the \eqn{K} mixture components. The output of all methods will be relabelled in a way that the resulting single best clusterings maximize their similarity with the ground truth. This option is very useful in simulation studies or in any other case that the cluster labels are known in order to perform comparisons between methods. 
}
  \item{thrECR}{
An (optional) positive number controlling the convergence criterion for \code{ecr.iterative.1} and \code{ecr.iterative.2}. Default value: 1e-6.
}
  \item{thrSTE}{
An (optional) positive number controlling the convergence criterion for \code{stephens}. Default value: 1e-6.
}
  \item{thrSJW}{
An (optional) positive number controlling the convergence criterion for \code{sjw}. Default value: 1e-6.
}
  \item{maxECR}{
An (optional) integer controlling the max number of iterations for \code{ecr.iterative.1} and \code{ecr.iterative.2}. Default value: 100.
}
  \item{maxSTE}{
An (optional) integer controlling the max number of iterations for \code{stephens}. Default value: 100.
}
  \item{maxSJW}{
An (optional) integer controlling the max number of iterations for \code{sjw}. Default value: 100.
}
  \item{userPerm}{
An (optional) list with user-defined permutations. It is required only if "USER-PERM" has been chosen in \code{method}. In this case, \code{userPerm[[i]]} is an \eqn{m\times K} array of permutations for all \eqn{i = 1,\ldots,S}, where \eqn{S} denotes the number of permutation arrays. This is useful in case that the user wants to compare his/hers own relabelling method with the available ones. 
}
}
\value{
\item{permutations }{an \eqn{m\times K} array of permutations per method.}
\item{clusters }{an \eqn{n} dimensional vector of best clustering of the the observations for each method.}
\item{timings }{CPU time needed for each relabelling method.}
\item{similarity }{correlation matrix between the label switching solving methods in terms of their matching best-clustering allocations.}
}
\author{
Panagiotis Papastamoulis
}
\note{
If the ground truth is not given, all methods are reordered using the estimated single best clustering of the first provided method. The methods \code{\link{sjw}} and \code{\link{pra}} are not suggested for large number of components. Also note that \code{\link{sjw}} might be quite slow even for small number of components. In this case try adjusting \code{thrSJW} or \code{maxSJW} to smaller values the default ones.
}
\references{
Papastamoulis P. (2016). label.switching: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. Journal of Statistical Software, Code Snippets, 69(1): 1-24.

}


\seealso{
\code{\link{ecr}}, \code{\link{ecr.iterative.1}}, \code{\link{ecr.iterative.2}}, \code{\link{stephens}}, \code{\link{pra}}, \code{\link{sjw}}, \code{\link{dataBased}}, \code{\link{aic}}
}
\examples{
# We will apply four methods:
# ECR, ECR-ITERATIVE-1, PRA, AIC, STEPHENS and DATA-BASED.
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
x<-data_list$"x"
mapindex<-data_list$"mapindex"
mapindex.non<-data_list$"mapindex.non"
# The PRA method will use as pivot the iteration that corresponds to
# the non-complete MAP estimate (mapindex). 

#Apply the six methods by typing:

ls<-label.switching(method=c("ECR","ECR-ITERATIVE-1","PRA","STEPHENS","AIC","DATA-BASED"),
zpivot=z[c(mapindex,mapindex.non),],z = z,K = K, data = x,
prapivot = mcmc.pars[mapindex,,],p=p,mcmc = mcmc.pars)

#plot the raw and reordered means of the K=2 normal mixture components for each method
par(mfrow = c(2,4))
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
matplot(permute.mcmc(mcmc.pars,ls$permutations$"AIC")$output[,,1],type="l",
xlab="iteration",main="AIC",ylab = "means")
matplot(permute.mcmc(mcmc.pars,ls$permutations$"DATA-BASED")$output[,,1],type="l",
xlab="iteration",main="DATA-BASED",ylab = "means")

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
#"PRA","STEPHENS","SJW","AIC","DATA-BASED"),
#zpivot=z[c(mapindex,mapindex.non),],z = z,
#K = K,prapivot = mcmc.pars[mapindex,,],p=p,
#complete = complete.normal.loglikelihood,mcmc.pars,
#data = x)

}
