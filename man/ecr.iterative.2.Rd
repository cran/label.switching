\name{ecr.iterative.2}
\alias{ecr.iterative.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
ECR algorithm (iterative version 2)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function applies the second iterative version of Equivalence Classes Representatives (ECR) algorithm (Papastamoulis and Iliopoulos, 2010, Rodriguez and Walker, 2012). The set of all allocation variables is partitioned into equivalence classes and exactly one representative is chosen from each class. In this version the \eqn{m\times n \times K} of allocation probabilities should be given as input as well.   
}
\usage{
ecr.iterative.2(z, K, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{\eqn{m\times n} integer array of the latent allocation vectors generated from an MCMC algorithm.
%%     ~~Describe \code{z} here~~
}
  \item{K}{the number of mixture components (at least equal to 2).
%%     ~~Describe \code{K} here~~
}
  \item{p}{\eqn{m\times n \times K} dimensional array of allocation probabilities of the \eqn{n} observations among the \eqn{K} mixture components, for each iteration \eqn{t = 1,\ldots,m} of the MCMC algorithm.
%%     ~~Describe \code{p} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For a given MCMC iteration \eqn{t=1,\ldots,m}, let \eqn{w_k^{(t)}} and \eqn{\theta_k^{(t)}}, \eqn{k=1,\ldots,K} denote the simulated mixture weights and component specific parameters respectively. Then, the \eqn{(t,i,k)} element of \code{p} corresponds to the conditional probability that observation \eqn{i=1,\ldots,n} belongs to component \eqn{k} and is proportional to \eqn{p_{tik} \propto w_k^{(t)} f(x_i|\theta_k^{(t)}), k=1,\ldots,K}, where \eqn{f(x_i|\theta_k)} denotes the density of component \eqn{k}. This means that:
\deqn{p_{tik} = \frac{w_k^{(t)} f(x_i|\theta_k^{(t)})}{w_1^{(t)} f(x_i|\theta_1^{(t)})+\ldots + w_K^{(t)} f(x_i|\theta_K^{(t)})}.}
In case of hidden Markov models, the probabilities \eqn{w_k} should be replaced with the proper left (normalized) eigenvector of the state-transition matrix.
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
Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.

Rodriguez C.E. and Walker S. (2012). Label Switching in Bayesian Mixture Models: Deterministic relabeling strategies. Journal of Computational and Graphical Statistics. DOI:10.1080/10618600.2012.735624
}
\author{
%%  ~~who you are~~
Panagiotis Papastamoulis
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{permute.mcmc}}, \code{\link{label.switching}}, \code{\link{ecr}}, \code{\link{ecr.iterative.1}}, \code{\link{stephens}} 
}
\examples{
#load a toy example: MCMC output consists of the random beta model
# applied to a normal mixture of \code{K=2} components. The number of
# observations is equal to \code{n=5}. The number of MCMC samples is
# equal to \code{m=1000}. The 300 simulated allocations are stored to
# array \code{z}. The matrix of allocation probabilities is stored to
# array \code{p}. 
data("mcmc_output")
z<-data_list$"z"
K<-data_list$"K"
p<-data_list$"p"
# mcmc parameters are stored to array \code{mcmc.pars}
mcmc.pars<-data_list$"mcmc.pars"
# mcmc.pars[,,1]: simulated means of the two components
# mcmc.pars[,,2]: simulated variances 
# mcmc.pars[,,3]: simulated weights
# the relabelling algorithm will run with the default initialization
# (no opt_init is specified)
run<-ecr.iterative.2(z = z, K = 2, p = p)
# apply the permutations returned by typing:
reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
# reordered.mcmc[,,1]: reordered means of the two mixture components
# reordered.mcmc[,,2]: reordered variances of the two components
# reordered.mcmc[,,3]: reordered weights of the two components
}
