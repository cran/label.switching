\name{label.switching-package}
\alias{label.switching-package}
\docType{package}
\title{
Algorithms for solving the label switching problem
}
\description{
This package can be used to reorder MCMC outputs of parameters of mixture models (or more general ones, like hidden Markov). The label switching phenomenon is a fundamental problem to MCMC estimation of the parameters of such models. This package contains eight label switching solving algorithms: the default and iterative versions of ECR algorithm (Papastamoulis and Iliopoulos, 2010, 2013, Rodriguez and Walker, 2014, Papastamoulis, 2014), the data-based algorithm (Rodriguez and Walker, 2014), the Kullback-Leibler based algorithm of Stephens (2000), the probabilistic relabelling algorithm of Sperrin et al (2010), the artificial identifiability constraints method and the PRA algorithm (Marin et al, 2005, Marin and Robert, 2007). The user input depends on each method. Each algorithm returns a list of permutations. For comparison purposes, the user can also provide his/hers own set of permutations.
}
\details{
\tabular{ll}{
Package: \tab label.switching\cr
Type: \tab Package\cr
Version: \tab 1.8\cr
Date: \tab 2019-07-01\cr
License: \tab GPL-2\cr
}
This is NOT a package to simulate MCMC samples from the posterior distribution of mixture models. MCMC output and related information serves as input to the available methods. There are eight functions that can be used to post-process the MCMC output:
\tabular{lll}{
Function \tab Method \tab Input \cr
---------------------------- \tab ---------------------------- \tab ----------------------------------- \cr
\code{aic} \tab ordering constraints \tab \code{mcmc,constraint}\cr
\code{dataBased} \tab data based \tab \code{x,K,z}\cr
\code{ecr} \tab ECR (default) \tab \code{zpivot, z, K}\cr
\code{ecr.iterative.1} \tab ECR (iterative vs. 1)\tab \code{z, K}\cr
\code{ecr.iterative.2} \tab ECR (iterative vs. 2)\tab\code{z, K, p}\cr
\code{pra} \tab PRA \tab \code{mcmc, pivot} \cr
\code{stephens} \tab Stephens\tab \code{p}\cr
\code{sjw} \tab Probabilistic\tab \code{mcmc, z, complete, x}
}
Each function returns an \eqn{m\times K} array of permutations, where \eqn{m} and \eqn{K} denote the MCMC sample size and number of mixture components, respectively. Next, these permutations can be applied to reorder the MCMC sample by applying the function \code{permute.mcmc}. The useR can call any of the above functions simultaneously using the main function of the package: \code{\link{label.switching}}. 
}
\note{
The most common method is to impose an identifiability constraint \code{aic}, however this approach has been widely criticized in the literature. The methods \code{ecr, ecr.iterative.1,ecr.iterative.2}, \code{stephens, dataBased} are solving the label switching problem using the function \code{lpAssign} of the package \code{lpSolve}. This is an integer programming algorithm for the solution of the assignment problem. Hence, these functions are computationally efficient even in cases where the number of components is quite large. On the other hand, methods \code{pra} and \code{sjw} are not designed in this way, so they are not suggested for large \eqn{K}. 
}

\author{
Panagiotis Papastamoulis

Maintainer: <papapast@yahoo.gr>
}
\references{
Marin, J.M., Mengersen, K. and Robert, C.P. (2005). Bayesian modelling and inference on mixtures of distributions. Handbook of Statistics (25), D. Dey and C.R. Rao (eds). Elsevier-Sciences.

Marin, J.M. and Robert, C.P. (2007). Bayesian Core: A Practical Approach to Computational Bayesian Statistics, Springer-Verlag, New York.

Papastamoulis P. and Iliopoulos G. (2010). An artificial allocations based solution to the label switching problem in Bayesian analysis of mixtures of distributions. Journal of Computational and Graphical Statistics, 19: 313-331.

Papastamoulis P. and Iliopoulos G. (2013). On the convergence rate of Random Permutation Sampler and ECR algorithm in missing data models. Methodology and Computing in Applied Probability, 15(2): 293-304.

Papastamoulis P. (2014). Handling the label switching problem in latent class models via the ECR algorithm. Communications in Statistics, Simulation and Computation, 43(4): 913-927.

Papastamoulis P. (2016). label.switching: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. Journal of Statistical Software, Code Snippets, 69(1): 1-24.

Rodriguez C.E. and Walker S. (2014). Label Switching in Bayesian Mixture Models: Deterministic relabeling strategies. Journal of Computational and Graphical Statistics. 23:1, 25-45

Sperrin M, Jaki T and Wit E (2010). Probabilistic relabelling strategies for the label switching problem in Bayesian mixture models. Statistics and Computing, 20(3), 357-366.

Stephens, M. (2000). Dealing with label Switching in mixture models. Journal of the Royal Statistical Society Series B, 62, 795-809.
}
\keyword{ package }
\seealso{
\code{\link{label.switching}}
}
