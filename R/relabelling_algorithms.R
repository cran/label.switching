# ECR algorithm:
# zpivot: the pivot
# z: mcmc output of allocation variables with dimension equal to mxn (m denotes MCMC iterations)
ecr<-function(zpivot,z,K){
if(K<max(z)){stop("K should be at least equal to max(z)")}
if(dim(z)[2]!=length(zpivot)){if(K<max(z)){stop("zpivot has not the same length with simulated z's")}}
n<-length(zpivot)
m<-dim(z)[1]
k<-K
# ECR algorithm
# By considering the problem of maxizing the S similarity between z[,iter] and zpivot
# as a special case of the assignment problem, it minimizes the corresponding cost matrix 
# between z[,iter] and zpivot
#library(lpSolve);				# contains the routine lp.assign for the solution of the assignment problem
st <- 1:k;							# the set {1,...,k}
perm <- array(data = NA, dim=c(m,K));	# the permutation of st that will be used to reorder the parameter values at each iteration
cost.matrix <- matrix(numeric(k*k),nrow=k,ncol=k); # (k times k) cost matrix of the assignment problem
s <- 1:n
#ptm<-proc.time()
for(iter in 1:m){
	alloc <- z[iter,]
	for(i in 1:k){
		so <- zpivot[s[alloc==i]];	# finding the indices of zpivot that correspond to z[,iter]==i
		l <- length(alloc[alloc==i]);	# the length of the vector above
		for(j in 1:k) cost.matrix[i,j] <- l-length(so[so==j])		# the cost of assigning to index i the permutation j
		}
	matr <- lp.assign(cost.matrix)$solution ;	# solution: matr[j,i] = 1 <=> index i assigned to index j
	for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]}	# the optimal permutation for the current iteration
	}
	#time<-proc.time() - ptm
#	print(paste("ECR algorithm time:", as.numeric(time[3]),"seconds."))
	results<-list(perm)
	names(results)<-c("permutations")
	return(results)
}
#################################################################################
#################################################################################
ecr.iterative.1<-function (z, K, opt_init) {
    if (K < max(z)) {
        stop("K should be at least equal to max(z)")
    }
    m <- dim(z)[1]
    k <- K
    st <- 1:k
    perm <- array(data = NA, dim = c(m, K))
    cost.matrix <- matrix(numeric(k * k), nrow = k, ncol = k)
    n <- dim(z)[2]
    s <- 1:n
    if (missing(opt_init)) {
        for (j in 1:K) {
            perm[, j] <- j
        }
    }
    else {
        perm <- opt_init
    }
    permz <- z
    for(iter in 1:m){
	permz[iter, ] <- perm[iter, ][z[iter, ]]
    }	
    zpivot <- numeric(n)
    criterion <- 99
    threshold <- 10^(-6)
    t <- 1
    zpivot <- apply(permz, 2, function(y) order(table(y), decreasing = T)[1])
    cf <- 0
    for (iter in 1:m) {
        alloc <- z[iter, ]
        for (i in 1:k) {
            so <- zpivot[s[alloc == i]]
            l <- length(alloc[alloc == i])
            for (j in 1:k) cost.matrix[i, j] <- l - length(so[so == 
                j])
        }
        matr <- lp.assign(cost.matrix)$solution
        for (i in 1:k) {
            perm[iter, i] <- st[matr[, i] > 0]
        }
        cf <- cf + sum(cost.matrix * matr)
        permz[iter, ] <- perm[iter, ][z[iter, ]]
    }
    previous <- cf
    #print(paste("t = ",t,"cost function = ", cf))
    while (criterion > threshold) {
	previousperm<-perm
  #   while (t < 4) {
        t <- t + 1
        zpivot <- apply(permz, 2, function(y) order(table(y), 
            decreasing = T)[1])
        cf <- 0
        for (iter in 1:m) {
            alloc <- z[iter, ]
            for (i in 1:k) {
                so <- zpivot[s[alloc == i]]
                l <- length(alloc[alloc == i])
                for (j in 1:k) cost.matrix[i, j] <- l - length(so[so == 
                  j])
            }
            matr <- lp.assign(cost.matrix)$solution
            for (i in 1:k) {
                perm[iter, i] <- st[matr[, i] > 0]
            }
            cf <- cf + sum(cost.matrix * matr)
            permz[iter, ] <- perm[iter, ][z[iter, ]]
        }
        current <- cf
        criterion <- previous - current
        previous <- cf
	if (criterion<0){
		perm<-previousperm
	}
        #print(paste("t = ",t,"cost function = ", cf))
    }
    results <- list(perm, t)
    names(results) <- c("permutations", "iterations")
    return(results)
}









###################################################################################
#######################################################################################
#######################################################################################
#######################################################################################





ecr.iterative.2<-function(z,K,p){
if(K<max(z)){stop("K should be at least equal to max(z)")}
n<-dim(z)[2]
m<-dim(z)[1]
k<-K
#library(lpSolve);				# contains the routine lp.assign for the solution of the assignment problem
st <- 1:k;							# the set {1,...,k}
perm <- array(data = NA, dim=c(m,K));	# the permutation of st that will be used to reorder the parameter values at each iteration
cost.matrix <- matrix(numeric(k*k),nrow=k,ncol=k); # (k times k) cost matrix of the assignment problem
s <- 1:n
#ptm<-proc.time()

#initial permutations: identity
for (j in 1:K){
perm[,j]<-j
}


zpivot<-numeric(n)

criterion<-99
threshold<-10**(-6)
t<-1
#estimating pivot

q<-array(data = 0, dim =c(n,k))
for (j in 1:k){
	for(iter in 1:m){q[,j]<-q[,j] + p[iter,,perm[iter,j]]}
}
q<-q/m

for(i in 1:n){zpivot[i]<-order(q[i,])[K]}
cf<-0
for(iter in 1:m){
	alloc <- z[iter,]
	for(i in 1:k){
		so <- zpivot[s[alloc==i]];	# finding the indices of zpivot that correspond to z[,iter]==i
		l <- length(alloc[alloc==i]);	# the length of the vector above
		for(j in 1:k)cost.matrix[i,j] <- l-length(so[so==j])
	}
	matr <- lp.assign(cost.matrix)$solution ;	# solution: matr[j,i] = 1 <=> index i assigned to index j
	for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]}	# the optimal permutation for the current iteration
	cf<-cf + sum(cost.matrix*matr)
}
previous<-cf
#print(paste("iteration 1, cost function =",cf))

while(criterion>threshold){
	t<-t + 1
	#estimating pivot
	q<-array(data = 0, dim =c(n,k))
	for (j in 1:k){
		for(iter in 1:m){q[,j]<-q[,j] + p[iter,,perm[iter,j]]}
	}
	q<-q/m

	for(i in 1:n){zpivot[i]<-order(q[i,])[K]}
	cf<-0
	for(iter in 1:m){
		alloc <- z[iter,]
		for(i in 1:k){
			so <- zpivot[s[alloc==i]];	# finding the indices of zpivot that correspond to z[,iter]==i
			l <- length(alloc[alloc==i]);	# the length of the vector above
			for(j in 1:k)cost.matrix[i,j] <- l-length(so[so==j])
		}
		matr <- lp.assign(cost.matrix)$solution ;	# solution: matr[j,i] = 1 <=> index i assigned to index j
		for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]}	# the optimal permutation for the current iteration
		cf<-cf + sum(cost.matrix*matr)
	}
	current<-cf
	criterion<-abs(previous - current)
	previous<-cf
	#print(paste("iteration", t, "criterion =",cf))


}
	#time<-proc.time() - ptm
	#print(paste("Iterative ECR algorithm 2 converged at", t,"iterations"))
	#print(paste("Iterative ECR algorithm 2 time:", as.numeric(time[3]),"seconds."))
	results<-list(perm,t)
	names(results)<-c("permutations","iterations")
	return(results)
}












###################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
# p: mxnxK matrix of allocation proababilities

stephens<-function(p){
n<-dim(p)[2]
m<-dim(p)[1]
k<-dim(p)[3]
burnin<-0
K<-k
#library(lpSolve);				# contains the routine lp.assign for the solution of the assignment problem
st <- 1:k;							# the set {1,...,k}
perm <- c(numeric(k));	# the permutation of st that will be used to reorder the parameter values at each iteration
cost.matrix <- matrix(numeric(k*k),nrow=k,ncol=k); # (k times k) cost matrix of the assignment problem
s <- 1:n
###

up.threshold<-1-10**(-6)
down.threshold<-10**(-6)
for(k in 1:K){
for(i in 1:n){
up.index<-which(p[,i,k]>up.threshold)
down.index<-which(p[,i,k]<down.threshold)
if(length(up.index)>0){p[up.index,i,k]<-rep(up.threshold,length(up.index))}
if(length(down.index)>0){p[down.index,i,k]<-rep(down.threshold,length(down.index))}
}
}
for(iter in 1:m){
p[iter,,]<-p[iter,,]/rowSums(p[iter,,])
}
##
#print(paste("finished smoothing"))
#step 0.
perm <- array(data = NA, dim = c(m,k))
for(j in 1:k){perm[,j]<-j}
q<-array(data = 0, dim =c(n,k))
previous<- -99
criterion<-99
threshold<-10**(-6)
maxiter<-99
#ptm<-proc.time()
# t denotes stephens algorithm iterations
t<-0
while((criterion>threshold)&&(t<maxiter)){
	t<-t+1
#compute q matrix
	q<-array(data = 0, dim =c(n,k))
	for (j in 1:k){
		for(iter in 1:m){q[,j]<-q[,j] + p[iter,,perm[iter,j]]}
	}
	q<-q/m
	
	for(iter in 1:m){
		for(j in 1:k){		
			temp<-p[iter,,]*(log(p[iter,,]) - log(q[,j]))
			cost.matrix[j,]<-colSums(temp)
		}
		matr <- lp.assign(cost.matrix)$solution
		for(i in 1:k){perm[iter,i] <- st[matr[,i]>0]
		}
		perm[iter,]<-order(perm[iter,])
	}
	current<-cost.function<-sum(cost.matrix)
	criterion<-abs(previous - current)
	previous<-current
	#if(t>1){print(paste("iteration ", t,", criterion =",criterion))}
#this ends t
#
}
	#if(t<maxiter){print(paste("Stephens' algorithm converged at", t,"iterations"))}else{
	#	print(paste("Stephens' algorithm didn't converged at", maxiter,"iterations"))
	#}
	if(t>maxiter){print(paste("Stephens' algorithm didn't converged at", maxiter,"iterations"))}
	#time<-proc.time() - ptm
	#print(paste("Stephens' algorithm time:", as.numeric(time[3]),"seconds."))
	results<-list(perm,t)
	names(results)<-c("permutations","iterations")
	return(results)
}
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


sjw<-function(mcmc.pars,z,complete,x,init){
K<-dim(mcmc.pars)[2]
J<-dim(mcmc.pars)[3]
m<-dim(mcmc.pars)[1]

ptm <- proc.time()

#initial estimate
# if init = 0 algorithm starts from the overall mean.
if(init <1){
#print(paste("Initialization using the overall mean"))
theta.hat<-array(data = NA, dim =c(K,J))
for(k in 1:K){
for(j in 1:J){
theta.hat[k,j]<-mean(mcmc.pars[,k,j])
}
}
}else{
init = floor(init + 0.5)
#print(paste("Initialization using MCMC iteration: ",init))
theta.hat<-mcmc.pars[init,,]}
#print(paste(theta.hat))
previous<-theta.hat
#library(combinat)
permutations<-permn(1:K)
kfactorial<-gamma(K+1)
#permutation probabilities
pp<-array(data = NA, dim =c(m,kfactorial))
#t denotes EM iterations
maxiter<-100
t<-0
threshold<-10**(-6)
criterion <- 99
while((criterion>threshold)&&t<maxiter){
t<-t+1
temp<-numeric(kfactorial)
#E-step: compute the permutation probabilities
for(iter in 1:m){

for(k in 1:kfactorial){
perm<-permutations[[k]]
temp[k]<-complete(x,perm[z[iter,]],theta.hat)
}

for(k in 1:(kfactorial-1)){
pp[iter,k]<-1/sum(exp(temp-temp[k]))
}
pp[iter,kfactorial]<- 1 - sum(pp[iter,1:(kfactorial-1)])
if(is.na(max(pp[iter,]))==TRUE){pp[iter,]=rep(0,kfactorial);ind<-order(temp)[kfactorial];pp[iter,ind]<-1}
index<-which(pp[iter,]<0)
if(length(index)>0){pp[iter,index]=rep(0,length(index))}

}

#M-step: update parameters estimate
theta.hat<-array(data = 0, dim =c(K,J))
for(iter in 1:m){
for(k in 1:kfactorial){
for(j in 1:J){
perm<-permutations[[k]]
theta.hat[,j]<-theta.hat[,j] + pp[iter,k]*mcmc.pars[iter,perm,j]
}
}
}

theta.hat<-theta.hat/(m)
current<-theta.hat
criterion<-max(abs(previous - current))
previous<-theta.hat
#print(paste("iter = ",t, ", criterion = ",criterion))
#print(paste(theta.hat))


#this ends t

}
perm <- array(data = NA, dim=c(m,K));
#sample of permutations per iteration
for (i in 1:m){
k<-sample(1:kfactorial,1,prob = pp[i,])
perm[i,]<-order(permutations[[k]])
}

#time<-proc.time() - ptm
#print(paste("SJW algorithm time:", as.numeric(time[3]),"seconds."))

results<-list(perm,t)
names(results)<-c("permutations","iterations")
return(results)

}



##########################################################################################
##########################################################################################
#		PRA
pra<-function(mcmc.pars,pivot){
K<-dim(mcmc.pars)[2]
J<-dim(mcmc.pars)[3]
m<-dim(mcmc.pars)[1]
l1<-dim(pivot)
if((l1[1]!=K)||(l1[2]!=J)){
stop("Pivot and MCMC samples should have the same length")
}
perms <- array(data = NA, dim = c(m,K))
#library(combinat)
#ptm <- proc.time()
permutations<-permn(1:K)
kfactorial<-gamma(K+1)
for(iter in 1:m){
	t<-1
	perm<-permutations[t][[1]]
	inner.product<-sum(pivot*mcmc.pars[iter,perm,])
	max.val<-inner.product
	perms[iter,]<-perm
	for(t in 2:kfactorial){
		perm<-permutations[t][[1]]
		inner.product<-sum(pivot*mcmc.pars[iter,perm,])
		if(inner.product>max.val){perms[iter,]<-perm;max.val<-inner.product}
	}
}
#	time<-proc.time() - ptm
#	print(paste("PRA time:", as.numeric(time[3]),"seconds."))

	results<-list(perms)
	names(results)<-c("permutations")
	return(results)
}


############################################################################
############################################################################
############################################################################
###########################################################################





############################################################################


permute.mcmc<-function(mcmc,permutations){
m<-dim(permutations)[1]
K<-dim(permutations)[2]
J<-dim(mcmc)[3]
mcmc.permuted<-mcmc
	for(iter in 1:m){
		for(j in 1:J){
			mcmc.permuted[iter,,j]<-mcmc[iter,permutations[iter,],j]
		}	
	}
	results<-list(mcmc.permuted)
	names(results)<-c("output")
	return(results)
}

################################################################################################
################################################################################################
################################################################################################



label.switching<-function(method, zpivot, z, K, prapivot, p, complete, mcmc, sjwinit, data){
    L <- length(method)
    for (l in 1:L) {
        if ((method[l] %in% c("ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", 
            "STEPHENS", "SJW", "PRA")) == FALSE) {
            stop(paste("method:", method[l], "is not recognised"))
        }
    }
    if (missing(z) == FALSE) {
        m <- dim(z)[1]
    }
    else {
        if (missing(mcmc) == FALSE) {
            m <- dim(mcmc)[1]
        }
        else {
            if (missing(p) == FALSE) {
                m <- dim(p)[1]
            }
            else {
                stop(paste("At least one of z, mcmc, or p should be provided"))
            }
        }
    }
    fr <- 1
    dimname <- L
    if ((is.array(zpivot) == TRUE) && (("ECR" %in% method) == 
        TRUE)) {
        fr <- dim(zpivot)[1]
        dimname <- L + fr - 1
    }
    if ("ECR" %in% method) {
        ind <- which(method == "ECR")
        if (length(ind) > 1) {
            stop(paste("ECR appearing more than 1 times"))
        }
        nam <- numeric(dimname)
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind - 1)) {
                nam[i] <- method[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- paste("ECR", i, sep = "-")
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- method[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "ECR"
        }
    }
    else {
        nam <- method
    }
    permutations <- vector("list", length = dimname)
    names(permutations) <- nam
    timings <- numeric(dimname)
    names(timings) <- nam
    f <- 0
    t <- 1
    cat(paste("Method                        Time (sec)\n"))
    for (l in 1:L) {
        cat(paste("   ", method[l]))
        if (method[l] == "ECR") {
            if (is.array(zpivot) == TRUE) {
                while (f < dim(zpivot)[1]) {
                  f <- f + 1
                  if (f == 1) {
                    cat(paste(" (using pivot", f, "of", dim(zpivot)[1]), 
                      ")", sep = "")
                  }
                  else {
                    cat(paste("        (using pivot", f, "of", 
                      dim(zpivot)[1]), ")", sep = "")
                  }
                  if (dim(zpivot)[2] != dim(z)[2]) {
                    stop(paste("length(zpivot) and number of columns of z are not equal"))
                  }
                  if (K < max(z)) {
                    stop(paste("K should be at least equal to", 
                      max(z)))
                  }
                  tpm <- proc.time()
                  permutations[[t]] <- ecr(zpivot[f, ], z, K)$permutations
                  time <- proc.time() - tpm
                  time <- as.numeric(time[3])
                  cat(paste("   ", time, "\n"))
                  timings[t] <- time
                  t <- t + 1
                }
            }
            else {
                if (missing(zpivot)) {
                  stop(paste("zpivot is missing"))
                }
                else {
                  if (length(zpivot) != dim(z)[2]) {
                    stop(paste("length(zpivot) and number of columns of z are not equal"))
                  }
                  if (K < max(z)) {
                    stop(paste("K should be at least equal to", 
                      max(z)))
                  }
                  tpm <- proc.time()
                  permutations[[t]] <- ecr(zpivot, z, K)$permutations
                  time <- proc.time() - tpm
                  time <- as.numeric(time[3])
                  cat(paste("                        ", time, 
                    "\n"))
                  timings[t] <- time
                  t <- t + 1
                }
            }
        }
        if (method[l] == "ECR-ITERATIVE-1") {
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- ecr.iterative.1(z, K)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- as.numeric(time[3])
            cat(paste("            ", time, "  (converged at", 
                hold$iterations, "iterations)\n"))
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "ECR-ITERATIVE-2") {
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(p)) {
                stop(paste("p is missing"))
            }
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- ecr.iterative.2(z, K, p)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- as.numeric(time[3])
            cat(paste("            ", time, "  (converged at", 
                hold$iterations, "iterations)\n"))
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "PRA") {
            tpm <- proc.time()
            permutations[[t]] <- pra(mcmc, prapivot)$permutations
            time <- proc.time() - tpm
            time <- as.numeric(time[3])
            cat(paste("                        ", time, "\n"))
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "STEPHENS") {
            if (missing(p)) {
                stop(paste("p is missing"))
            }
            tpm <- proc.time()
            hold <- stephens(p)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- as.numeric(time[3])
            cat(paste("                   ", time, "  (converged at", 
                hold$iterations, "iterations)\n"))
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "SJW") {
            if (missing(mcmc)) {
                stop(paste("mcmc is missing"))
            }
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(complete)) {
                stop(paste("Complete log-likelihood function is missing"))
            }
            if (missing(data)) {
                stop(paste("Data is missing"))
            }
            if (missing(sjwinit)) {
                sjwinit = 0
            }
            tpm <- proc.time()
            hold <- sjw(mcmc, z, complete, x = data, sjwinit)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- as.numeric(time[3])
            cat(paste("                        ", time, "  (converged at", 
                hold$iterations, "iterations)\n"))
            timings[t] <- time
            t <- t + 1
        }
    }
    cat(paste("\n\n"))

    if (missing(z) == FALSE){
        n <- dim(z)[2]
	best.clusterings <- array(data = NA, dim = c(dimname, n))
        st <- 1:K
        m <- dim(z)[1]
        perm <- numeric(K)
        z1 <- z2 <- numeric(n)
        temp <- z
        for (i in 1:m) {
            temp[i, ] <- permutations[[1]][i, ][z[i, ]]
        }
        z1 <- apply(temp, 2, function(y) order(table(y), decreasing = T)[1])
        zpivot <- z1
        best.clusterings[1, ] <- zpivot
    }



    if ((dimname > 1) && (missing(z) == FALSE)) {
        cat(paste("Reordering outputs to maximize their correlation with method:", 
            method[1], "..."))
        for (l in 2:dimname) {
            for (i in 1:m) {
                temp[i, ] <- permutations[[l]][i, ][z[i, ]]
            }
            z2 <- apply(temp, 2, function(y) order(table(y), 
                decreasing = T)[1])
            st <- 1:K
            k <- K
            cost.matrix <- matrix(numeric(k * k), nrow = k, ncol = k)
            s <- 1:n
            alloc <- z2
            for (i in 1:K) {
                so <- zpivot[s[alloc == i]]
                l2 <- length(alloc[alloc == i])
                for (j in 1:k) cost.matrix[i, j] <- l2 - length(so[so == 
                  j])
            }
            matr <- lp.assign(cost.matrix)$solution
            for (i in 1:k) {
                perm[i] <- st[matr[, i] > 0]
            }
            for (i in 1:m) {
                permutations[[l]][i, ] <- permutations[[l]][i, 
                  ][perm]
		temp[i, ] <- permutations[[l]][i, ][z[i, ]]
            }
            best.clusterings[l, ] <- apply(temp, 2, function(y) order(table(y), 
                decreasing = T)[1])

        }
        rownames(best.clusterings) <- nam
        similarity.matrix <- array(data = 0, dim = c(dimname, 
            dimname))
        colnames(similarity.matrix) <- nam
        rownames(similarity.matrix) <- nam
        for (i in 1:dimname) {
            for (j in 1:i) {
                similarity.matrix[i, j] <- length(which(t(best.clusterings)[, 
                  i] == t(best.clusterings)[, j]))
            }
        }
        similarity.matrix = (similarity.matrix + t(similarity.matrix))/n
        diag(similarity.matrix) <- rep(1, dimname)
    }
    cat(paste(" done!\n"))
    cat(paste("Retrieve the", dimname, "permutation arrays by typing:\n"))
    for (i in 1:dimname) {
        cat(paste("   [...]$permutations$\"", nam[i], "\"", sep = "", 
            "\n"))
    }
    cat(paste("Retrieve the", dimname, "best clusterings by typing: [...]$clusters\n"))
    cat(paste("Retrieve the", dimname, "CPU times by typing: [...]$timings\n"))
    cat(paste("Retrieve the", dimname, "X", dimname, "similarity matrix: [...]$similarity\n"))
    results <- list(permutations, best.clusterings, timings, 
        similarity.matrix)
    names(results) <- c("permutations", "clusters", "timings", 
        "similarity")
    return(results)
}


compare.clust<-function(pivot.clust,perms,z,K){
	if(is.list(perms)==FALSE){rrr<-list(perms);names(rrr)<-deparse(substitute(perms));perms<-rrr}
	if(length(pivot.clust)!=dim(z)[2]){stop("pivot.clust should have the same length with dim(z)[2]")}
	outdim<-length(perms) + 1
	n<-length(pivot.clust)
        best.clusterings <- array(data = NA, dim = c(outdim - 1,n))
        st <- 1:K
        m <- dim(z)[1]
	temp <- z
        perm <- numeric(K)
        z1 <- z2 <- numeric(n)
        zpivot <- pivot.clust
	permutations<-vector("list",length=outdim - 1)
	for(l in 1:(outdim - 1)){permutations[[l]]<-perms[[l]]}
        for (l in 1:(outdim - 1)) {
            for (i in 1:m) {
                temp[i, ] <- permutations[[l]][i, ][z[i, ]]
            }
            z2 <- apply(temp, 2, function(y) order(table(y), 
                decreasing = T)[1])
            st <- 1:K
            k <- K
            cost.matrix <- matrix(numeric(k * k), nrow = k, ncol = k)
            s <- 1:n
            alloc <- z2
            for (i in 1:K) {
                so <- zpivot[s[alloc == i]]
                l2 <- length(alloc[alloc == i])
                for (j in 1:k) cost.matrix[i, j] <- l2 - length(so[so == j])
            }
            matr <- lp.assign(cost.matrix)$solution
            for (i in 1:k) {perm[i] <- st[matr[, i] > 0]}
            for (i in 1:m) {
                permutations[[l]][i, ] <- permutations[[l]][i, ][perm]
		temp[i, ] <- permutations[[l]][i, ][z[i, ]]
            }
	    best.clusterings[l, ] <- apply(temp, 2, function(y) order(table(y), decreasing = T)[1])

        }
        similarity.matrix <- array(data = 0, dim = c(outdim,outdim))
        for (i in 1:(outdim - 1)) {
            for (j in 1:i) {
                similarity.matrix[i, j] <- length(which(t(best.clusterings)[,i] == t(best.clusterings)[, j]))
            }
        }
	i<-outdim
        for (j in 1:(outdim - 1)) {
              similarity.matrix[i, j] <- length(which(zpivot == t(best.clusterings)[, j]))}
	
        similarity.matrix = (similarity.matrix + t(similarity.matrix))/n
        diag(similarity.matrix) <- rep(1, outdim)
	rownames(similarity.matrix)<-c(names(perms),deparse(substitute(pivot.clust)))
	colnames(similarity.matrix)<-c(names(perms),deparse(substitute(pivot.clust)))
	results<-list(similarity.matrix,best.clusterings)
	names(results)<-c("similarity","clusters")
	return(results)
   
}




