sim.JS.SA <- function(lambda1=NA,phi=NA,gamma=NA,
                      p=NA,n.primary=NA,K=NA,M=NA,tau=NA){
  if(length(phi)!=(n.primary-1))stop("phi must be of length n.primary-1")
  if(length(gamma)!=(n.primary-1))stop("gamma must be of length n.primary-1")
  if(length(tau)!=(n.primary-1))stop("tau must be of length n.primary-1")
  if(any(phi <= 0 | phi >= 1))stop("phi must be between 0 and 1")
  if(any(gamma <= 0))stop("gamma must be > 0")
  if(lambda1 <= 0)stop("lambda1 must be > 0")

  #absolute expected abundance and expected entries
  lambda <- EN <- rep(NA,n.primary)
  lambda[1] <- EN[1] <- lambda1
  phi.int <- phi^tau
  for(g in 1:(n.primary-1)){
    #linear scaling of per-capita recruitment with interval length
    lambda[g+1] <- EN[g]*gamma[g]*tau[g]
    EN[g+1] <- EN[g]*phi.int[g] + lambda[g+1]
  }

  lambda.super <- sum(lambda)
  if(lambda.super >= M)stop("Expected N.super must be less than M")

  psi <- lambda.super/M
  pi <- lambda/lambda.super

  #####Population Dynamics############
  z.super <- rbinom(M,1,psi) #superpopulation membership
  N.super <- sum(z.super)
  if(N.super==M)stop("Maxed out M")

  e <- rep(0,M) #entry occasion, 0 = never a member
  e[z.super==1] <- sample(1:n.primary,N.super,replace=TRUE,prob=pi)

  z <- matrix(0,M,n.primary)
  z[cbind(which(z.super==1),e[z.super==1])] <- 1 #alive at entry
  for(g in 2:n.primary){
    surv <- rbinom(M,1,phi.int[g-1])*z[,g-1] #survivors from g-1
    z[,g] <- pmax(z[,g],surv)
  }

  N <- colSums(z)
  B <- tabulate(e[z.super==1],nbins=n.primary)

  #detection
  y <- z*0
  for(g in 1:n.primary){
    y[,g] <- rbinom(M,K[g],p[g]*z[,g])
  }

  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,,drop=FALSE]
  z <- z[keep.idx,,drop=FALSE]

  return(list(y=y,N=N,B=B,z=z,N.super=N.super,n.primary=n.primary,K=K,
              EN=EN,lambda=lambda,lambda.super=lambda.super,
              phi.int=phi.int,psi=psi,pi=pi))
}
