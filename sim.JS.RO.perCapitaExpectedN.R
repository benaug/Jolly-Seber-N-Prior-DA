sim.JS.RO.perCapitaExpectedN <- function(lambda1=NA,gamma=NA,phi=NA,
                       p=NA,n.primary=NA,K=NA,M=NA,tau=NA){
  if(length(gamma)!=(n.primary-1))stop("gamma must be of length n.primary-1")
  if(length(phi)!=(n.primary-1))stop("phi must be of length n.primary-1")
  if(length(tau)!=(n.primary-1))stop("tau must be of length n.primary-1")
  if(lambda1 <= 0)stop("lambda1 must be > 0")
  if(any(gamma <= 0))stop("gamma must be > 0")
  if(any(phi <= 0 | phi >= 1))stop("phi must be between 0 and 1")

  phi.int <- phi^tau
  lambda <- EN <- rep(0,n.primary)
  lambda[1] <- EN[1] <- lambda1
  for(g in 2:n.primary){
    lambda[g] <- EN[g-1]*gamma[g-1]*tau[g-1] #linear time scaling
    EN[g] <- EN[g-1]*phi.int[g-1] + lambda[g]
  }

  lambda.super <- sum(lambda)
  if(lambda.super >= M)stop("Expected N.super must be less than M")

  psi <- lambda[1]/M
  gamma.RO <- rep(0,n.primary-1)
  cum.lambda <- lambda[1]
  for(g in 1:(n.primary-1)){
    gamma.RO[g] <- lambda[g+1]/(M-cum.lambda)
    cum.lambda <- cum.lambda+lambda[g+1]
  }

  #####Population Dynamics############
  z <- a <- matrix(0,M,n.primary)
  z[,1] <- rbinom(M,1,psi)
  a[,1] <- z[,1]

  for(g in 2:n.primary){
    Ez <- z[,g-1]*phi.int[g-1] + (1-a[,g-1])*gamma.RO[g-1]
    z[,g] <- rbinom(M,1,Ez)
    a[,g] <- pmax(a[,g-1],z[,g])
  }

  N <- colSums(z)
  Acum <- colSums(a)
  B <- diff(Acum)
  N.super <- Acum[n.primary]

  #detection
  y <- z*0
  for(g in 1:n.primary){
    y[,g] <- rbinom(M,K[g],p[g]*z[,g])
  }

  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,,drop=FALSE]
  z <- z[keep.idx,,drop=FALSE]
  return(list(y=y,N=N,B=B,z=z,N.super=N.super,EN=EN,n.primary=n.primary,K=K,tau=tau,
              lambda=lambda,lambda.super=lambda.super,psi=psi,gamma.RO=gamma.RO,
              gamma=gamma,phi=phi,phi.int=phi.int))
}
