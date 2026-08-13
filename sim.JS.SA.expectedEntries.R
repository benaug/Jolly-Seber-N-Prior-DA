sim.JS.SA.expectedEntries <- function(lambda=NA,beta0.phi=NA,beta1.phi=NA,
                      p=NA,n.primary=NA,K=NA,M=NA,tau=NA){
  if(length(lambda)!=n.primary)stop("lambda must be of length n.primary")
  if(any(lambda <= 0))stop("lambda must be > 0")
  
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
  phi.cov <- rnorm(M,0,1) #simulate ind survival covariate
  phi <- plogis(beta0.phi + phi.cov * beta1.phi) #survival over one unit of time
  
  z <- matrix(0,M,n.primary)
  z[cbind(which(z.super==1),e[z.super==1])] <- 1 #alive at entry
  for(g in 2:n.primary){
    phi.int <- phi^tau[g-1]
    surv <- rbinom(M,1,phi.int)*z[,g-1]
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
  y <- y[keep.idx,]
  phi.cov <- phi.cov[keep.idx]
  z <- z[keep.idx,]
  return(list(y=y,phi.cov=phi.cov,N=N,B=B,z=z,N.super=N.super,
              lambda=lambda,lambda.super=lambda.super,psi=psi,pi=pi))
}