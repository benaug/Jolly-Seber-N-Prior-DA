sim.JS.typical <- function(psi=NA,gamma=NA,beta0.phi=NA,beta1.phi=NA,
                   p=NA,n.year=NA,K=NA,M=NA){
  #####Population Dynamics############
  N1 <- rbinom(1,M,psi)
  z <- a <- matrix(NA,M,n.year)
  z[1:N1,1] <- 1
  z[(N1+1):M] <- 0
  a[,1] <-1-z[,1] # Available to be recruited?
  gamma.prime <- rep(NA,n.year-1)
  phi.cov <- rnorm(M,0,1) #simulate ind survival covariate
  phi <- matrix(NA,M,n.year-1) #M x n.year-1 survival probs
  
  for(g in 2:n.year) {
    ER <- sum(z[,g-1])*gamma[g-1] # Expected number of recruits
    A <- sum(a[,g-1]) # nAvailable to be recruited
    if(ER>A){
      stop("M is too low. There aren't any individuals left to be recruited")
    }
    gamma.prime[g-1] <- ER/A # individual-level recruitment *probability*
    phi[,g-1] <- plogis(beta0.phi+phi.cov*beta1.phi) #individual by year survival probability
    Ez <- z[,g-1]*phi[,g-1] + a[,g-1]*gamma.prime[g-1]
    z[,g] <- rbinom(M, 1, Ez)
    a[,g] <- apply(z[,1:g]<1, 1, all)
  }
  N <- colSums(z)
  
  #detection
  y <- z*0
  for(g in 1:n.year){
      y[,g] <- rbinom(nrow(z),K[g],p[g]*z[,g])
  }
  
  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,]
  phi.cov <- phi.cov[keep.idx]
  z <- z[keep.idx,]
  return(list(y=y,phi.cov=phi.cov,N=N,z=z))
}