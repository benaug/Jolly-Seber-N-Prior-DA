sim.JS.typical <- function(psi=NA,gamma=NA,beta0.phi=NA,beta1.phi=NA,
                   p=NA,n.year=NA,K=NA,M=NA){
  #####Population Dynamics############
  N1 <- rbinom(1,M,psi)
  if(N1==M)stop("Maxed out M in year 1")
  z <- a <- matrix(NA,M,n.year)
  z[1:N1,1] <- 1
  z[(N1+1):M,1] <- 0
  a[,1] <- z[,1] # Available to be recruited?
  gamma.prime <- rep(NA,n.year-1)
  phi.cov <- rnorm(M,0,1) #simulate ind survival covariate
  phi <- plogis(beta0.phi + phi.cov * beta1.phi)  # individual survival probs, constant across time
  
  for(g in 2:n.year){
    ER <- sum(z[,g-1])*gamma[g-1]
    A <- M - sum(a[,g-1])
    if(A <= 0)stop("Ran out of recruits, raise M")
    gamma.prime[g-1] <- ER/A
    if(gamma.prime[g-1]>1)stop("Recruitment probability >1 in one year, raise M")
    Ez <- z[,g-1]*phi + (1 - a[,g-1])*gamma.prime[g-1]
    z[,g] <- rbinom(M, 1, Ez)
    a[,g] <- apply(z[,1:g], 1, max)  # a=1 if ever been alive
  }
  N <- colSums(z)
  
  if(any(gamma.prime>0.1))warning("At least one recruitment probability >0.1, Poisson approx. likely poor, recommend raising M.")
  
  #detection
  y <- z*0
  for(g in 1:n.year){
    y[,g] <- rbinom(M,K[g],p[g]*z[,g])
  }
  
  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,]
  phi.cov <- phi.cov[keep.idx]
  z <- z[keep.idx,]
  return(list(y=y,phi.cov=phi.cov,N=N,z=z,N.super=nrow(z),gamma.prime=gamma.prime))
}