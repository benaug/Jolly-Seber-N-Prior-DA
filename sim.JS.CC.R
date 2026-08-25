sim.JS.CC <- function(psi=NA,gamma=NA,tau=NA,beta0.phi=NA,beta1.phi=NA,
                   p=NA,n.primary=NA,K=NA,M=NA){
  if(length(tau)!=(n.primary-1)){
    stop("tau must have length n.primary - 1")
  }
  #####Population Dynamics############
  N1 <- rbinom(1,M,psi)
  if(N1==M)stop("Maxed out M in primary session 1")
  z <- a <- matrix(NA,M,n.primary)
  z[1:N1,1] <- 1
  z[(N1+1):M,1] <- 0
  a[,1] <- z[,1] #1 if individual has entered the population
  gamma.prime <- rep(NA,n.primary-1)
  phi.cov <- rnorm(M,0,1) #simulate ind survival covariate
  phi <- plogis(beta0.phi + phi.cov * beta1.phi)  # individual survival probs, constant across time
  
  for(g in 2:n.primary){
    ER <- sum(z[,g-1])*gamma[g-1]*tau[g-1]
    A <- M - sum(a[,g-1])
    if(A <= 0)stop("Ran out of recruits, raise M")
    gamma.prime[g-1] <- ER/A
    if(gamma.prime[g-1] >= 1){
      stop("Recruitment probability >= 1; raise M")
    }
    if(any(gamma.prime>0.1,na.rm=TRUE))warning("At least one gamma.prime > 0.1; conditional recruitment variance is >10% below that of a Poisson distribution with the same mean. Consider raising M.")
    phi.int <- phi^tau[g-1]
    Ez <- z[,g-1]*phi.int + (1-a[,g-1])*gamma.prime[g-1]
    z[,g] <- rbinom(M, 1, Ez)
    a[,g] <- pmax(a[,g-1], z[,g])
  }
  N <- colSums(z)
  Acum <- colSums(a)
  B <- diff(Acum) #B[1:(n.primary-1)], recruits entering at g+1
  N.super <- Acum[n.primary]
  
  if(any(gamma.prime>0.1))warning("At least one recruitment probability >0.1, Poisson approx. likely poor, recommend raising M.")
  
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
  return(list(y=y,phi.cov=phi.cov,N=N,B=B,z=z,N.super=N.super,gamma.prime=gamma.prime,
              n.primary=n.primary,K=K,tau=tau))
}
