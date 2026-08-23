NimModel <- nimbleCode({
  psi ~ dbeta(1,1) #uniform, but get conjugate sampler
  for(g in 1:n.primary){
    alpha[g] <- 1
  }
  for(g in 1:(n.primary-1)){
    eta[g] ~ dbeta(alpha[g], sum(alpha[(g+1):n.primary]))
  }
  eta[n.primary] <- 1
  #derive pi
  pi[1] <- eta[1]
  for(g in 2:(n.primary-1)){
    pi[g] <- eta[g] * prod(1 - eta[1:(g-1)])
  }
  pi[n.primary] <- prod(1 - eta[1:(n.primary-1)])
  
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10) #individual covariate effect
  phi.cov.mu ~ dunif(-10, 10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  
  #population dynamics
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu,sd=phi.cov.sd) #individual survival covs
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi
    z.super[i] ~ dbern(psi)
    z[i,1] ~ dbern(eta[1])
    a[i,1] <- 1 - z[i,1]
    for(g in 2:n.primary){
      mu[i,g] <- phi[i]*z[i,g-1] + eta[g]*a[i,g-1]
      z[i,g] ~ dbern(mu[i,g])
      a[i,g] <- a[i,g-1]*(1 - z[i,g])
    }
  }
  #observation model
  for(g in 1:n.primary){
    p[g] ~ dbeta(1,1)
    for(i in 1:M){
      y[i,g] ~ dbinom(p=p[g],size=K[g]*z[i,g]*z.super[i])
    }
  }
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g]*z.super[1:M])
    Acum[g] <- sum((1-a[1:M,g])*z.super[1:M])
  }
  B[1] <- Acum[1]
  for(g in 2:n.primary){
    B[g] <- Acum[g] - Acum[g-1]
  }
  N.super <- sum(z.super[1:M])
})
