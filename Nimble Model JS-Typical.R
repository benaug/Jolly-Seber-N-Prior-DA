NimModel <- nimbleCode({
  #First year population size
  psi ~ dunif(0,1)
  beta0.phi ~ dlogis(0, 1)
  beta1.phi ~ dnorm(0,sd=10) #individual covariate effect
  #population dynamics
  for(g in 1:n.year){
    N[g] <- sum(z[1:M,g])
    gamma[g] ~ dunif(0,10)
    ER[g] <- N[g]*gamma[g] #expected recruits
    A[g] <- max(M - sum(a[1:M,g]),0.01) #available recruits
    gamma.prime[g] <- min(ER[g]/A[g],0.999) #individual recruitment prob
  }
  for(i in 1:M){
    phi.cov[i] ~ dnorm(0,1) #individual survival covs
    z[i,1] ~ dbern(psi) #typical way to do first year population size
    a[i,1] <- z[i,1]
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi #individual survival prob
    for(g in 2:n.year){
      Ez[i,g-1] <- z[i,g-1]*phi[i] + (1 - a[i,g-1])*gamma.prime[g-1]
      z[i,g] ~ dbern(Ez[i,g-1])
      a[i,g] <- max(z[i,1:g]) #not available to recruit if previously alive
    }
  }
  #observation model
  for(g in 1:n.year){
    p[g] ~ dunif(0,1)
    for(i in 1:M){
      y[i,g] ~ dbinom(p=p[g]*z[i,g],size=K[g])
    }
  }
})