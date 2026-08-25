NimModel <- nimbleCode({
  #First year population size
  # psi ~ dunif(0,1)
  psi ~ dbeta(1,1) #same as above, but get conjugate sampler
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10) #individual covariate effect
  #population dynamics
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g])
    Acum[g] <- sum(a[1:M,g])
  }
  gamma.fixed ~ dunif(0,2) #fixed gamma, must keep individual gamma[g] nodes for custom z updates to work
  for(g in 1:(n.primary-1)){
    gamma[g] <- gamma.fixed #fixed gamma
    # gamma[g] ~ dunif(0,2) #gamma varies by primary occasion
    ER[g] <- N[g]*gamma[g]*tau[g] #expected recruits
    A.raw[g] <- M - Acum[g] #available recruits
    A[g] <- max(A.raw[g],0.01) #trick to prevent model from crashing, but can bias estimates if it happens
    gamma.prime.raw[g] <- ER[g]/A[g] #individual recruitment prob
    gamma.prime[g] <- min(gamma.prime.raw[g],0.999) #trick to prevent model from crashing, but can bias estimates if it happens
  }
  #compute realized recruits
  for(g in 2:n.primary){
    B[g-1] <- Acum[g] - Acum[g-1]
  }
  phi.cov.mu ~ dunif(-10,10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0,sigma=1,df=7), 0, Inf) #phi individual covariate sd prior
  for(i in 1:M){
    phi.cov[i] ~  dnorm(phi.cov.mu, sd = phi.cov.sd) #individual survival covs
    z[i,1] ~ dbern(psi)
    a[i,1] <- z[i,1]
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi #individual survival prob
    for(g in 2:n.primary){
      phi.int[i,g-1] <- phi[i]^tau[g-1]
      Ez[i,g-1] <- z[i,g-1]*phi.int[i,g-1] + (1-a[i,g-1])*gamma.prime[g-1]
      z[i,g] ~ dbern(Ez[i,g-1])
      a[i,g] <- max(a[i,g-1],z[i,g]) #not available to recruit if previously alive
    }
  }
  #observation model
  for(g in 1:n.primary){
    p[g] ~ dbeta(1,1) #uniform
    for(i in 1:M){
      y[i,g] ~ dbinom(p=p[g],size=K[g]*z[i,g])
    }
  }
})
