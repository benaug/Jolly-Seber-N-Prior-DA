NimModel <- nimbleCode({
  #Royle and Dorazio 2008 priors
  #First year expected population size
  psi ~ dbeta(1,1) #uniform, but get conjugate sampler
  #entry priors
  for(g in 1:(n.primary-1)){
    gamma[g] ~ dunif(0,1)
  }
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10) #individual covariate effect
  #population dynamics
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g])
    Acum[g] <- sum(a[1:M,g])
  }
  phi.cov.mu ~ dunif(-10,10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  for(i in 1:M){
    phi.cov[i] ~  dnorm(phi.cov.mu, sd = phi.cov.sd) #individual survival covs
    z[i,1] ~ dbern(psi)
    a[i,1] <- z[i,1]
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi #individual survival prob
    for(g in 2:n.primary){
      Ez[i,g-1] <- z[i,g-1]*phi[i] + (1-a[i,g-1])*gamma[g-1] #variable gamma
      z[i,g] ~ dbern(Ez[i,g-1])
      a[i,g] <- max(a[i,g-1],z[i,g]) #not available to recruit if previously alive
    }
  }
  #compute realized recruits
  for(g in 2:n.primary){
    B[g-1] <- Acum[g] - Acum[g-1]
  }
  #observation model
  for(g in 1:n.primary){
    p[g] ~ dbeta(1,1)
    for(i in 1:M){
      y[i,g] ~ dbinom(p=p[g],size=K[g]*z[i,g])
    }
  }
})
