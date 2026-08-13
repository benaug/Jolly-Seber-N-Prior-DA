NimModel <- nimbleCode({
  #expected initial population
  lambda[1] ~ dunif(0,1000)
  #expected recruitment per unit time
  rho ~ dunif(0,1000)
  #expected entries during each interval
  for(g in 2:n.primary){
    lambda[g] <- rho*tau[g-1]
  }
  lambda.super <- sum(lambda[1:n.primary])
  # Required because psi must be <= 1
  lambda.valid ~ dbern(step(M - lambda.super))
  # Derived DA parameters
  psi <- lambda.super/M
  for(g in 1:n.primary){
    pi[g] <- lambda[g]/lambda.super
  }
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10)
  
  phi.cov.mu ~ dunif(-10,10)
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf)
  
  # population dynamics
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu, sd=phi.cov.sd)
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi
    z.super[i] ~ dbern(psi)
    e[i] ~ dcat(pi[1:n.primary])
    z[i,1] <- equals(e[i],1)
    recruit[i,1] <- z.super[i]*z[i,1]
    for(g in 2:n.primary){
      surv.p[i,g] <- phi[i]*z[i,g-1]
      surv[i,g] ~ dbern(surv.p[i,g])
      z[i,g] <- max(surv[i,g],equals(e[i],g))
      recruit[i,g] <- z.super[i]*equals(e[i],g)
    }
  }
  
  #observation model
  for(g in 1:n.primary){
    p[g] ~ dbeta(1,1)
    for(i in 1:M){
      y.size[i,g] <- K[g]*z[i,g]*z.super[i] #keep separate from line below for custom updates
      y[i,g] ~ dbinom(p=p[g],size=y.size[i,g])
    }
  }
  
  # derived variables
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g]*z.super[1:M])
    B[g] <- sum(recruit[1:M,g])
  }
  N.super <- sum(z.super[1:M])
})