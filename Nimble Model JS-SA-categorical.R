NimModel <- nimbleCode({
  psi ~ dbeta(1,1) #uniform, but get conjugate sampler
  for(g in 1:n.primary){
    alpha[g] <- 1
  }
  pi[1:n.primary] ~ ddirch(alpha[1:n.primary]) #entry probabilities
  
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10) #individual covariate effect
  phi.cov.mu ~ dunif(-10, 10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  
  #population dynamics
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu,sd=phi.cov.sd) #individual survival covs
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi
    z.super[i] ~ dbern(psi)
    e[i] ~ dcat(pi[1:n.primary]) #entry occasion
    z[i,1] <- equals(e[i],1)
    recruit[i,1] <- z.super[i]*z[i,1]
    for(g in 2:n.primary){
      surv.p[i,g] <- phi[i]*z[i,g-1] #keep separate from line below for custom updates
      surv[i,g] ~ dbern(surv.p[i,g]) #survival outcome
      z[i,g] <- max(surv[i,g],equals(e[i],g)) #alive if survived or entering now
      recruit[i,g] <- z.super[i]*equals(e[i],g) #used to compute B
    }
  }
  #observation model
  for(g in 1:n.primary){
    p[g] ~ dunif(0,1)
    for(i in 1:M){
      y.p[i,g] <- p[g]*z[i,g]*z.super[i] #keep separate from line below for custom updates
      y[i,g] ~ dbinom(p=y.p[i,g],size=K[g])
    }
  }
  #derived variables
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g]*z.super[1:M])
    B[g] <- sum(recruit[1:M,g])
  }
  N.super <- sum(z.super[1:M])
})