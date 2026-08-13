NimModel <- nimbleCode({
  psi ~ dbeta(1,1) #uniform, but get conjugate sampler
  for(g in 1:n.primary){
    alpha[g] <- 1
  }
  pi[1:n.primary] ~ ddirch(alpha[1:n.primary]) #entry probabilities
  #conditional entry probabilities. eta[n.primary] is 1 by construction
  eta[1] <- pi[1]
  for(g in 2:n.primary){
    eta[g] <- pi[g]/(1 - sum(pi[1:(g-1)]))
  }
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10)
  phi.cov.mu ~ dunif(-10,10)
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf)
  for(g in 1:n.primary){
    p[g] ~ dunif(0,1)
  }
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu,sd = phi.cov.sd)
    logit(phi[i]) <- beta0.phi + phi.cov[i]*beta1.phi
    y[i,1:n.primary] ~ dJS(eta=eta[1:n.primary],p=p[1:n.primary],
                           K=K[1:n.primary],phi=phi[i],psi=psi)
  }
  #placeholder priors never used because NBSampler simulates these directly from their
  #full conditionals. do not assign them default samplers
  for(g in 1:n.primary){
    N[g] ~ dpois(1)
    B[g] ~ dpois(1)
  }
  N.super ~ dpois(1)
})