NimModel <- nimbleCode({
  ##Abundance##
  lambda.y1 ~ dunif(0,1000) #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.primary){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #abundance by primary occasion
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.primary-1)]) #size of superpopulation
  
  gamma ~ dunif(0,2) #fixed recruitment parameter
  for(g in 1:(n.primary-1)){
    # gamma[g] ~ dunif(0,2) #recruitment priors by primary occasion
    # ER[g] <- N[g]*gamma[g]*tau[g] #expected recruits, variable gamma
    ER[g] <- N[g]*gamma*tau[g] #expected recruits, if gamma fixed
    N.recruit[g] ~ dpois(ER[g]) #realized recruits
  }
  
  #Individual covariates
  phi.cov.mu ~ dunif(-10, 10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  sigma.move ~ dunif(0,2)
  for(g in 1:(n.primary-1)){#time scaled sigma move
    sigma.move.int[g] <- sigma.move*sqrt(tau.move[g])
  }
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu,sd=phi.cov.sd)
    #1st primary occasion ACs
    # s[i,1,1] ~ dunif(xlim[1],xlim[2])
    # s[i,1,2] ~ dunif(ylim[1],ylim[2])
    #same as above, but vectorized and gated by z.super
    #all s set to 0 if not in population, z.super[i]=0
    s[i,1,1:2] ~ dunif2D(xlim=xlim[1:2],ylim=ylim[1:2],z.super=z.super[i])
    for(g in 2:n.primary){
      # s[i,g,1] ~ T(dnorm(s[i,g-1,1],sd=sigma.move.int[g-1]),xlim[1],xlim[2])
      # s[i,g,2] ~ T(dnorm(s[i,g-1,2],sd=sigma.move.int[g-1]),ylim[1],ylim[2])
      #same as above, but vectorized and gated by z.super
      s[i,g,1:2] ~ dTruncNorm(s.prev=s[i,g-1,1:2],sigma.move=sigma.move.int[g-1], 
                                       xlim=xlim[1:2],ylim=ylim[1:2],
                                       z.super=z.super[i])
    }
  }
  
  #Survival (phi must have M x n.primary - 1 dimension for custom updates to work)
  #without individual or primary occasion effects, use for loop to plug into phi[i,g]
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0, sd=10) #individual covariate effect on survival
  for(i in 1:M){
    #unit-time individual survival probability
    logit(phi.unit[i]) <- beta0.phi + beta1.phi*phi.cov[i]
    #survival over each actual interval
    for(g in 1:(n.primary-1)){#plugging same individual phi's into each primary occasion for custom update
      phi[i,g] <- phi.unit[i]^tau[g]
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.primary] ~ dSurvival(phi=phi[i,1:(n.primary-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
  }
  
  ##Detection##
  sigma ~ dunif(0,10) #fixing sigma across primary occasions
  for(g in 1:n.primary){
    p0[g] ~ dunif(0,1) #p0 varies by primary occasion
    for(i in 1:M){ #only compute pd and y when z.super[i]=1&z[i,g]=1
      pd[i,g,1:J[g]] <- GetDetectionProb(s = s[i,g,1:2], X = X[g,1:J[g],1:2],
                                         J=J[g], sigma=sigma, p0=p0[g], z=z[i,g], z.super=z.super[i])
      y[i,g,1:J[g]] ~ dBinomialVector(pd = pd[i,g,1:J[g]], K = K1D[g,1:J[g]],
                                      z = z[i,g], z.super=z.super[i]) #vectorized obs mod
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update
#4) s update for z.super=1 only
