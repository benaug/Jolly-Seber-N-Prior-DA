NimModel <- nimbleCode({
  ##Abundance##
  lambda.y1.M ~ dunif(0,1000) #Expected starting population size, males
  lambda.y1.F ~ dunif(0,1000) #Expected starting population size, females
  N.M[1] ~ dpois(lambda.y1.M) #Realized starting population size
  N.F[1] ~ dpois(lambda.y1.F) #Realized starting population size
  N[1] <- N.M[1] + N.F[1]
  for(g in 2:n.primary){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #total abundance by primary occasion
    N.M[g] <- N.survive.M[g-1] + N.recruit.M[g-1] #male abundance by primary occasion
    N.F[g] <- N.survive.F[g-1] + N.recruit.F[g-1] #female abundance by primary occasion
    #sex-specific N.recruit and N.survive information also contained in z, z.start, z.stop, and sex
    #sex-specific N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.primary-1)]) #size of superpopulation
  
  #Sex-specific Recruitment
  #fixed male and female recruitment priors
  gamma.sex[1] ~ dunif(0,2)
  gamma.sex[2] ~ dunif(0,2)
  for(g in 1:(n.primary-1)){
    #fixed gamma.sex
    ER.M[g] <- N[g]*gamma.sex[1] #male expected recruits per total N
    ER.F[g] <- N[g]*gamma.sex[2] #female expected recruits per total N
    #gamma.sex by primary occasion
    # gamma.sex[1,g] ~ dunif(0,2) #male recruitment priors by primary occasion
    # gamma.sex[2,g] ~ dunif(0,2) #female recruitment priors by primary occasion
    # ER.M[g] <- N[g]*gamma.sex[1,g] #male expected recruits per total N
    # ER.F[g] <- N[g]*gamma.sex[2,g] #female expected recruits per total N
    N.recruit.M[g] ~ dpois(ER.M[g]) #male realized recruits
    N.recruit.F[g] ~ dpois(ER.F[g]) #female realized recruits
  }
  
  #Individual covariates. individual sex info derived from recruitment model
  sigma.move.sex[1] ~ dunif(0,5) #male movement sigma
  sigma.move.sex[2] ~ dunif(0,5) #female movement sigma
  for(i in 1:M){
    #1st primary occasion ACs
    # s[i,1,1] ~ dunif(xlim[1],xlim[2])
    # s[i,1,2] ~ dunif(ylim[1],ylim[2])
    #same as above, but vectorized and gated by z.super
    #all s set to 0 if not in population, z.super[i]=0
    s[i,1,1:2] ~ dunif2D(xlim=xlim[1:2],ylim=ylim[1:2],z.super=z.super[i])
    for(g in 2:n.primary){
      #If you change movement distribution, need to modify custom s updates (3rd one for z.super=0)
      # s[i,g,1] ~ T(dnorm(s[i,g-1,1],sd=sigma.move.sex[sex[i]+1]),xlim[1],xlim[2])
      # s[i,g,2] ~ T(dnorm(s[i,g-1,2],sd=sigma.move.sex[sex[i]+1]),ylim[1],ylim[2])
      #same as above, but vectorized and gated by z.super
      s[i,g,1:2] ~ dTruncNorm(s.prev=s[i,g-1,1:2],sigma.move=sigma.move.sex[sex[i]+1], 
                              xlim=xlim[1:2],ylim=ylim[1:2],
                              z.super=z.super[i])
    }
  }

  #Survival (phi must have M x n.primary - 1 dimension for custom updates to work)
  #fixed sex-specific survival
  phi.sex[1] ~ dunif(0,1)
  phi.sex[2] ~ dunif(0,1)
  #sex-specific survival by primary occasion
  # for(g in 1:(n.primary-1)){
  #   phi.sex[1,g] ~ dunif(0,1)
  #   phi.sex[2,g] ~ dunif(0,1)
  # }
  for(i in 1:M){
    for(g in 1:(n.primary-1)){ #plugging same individual phi's into each primary occasion (phi: M x n.primary-1 expected by custom update)
      phi[i,g] <- phi.sex[sex[i]+1] #if fixed
      # phi[i,g] <- phi.sex[sex[i]+1,g] #if occasion-specific
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.primary] ~ dSurvival(phi=phi[i,1:(n.primary-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
  }

  ##Detection##
  #fixing sex-specific sigma across primary occasions
  sigma.sex[1] ~ dunif(0,10) #male sigma
  sigma.sex[2] ~ dunif(0,10) #female sigma
  for(g in 1:n.primary){
    #sex-specific p0 varies by primary occasions
    p0.sex[1,g] ~ dunif(0,1) #male p0
    p0.sex[2,g] ~ dunif(0,1) #female p0
    for(i in 1:M){ #only compute pd and y likelihood when z.super[i]=1&z[i,g]=1
      pd[i,g,1:J[g]] <- GetDetectionProb(s = s[i,g,1:2], X = X[g,1:J[g],1:2],
                                         J=J[g], sigma=sigma.sex[sex[i]+1], p0=p0.sex[sex[i]+1,g], z=z[i,g], z.super=z.super[i])
      y[i,g,1:J[g]] ~ dBinomialVector(pd = pd[i,g,1:J[g]], K = K1D[g,1:J[g]],
                                      z = z[i,g], z.super=z.super[i]) #vectorized obs mod
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update
#4) sex update: update sex-specific N structures with sex[i] update
#5) s update for z.super=1 only
