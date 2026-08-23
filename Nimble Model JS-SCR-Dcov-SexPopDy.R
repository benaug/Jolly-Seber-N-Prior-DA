NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  # D.intercept <- exp(D.beta0)*cellArea
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) #separate this component so s's do not depend on D.intercept
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  ##Abundance##
  p.sex ~ dunif(0,1) #proportion of expected N in primary occasion 1 that is female
  lambda.y1 <- D.intercept*pi.denom #Expected starting population size
  lambda.y1.M <- lambda.y1*(1-p.sex) #Expected starting population size, males
  lambda.y1.F <- lambda.y1*p.sex #Expected starting population size, females
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
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]]) #categorical likelihood for this cell, equivalent to zero's trick
  }

  ##Survival##
  #phi must have M x (n.primary-1) dimension for custom updates to work
  
  #OPTION 1: one phi shared by sexes and primary occasions
  #phi.sex[1,1] ~ dbeta(1,1)
  #phi.sex[2,1] <- phi.sex[1,1]
  #for(g in 2:(n.primary-1)){
  #  phi.sex[1,g] <- phi.sex[1,1]
  #  phi.sex[2,g] <- phi.sex[1,1]
  #}
  #OPTION 2: occasion-specific phi shared by sexes
  #for(g in 1:(n.primary-1)){
  #  phi.sex[1,g] ~ dbeta(1,1)
  #  phi.sex[2,g] <- phi.sex[1,g]
  #}
  #OPTION 3: fixed sex-specific phi
  phi.sex[1,1] ~ dbeta(1,1) #male survival
  phi.sex[2,1] ~ dbeta(1,1) #female survival
  for(g in 2:(n.primary-1)){
    phi.sex[1,g] <- phi.sex[1,1]
    phi.sex[2,g] <- phi.sex[2,1]
  }
  #OPTION 4: sex-specific phi varies by primary occasion
  # for(g in 1:(n.primary-1)){
  #   phi.sex[1,g] ~ dbeta(1,1) #male survival
  #   phi.sex[2,g] ~ dbeta(1,1) #female survival
  # }
  for(i in 1:M){
    for(g in 1:(n.primary-1)){
      phi[i,g] <- phi.sex[sex[i]+1,g]
    }
    z[i,1:n.primary] ~ dSurvival(phi=phi[i,1:(n.primary-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
  }

  ##Detection##
  #fixing sex-specific sigma across primary occasions
  sigma.sex[1] ~ dunif(0,10) #male sigma
  sigma.sex[2] ~ dunif(0,10) #female sigma
  for(g in 1:n.primary){
    #sex-specific p0 varies by primary occasion
    p0.sex[1,g] ~ dunif(0,1) #male p0
    p0.sex[2,g] ~ dunif(0,1) #female p0
    for(i in 1:M){ #only compute pd and y likelihood when z.super[i]=1&z[i,g]=1
      pd[i,g,1:J[g]] <- GetDetectionProb(s = s[i,1:2], X = X[g,1:J[g],1:2],
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
