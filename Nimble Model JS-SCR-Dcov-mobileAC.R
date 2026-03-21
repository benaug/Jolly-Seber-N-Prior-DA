NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  # D.intercept <- exp(D.beta0)*cellArea
  D.intercept <- D0*cellArea
  #First year density model
  #separate this component so s's do not depend on D.intercept
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells] / pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  ##Abundance##
  lambda.y1 <- D.intercept*pi.denom #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.year){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly abundance
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.year-1)]) #size of superpopulation
  
  #Recruitment
  gamma.fixed ~ dunif(0,2) #sharing gamma across years
  for(g in 1:(n.year-1)){
    gamma[g] <- gamma.fixed
    ER[g] <- N[g]*gamma[g] #yearly expected recruits
    N.recruit[g] ~ dpois(ER[g]) #yearly realized recruits
  }
  
  #Survival individual covariates
  phi.cov.mu ~ dunif(-10, 10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  #RSF coefficient for activity center relocation
  rsf.beta ~ dnorm(0,sd=10)
  sigma.move ~ dunif(0,100) #activity center relocation spatial scale parameter (BVN sd)
  #Resource selection function evaluated across all cells for activity center relocation
  rsf[1:n.cells] <- InSS[1:n.cells]*exp(rsf.beta*D.cov[1:n.cells])
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu, sd=phi.cov.sd)
    #all s set to 0 if not in population, z.super[i]=0
    s[i,1,1:2] ~ dHabYear1(pi.cell=pi.cell[1:n.cells],cells=cells[1:n.cells.x,1:n.cells.y],
                           res=res,dSS=dSS[1:n.cells,1:2],xlim=xlim[1:2],ylim=ylim[1:2],z.super=z.super[i])
    #subsequent year activity center - movement with resource selection
    for(g in 2:n.year){
      #all avail.dist, use.dist, and s set to 0 if not in population, z.super[i]=0
      avail.dist[i,g-1,1:n.cells] <- getAvail(s=s[i,g-1,1:2],sigma=sigma.move,res=res,
                                              x.vals=x.vals[1:n.cells.x],y.vals=y.vals[1:n.cells.y],
                                              n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=z.super[i])
      use.dist[i,g-1,1:n.cells] <- getUse(rsf=rsf[1:n.cells], avail.dist=avail.dist[i,g-1,1:n.cells],z.super=z.super[i])
      s[i,g,1:2] ~ dHabMove(s.prev=s[i,g-1,1:2],use.dist=use.dist[i,g-1,1:n.cells],
                            dSS=dSS[1:n.cells,1:2],cells=cells[1:n.cells.x,1:n.cells.y],
                            res=res,sigma.move=sigma.move,z.super=z.super[i])
    }
  }
  
  #Survival (phi must have M x n.year - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0,sd=10) #individual covariate effect on survival
  for(i in 1:M){
    for(g in 1:(n.year-1)){#plugging same individual phi's into each year for custom update
      logit(phi[i,g]) <- beta0.phi + beta1.phi*phi.cov[i] #individual by year survival
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.year] ~ dSurvival(phi=phi[i,1:(n.year-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
  }
  ##Detection##
  sigma ~ dunif(0,10) #fixing sigma across years
  for(g in 1:n.year){
    p0[g] ~ dunif(0,1) #p0 varies by year
    for(i in 1:M){ #only compute pd and y when z.super[i]=1&z[i,g]=1
      pd[i,g,1:J[g]] <- GetDetectionProb(s=s[i,g,1:2],X=X[g,1:J[g],1:2],
                                         J=J[g],sigma=sigma,p0=p0[g],z=z[i,g],z.super=z.super[i])
      y[i,g,1:J[g]] ~ dBinomialVector(pd=pd[i,g,1:J[g]],K=K1D[g,1:J[g]],
                                      z=z[i,g],z.super=z.super[i]) #vectorized obs mod
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update
#4) s update for z.super=1 only