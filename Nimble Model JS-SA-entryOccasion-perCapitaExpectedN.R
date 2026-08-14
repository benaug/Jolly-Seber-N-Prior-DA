NimModel <- nimbleCode({
  #probability of being in the superpopulation
  psi.super ~ dbeta(1,1) #uniform, but get conjugate sampler
  
  #occasion-specific unit-time survival and per capita recruitment
  gamma ~ dunif(0,2) #fixed gamma
  for(g in 1:(n.primary-1)){
    phi[g] ~ dbeta(1,1)
    # gamma[g] ~ dunif(0,2) #gamma varies by primary occasion
    phi.int[g+1] <- phi[g]^tau[g] #survival over actual interval
  }
  
  #relative expected abundance and expected entries
  entry.rel[1] <- 1
  EN.rel[1] <- 1
  for(g in 1:(n.primary-1)){
    entry.rel[g+1] <- EN.rel[g]*gamma*tau[g] #expected entries, fixed gamma
    # entry.rel[g+1] <- EN.rel[g]*gamma[g]*tau[g] #expected entries, gamma varies by primary occasion
    EN.rel[g+1] <- EN.rel[g]*phi.int[g+1] + entry.rel[g+1]
  }
  
  #entry probabilities and absolute expectations
  entry.denom <- sum(entry.rel[1:n.primary])
  for(g in 1:n.primary){
    pi[g] <- entry.rel[g]/entry.denom
    EB[g] <- M*psi.super*pi[g]
    EN[g] <- M*psi.super*EN.rel[g]/entry.denom
  }
  
  #population dynamics
  for(i in 1:M){
    z.super[i] ~ dbern(psi.super)
    e[i] ~ dcat(pi[1:n.primary]) #entry occasion
    z[i,1] <- equals(e[i],1)
    recruit[i,1] <- z.super[i]*z[i,1]
    for(g in 2:n.primary){
      surv.p[i,g] <- phi.int[g]*z[i,g-1] #keep separate from line below for custom updates
      surv[i,g] ~ dbern(surv.p[i,g]) #survival outcome
      z[i,g] <- max(surv[i,g],equals(e[i],g)) #alive if survived or entering now
      recruit[i,g] <- z.super[i]*equals(e[i],g) #used to compute B
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
  
  #derived variables
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g]*z.super[1:M])
    B[g] <- sum(recruit[1:M,g])
  }
  N.super <- sum(z.super[1:M])
})
