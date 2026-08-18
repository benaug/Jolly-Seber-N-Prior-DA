NimModel <- nimbleCode({
  psi ~ dbeta(1,1) #uniform, but get conjugate sampler

  #occasion-specific unit-time survival and per capita recruitment
  gamma ~ dunif(0,2) #fixed gamma
  for(g in 1:(n.primary-1)){
    phi[g] ~ dbeta(1,1) #unit-time survival
    # gamma[g] ~ dunif(0,2) #gamma varies by primary occasion
    phi.int[g] <- phi[g]^tau[g] #survival over actual interval
  }
  
  #relative expected entry and abundance
  #the absolute scale is supplied separately by psi
  entry.rel[1] <- 1
  EN.rel[1] <- 1
  for(g in 1:(n.primary-1)){
    #linear scaling of per-capita recruitment with interval length
    entry.rel[g+1] <- EN.rel[g]*gamma*tau[g] #fixed gamma
    # entry.rel[g+1] <- EN.rel[g]*gamma[g]*tau[g] #gamma varies by primary occasion
    EN.rel[g+1] <- EN.rel[g]*phi.int[g] + entry.rel[g+1]
  }

  #entry probabilities conditional on superpopulation membership
  entry.denom <- sum(entry.rel[1:n.primary])
  for(g in 1:n.primary){
    pi[g] <- entry.rel[g]/entry.denom
    #absolute expected abundance and entries, derived from M*psi
    EN[g] <- M*psi*EN.rel[g]/entry.denom
    EB[g] <- M*psi*pi[g]
  }

  #conditional entry probabilities
  eta[1] <- pi[1]
  for(g in 2:(n.primary-1)){
    eta[g] <- pi[g]/(1-sum(pi[1:(g-1)]))
  }
  eta[n.primary] <- 1

  #population dynamics
  for(i in 1:M){
    z.super[i] ~ dbern(psi)
    z[i,1] ~ dbern(eta[1])
    a[i,1] <- 1-z[i,1]
    for(g in 2:n.primary){
      mu[i,g] <- phi.int[g-1]*z[i,g-1] + eta[g]*a[i,g-1]
      z[i,g] ~ dbern(mu[i,g])
      a[i,g] <- a[i,g-1]*(1-z[i,g])
    }
  }

  #observation model
  for(g in 1:n.primary){
    p[g] ~ dbeta(1,1)
    for(i in 1:M){
      y.size[i,g] <- K[g]*z[i,g]*z.super[i]
      y[i,g] ~ dbinom(p=p[g],size=y.size[i,g])
    }
  }

  #derived abundance and recruitment
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g]*z.super[1:M])
    Acum[g] <- sum((1-a[1:M,g])*z.super[1:M])
  }
  B[1] <- Acum[1]
  for(g in 2:n.primary){
    B[g] <- Acum[g]-Acum[g-1]
  }
  N.super <- sum(z.super[1:M])
})
