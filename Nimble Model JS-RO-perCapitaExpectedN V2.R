NimModel <- nimbleCode({
  #probability of ever entering during the study
  psi.super ~ dbeta(1,1) #uniform, but use custom conjugate sampler
  
  #unit-time survival and per capita recruitment
  gamma ~ dunif(0,2) #fixed gamma
  for(g in 1:(n.primary-1)){
    phi[g] ~ dbeta(1,1)
    phi.int[g] <- phi[g]^tau[g] #survival over actual interval
  }
  
  #relative expected entry and abundance
  #the absolute scale is supplied separately by psi.super
  entry.rel[1] <- 1
  EN.rel[1] <- 1
  for(g in 1:(n.primary-1)){
    entry.rel[g+1] <- EN.rel[g]*gamma*tau[g] #expected entries, fixed gamma
    EN.rel[g+1] <- EN.rel[g]*phi.int[g] + entry.rel[g+1]
  }
  
  #unconditional entry probabilities and absolute expectations
  entry.denom <- sum(entry.rel[1:n.primary])
  for(g in 1:n.primary){
    pi[g] <- entry.rel[g]/entry.denom
    beta[g] <- psi.super*pi[g] #unconditional probability of first entry at g
    EN[g] <- M*psi.super*EN.rel[g]/entry.denom
    EB[g] <- M*psi.super*pi[g]
  }
  
  #restricted-occupancy conditional entry probabilities
  psi <- beta[1] #p(present on primary occasion 1)
  cum.beta[1] <- beta[1]
  for(g in 1:(n.primary-1)){
    gamma.RO[g] <- beta[g+1]/(1-cum.beta[g])
    cum.beta[g+1] <- cum.beta[g]+beta[g+1]
  }
  
  #population dynamics
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g])
    Acum[g] <- sum(a[1:M,g])
  }
  for(i in 1:M){
    z[i,1] ~ dbern(psi)
    a[i,1] <- z[i,1]
    for(g in 2:n.primary){
      Ez[i,g-1] <- z[i,g-1]*phi.int[g-1] + (1-a[i,g-1])*gamma.RO[g-1]
      z[i,g] ~ dbern(Ez[i,g-1])
      a[i,g] <- max(a[i,g-1],z[i,g]) #not available to recruit if previously alive
    }
  }
  
  #compute realized recruits
  for(g in 2:n.primary){
    B[g-1] <- Acum[g]-Acum[g-1]
  }
  N.super <- Acum[n.primary]
  
  #observation model
  for(g in 1:n.primary){
    p[g] ~ dbeta(1,1)
    for(i in 1:M){
      y[i,g] ~ dbinom(p=p[g],size=K[g]*z[i,g])
    }
  }
})