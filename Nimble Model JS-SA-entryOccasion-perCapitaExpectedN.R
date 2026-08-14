NimModel <- nimbleCode({
  #expected initial population
  lambda[1] ~ dunif(0,1000)
  #occasion-specific unit-time survival and per capita recruitment
  gamma ~ dunif(0,2) #fixed gamma
  for(g in 1:(n.primary-1)){
    phi[g] ~ dbeta(1,1)
    # gamma[g] ~ dunif(0,2) #gamma varies by primary occasion
    phi.int[g+1] <- phi[g]^tau[g] #survival over actual interval
  }
  #expected abundance and expected entries
  EN[1] <- lambda[1]
  for(g in 2:n.primary){
    lambda[g] <- EN[g-1]*gamma*tau[g-1] #expected entries, fixed gamma
    # lambda[g] <- EN[g-1]*gamma[g-1]*tau[g-1] #expected entries, gamma varies by primary occasion
    EN[g] <- EN[g-1]*phi.int[g] + lambda[g]
  }
  lambda.super <- sum(lambda[1:n.primary])
  #required because psi must be <= 1
  lambda.valid ~ dbern(step(M-lambda.super))
  #derived DA parameters
  psi <- lambda.super/M
  for(g in 1:n.primary){
    pi[g] <- lambda[g]/lambda.super
  }
  #population dynamics
  for(i in 1:M){
    z.super[i] ~ dbern(psi)
    e[i] ~ dcat(pi[1:n.primary])
    z[i,1] <- equals(e[i],1)
    recruit[i,1] <- z.super[i]*z[i,1]
    for(g in 2:n.primary){
      surv.p[i,g] <- phi.int[g]*z[i,g-1]
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
  #derived variables
  for(g in 1:n.primary){
    N[g] <- sum(z[1:M,g]*z.super[1:M])
    B[g] <- sum(recruit[1:M,g])
  }
  N.super <- sum(z.super[1:M])
})
