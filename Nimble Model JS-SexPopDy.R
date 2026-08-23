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
  #phi.sex[1,1] ~ dbeta(1,1) #male survival
  #phi.sex[2,1] ~ dbeta(1,1) #female survival
  #for(g in 2:(n.primary-1)){
  #  phi.sex[1,g] <- phi.sex[1,1]
  #  phi.sex[2,g] <- phi.sex[2,1]
  #}
  #OPTION 4: sex-specific phi varies by primary occasion
  for(g in 1:(n.primary-1)){
    phi.sex[1,g] ~ dbeta(1,1) #male survival
    phi.sex[2,g] ~ dbeta(1,1) #female survival
  }
  for(i in 1:M){
    for(g in 1:(n.primary-1)){
      phi[i,g] <- phi.sex[sex[i]+1,g]
    }
    z[i,1:n.primary] ~ dSurvival(phi=phi[i,1:(n.primary-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
  }
  
  ##Observation model - 
  #These options are structured so that the custom Gibbs sampler for p will work in all cases.
  #OPTION 1: one p shared by sexes and primary occasions
  #p.sex[1,1] ~ dbeta(1,1)
  #p.sex[2,1] <- p.sex[1,1]
  #for(g in 2:n.primary){
  #  p.sex[1,g] <- p.sex[1,1]
  #  p.sex[2,g] <- p.sex[1,1]
  #}
  #OPTION 2: fixed sex-specific p
  #p.sex[1,1] ~ dbeta(1,1) #male p
  #p.sex[2,1] ~ dbeta(1,1) #female p
  #for(g in 2:n.primary){
  #  p.sex[1,g] <- p.sex[1,1]
  #  p.sex[2,g] <- p.sex[2,1]
  #}
  #OPTION 3: occasion-specific p shared by sexes
  #for(g in 1:n.primary){
  #  p.sex[1,g] ~ dbeta(1,1)
  #  p.sex[2,g] <- p.sex[1,g]
  #}
  #OPTION 4: sex-specific p varies by primary occasion
  for(g in 1:n.primary){
    p.sex[1,g] ~ dbeta(1,1) #male p
    p.sex[2,g] ~ dbeta(1,1) #female p
  }
  for(g in 1:n.primary){
    for(i in 1:M){
      #must use this custom distribution for custom updates
      y[i,g] ~ dbinomial2(p=p.sex[sex[i]+1,g],K=K[g],z=z[i,g],z.super=z.super[i])
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update
#4) sex update: update sex-specific N structures with sex[i] update
