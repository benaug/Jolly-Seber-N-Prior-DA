#This is the Schwarz-Arnason model from Royle and Dorazio 2008.
#Expected recruitment is per capita as a function of *expected* abundance,
#with linear scaling by interval length. The demographic recursion is on a
#relative scale, with psi.super supplying the absolute superpopulation scale.
#Survival and per capita recruitment are occasion-specific unit-time parameters.
#I removed the individual survival covariate here to avoid numerical integration.
#individual covariates on detection can still be used, meaning this would work for SCR.
#requires custom sampler to update latent variables z.super, e, and surv
#this per capita version takes longer to build the model due to more dependencies,
#but runs fast
#The JS models, using N-prior data augmentation, model recruitment as a function of *realized*
#abundance, which makes more sense to me. 

library(nimble)
library(coda)
source("sim.JS.SA.perCapitaExpectedN.R")
source("Nimble Model JS-SA-entryOccasion-perCapitaExpectedN.R")
source("Nimble Functions JS-SA-entryOccasion-perCapitaExpectedN.R") #contains required custom update for z.super/e/surv nodes

n.primary <- 4 #number of primary occasions
M <- 200 #data augmentation size

primary.time <- c(0,1,3,4) #time of each primary occasion
tau <- diff(primary.time) #length of intervals preceding occasions 2:n.primary

lambda1 <- 75 #expected initial population used to simulate data
#model file is set up for fixed gamma, so keep them the same here or change model file
gamma <- c(0.20,0.20,0.20) #per capita expected recruitment per unit time by interval
phi <- c(0.85,0.80,0.90) #unit-time survival by interval
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(10,n.primary) #sampling occasions by primary occasion

#simulate some data
data <- sim.JS.SA.perCapitaExpectedN(lambda1=lambda1,gamma=gamma,phi=phi,
                                     p=p,n.primary=n.primary,K=K,M=M,tau=tau)

data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size

##### Initialize e, surv, z.super
n.det <- nrow(data$y)
y <- rbind(data$y,matrix(0,M-n.det,n.primary))
first.det <- last.det <- rep(0,M)
e.init <- rep(NA,M)
surv.init <- matrix(NA,M,n.primary)
z.super.init <- rep(0,M)

#detected individuals: enter at first detection, alive through last
for(i in 1:n.det){
  det.idx <- which(y[i,] > 0)
  first.det[i] <- min(det.idx)
  last.det[i] <- max(det.idx)
  e.init[i] <- first.det[i]
  z.super.init[i] <- 1
  for(g in 2:n.primary){
    surv.init[i,g] <- as.numeric(g > first.det[i] & g <= last.det[i])
  }
}

#augmented individuals: spread entry occasions, alive only at entry
aug.idx <- seq_len(M - n.det) + n.det
for(j in seq_along(aug.idx)){
  i <- aug.idx[j]
  e.init[i] <- ((j - 1) %% n.primary) + 1
  for(g in 2:n.primary){
    surv.init[i,g] <- 0
  }
}

p.num <- rep(0,n.primary)
p.den <- rep(0,n.primary)
for(i in 1:n.det){
  for(g in first.det[i]:last.det[i]){
    p.num[g] <- p.num[g] + y[i,g]
    p.den[g] <- p.den[g] + K[g]
  }
}
p.init <- p.num/p.den
p.init <- pmin(pmax(p.init,0.01),0.99)

#time intervals
tau <- diff(primary.time)

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M,tau=tau)

#initial values
psi.super.init <- max(sum(z.super.init)/M,0.01)
gamma.init <- 0.1
phi.init <- rep(0.8,n.primary-1)

#inits for Nimble
Niminits <- list(z.super=z.super.init,e=e.init,surv=surv.init,
                 psi.super=psi.super.init,gamma=gamma.init,
                 phi=phi.init,p=p.init)

#data for Nimble
Nimdata <- list(y=y)

#set parameters to monitor
parameters <- c('psi.super','gamma','phi','EN','EB','N','pi','p','B','N.super')
nt <- 1 #thinning rate

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
#don't configure e, surv, z.super
config.nodes <- c('psi.super','gamma','phi','p')
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,nodes=config.nodes)

#add z.super/e/surv sampler
z.obs <- as.numeric(rowSums(y) > 0)
first.det <- last.det <- rep(0,M)
for(i in 1:M){
  if(z.obs[i]==1){
    dets <- which(y[i,] > 0)
    first.det[i] <- min(dets)
    last.det[i] <- max(dets)
  }
}
z.super.nodes <- Rmodel$expandNodeNames("z.super")
e.nodes <- Rmodel$expandNodeNames("e")
surv.nodes <- Rmodel$expandNodeNames("surv")
calcNodes <- Rmodel$getDependencies(c("z.super","e","surv"))

conf$addSampler(target=c(z.super.nodes,e.nodes,surv.nodes),type=eSampler,
                control=list(M=M,K=K,n.primary=n.primary,z.obs=z.obs,
                             first.det=first.det,last.det=last.det,
                             calcNodes=calcNodes))

#Build and compile
Rmcmc <- buildMCMC(conf)
#runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time
time2

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:200),]))

data$N #realized abundance
data$B #realized recruits
data$N.super #realized N.super
gamma #per capita expected recruitment rates
phi #unit-time survival probabilities

#ESS per time (check time units)
effectiveSize(mcmc(mvSamples[-c(1:200),]))/as.numeric(time2)

