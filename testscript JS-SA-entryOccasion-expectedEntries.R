#This is the Schwarz-Arnason model from Royle and Dorazio 2008,
#but parameterized in terms of expected entries. This example shows
#how expected entries during each interval scale with interval length
#while holding the expected number of entries per unit time, rho, fixed.
#requires custom samplers to update latent variables e and surv
library(nimble)
library(coda)
source("sim.JS.SA.expectedEntries.R")
source("Nimble Model JS-SA-entryOccasion-expectedEntries.R")
source("Nimble Functions JS-SA-entryOccasion.R") #contains required custom update for e/surv nodes

n.primary <- 4 #number of primary occasions
M <- 200 #data simulator simulates from Chandler-Clark model with M as a parameter

primary.time <- c(0,1,3,4) #time of each primary occasion
tau <- diff(primary.time) #length of intervals preceding occasions 2:n.primary

lambda1 <- 75 #expected initial population
rho <- 18.75 #expected entries per unit time
lambda <- c(lambda1,rho*tau) #expected entries by primary occasion

beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(10,n.primary) #sampling occasions by primary occasion

set.seed(33955)
data <- sim.JS.SA.expectedEntries(lambda=lambda,
                  beta0.phi=beta0.phi,beta1.phi=beta1.phi,
                  p=p,n.primary=n.primary,K=K,M=M)
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

#augment data
y <- rbind(data$y,matrix(0,M-n.det,n.primary))

#phi covariate data. nimble can init for undetected inds
phi.cov.data <- rep(NA,M)
phi.cov.data[1:n.det] <- data$phi.cov

#time intervals
tau <- diff(primary.time)

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M,tau=tau)

#initial expected entries
lambda.super.init <- sum(z.super.init)
lambda1.init <- lambda.super.init/2
rho.init <- (lambda.super.init-lambda1.init)/sum(tau)

#inits for Nimble
Niminits <- list(z.super=z.super.init,e=e.init,surv=surv.init,
                 lambda=lambda1.init,rho=rho.init,
                 beta0.phi=0,beta1.phi=0,p=p.init,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov))

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data,lambda.valid=1)

#set parameters to monitor
parameters <- c('lambda.super','rho','N','beta0.phi','beta1.phi','lambda','p',
                'phi.cov.mu','phi.cov.sd',"B","N.super")
nt <- 1 #thinning rate

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)

#don't configure e or surv
config.nodes <- c('lambda[1]','rho','p','beta0.phi','beta1.phi',
                  'phi.cov.mu','phi.cov.sd','phi.cov','z.super')
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,nodes=config.nodes)

#add e/surv sampler
z.obs <- as.numeric(rowSums(y) > 0)
first.det <- last.det <- rep(0,M)
for(i in 1:M){
  if(z.obs[i]==1){
    dets <- which(y[i,] > 0)
    first.det[i] <- min(dets)
    last.det[i] <- max(dets)
  }
}
e.nodes <- Rmodel$expandNodeNames("e")
z.nodes <- Rmodel$expandNodeNames("z")
surv.p.nodes <- Rmodel$expandNodeNames("surv.p")
surv.nodes <- Rmodel$expandNodeNames("surv")
y.size.nodes <- Rmodel$expandNodeNames("y.size")
y.nodes <- Rmodel$expandNodeNames("y")
N.B.nodes <- Rmodel$expandNodeNames(c("recruit","N","B","N.super"))
calcNodes <- Rmodel$getDependencies(c("e","surv"))

conf$addSampler(target=paste0("e[1:",M,"]"), type=eSampler,
                control=list(M=M,K=K,n.primary=n.primary,z.obs=z.obs,
                             first.det=first.det,last.det=last.det,
                             e.nodes=e.nodes,z.nodes=z.nodes,
                             surv.p.nodes=surv.p.nodes,surv.nodes=surv.nodes,
                             y.size.nodes=y.size.nodes,y.nodes=y.nodes,
                             N.B.nodes=N.B.nodes,
                             calcNodes=calcNodes))

#Correlated posteriors. AF_slice is relatively cheap here.
conf$removeSampler("lambda[1]","rho")
conf$addSampler(target = c("lambda[1]","rho"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)

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
rho #expected recruitment rate

