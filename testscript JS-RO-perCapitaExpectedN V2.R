#This is the Restricted Occupancy model from Royle and Dorazio 2008,
#but parameterized so expected recruitment is per capita as a function of *expected*
#abundance, with linear scaling by interval length.
#The fitted model uses relative expected abundance internally, with psi.super supplying
#the absolute scale as the probability of ever entering during the study.
#I removed the individual survival covariate here to avoid numerical integration.
#The N-prior data augmentation models use per capita recruitment as a function of
#*realized* abundance, which makes more sense to me.

library(nimble)
library(coda)
source("sim.JS.RO.perCapitaExpectedN.R")
source("Nimble Model JS-RO-perCapitaExpectedN V2.R")
source("Nimble Functions JS-RO-perCapitaExpectedN V2.R")

n.primary <- 4 #number of primary occasions
M <- 200 #data augmentation size

primary.time <- c(0,1,3,4) #time of each primary occasion
tau <- diff(primary.time) #length of intervals following occasions 1:(n.primary-1)

lambda1 <- 75 #expected initial population used to simulate data
#model file is set up for fixed gamma, so keep these the same or change model file
gamma <- c(0.20,0.20,0.20) #per capita expected recruitment per unit time by interval
phi <- c(0.85,0.80,0.90) #unit-time survival by interval
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(10,n.primary) #sampling occasions by primary occasion

set.seed(33955)
data <- sim.JS.RO.perCapitaExpectedN(lambda1=lambda1,gamma=gamma,phi=phi,
                                     p=p,n.primary=n.primary,K=K,M=M,tau=tau)
data$EN #expected abundance
data$lambda #expected entries
data$lambda.super/M #probability of ever entering
data$N #realized abundances
data$B #realized recruits entering occasions 2:n.primary
data$N.super #realized superpopulation size

##### Initialize z using observed data #####
n.det <- nrow(data$y)
y <- rbind(data$y,matrix(0,M-n.det,n.primary))
z.init <- matrix(0,M,n.primary)
#detected individuals: alive from first through last detection
for(i in 1:n.det){
  det.idx <- which(y[i,]>0)
  z.init[i,min(det.idx):max(det.idx)] <- 1
}

p.num <- rep(0,n.primary)
p.den <- rep(0,n.primary)
for(i in 1:n.det){
  det.idx <- which(y[i,]>0)
  for(g in min(det.idx):max(det.idx)){
    p.num[g] <- p.num[g]+y[i,g]
    p.den[g] <- p.den[g]+K[g]
  }
}
p.init <- p.num/p.den
p.init <- pmin(pmax(p.init,0.01),0.99)

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M,tau=tau)

#initial values
psi.super.init <- max(n.det/M,0.01)
gamma.init <- 0.1
phi.init <- rep(0.8,n.primary-1)

#inits for Nimble
Niminits <- list(z=z.init,psi.super=psi.super.init,gamma=gamma.init,
                 phi=phi.init,p=p.init)

#data for Nimble
Nimdata <- list(y=y)

#set parameters to monitor
parameters <- c('psi.super','gamma','phi','phi.int','pi','EN','EB','N','p',
                'B','N.super','psi','gamma.RO')
nt <- 1 #thinning rate

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt)

#z sampler: remove sequential sampler, replace with Wu block Gibbs sampler
z.nodes <- grep("^z\\[",Rmodel$getNodeNames(stochOnly=TRUE),value=TRUE)
conf$removeSamplers(z.nodes)

#summarize data for custom update
z.obs <- as.integer(rowSums(y)>0)
first.det <- max.col(y>0,ties.method="first")
last.det <- ncol(y)+1-max.col((y[,ncol(y):1,drop=FALSE]>0),ties.method="first")
first.det[z.obs==0] <- 0
last.det[z.obs==0] <- 0

conf$addSampler(target=z.nodes,type=zSampler,
                control=list(M=M,K=K,n.primary=n.primary,z.obs=z.obs,
                             first.det=first.det,last.det=last.det))

#psi.super is conjugate conditional on the complete entry histories, but NIMBLE
#does not recognize the conjugacy through the restricted-occupancy parameterization
conf$removeSamplers("psi.super")
conf$addSampler(target="psi.super",type=psi.superSampler,control=list(M=M))

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

data$EN #expected abundance
data$lambda #expected recruits
data$N #realized abundance
data$B #realized recruits
data$N.super #realized N.super
gamma #per capita expected recruitment rates
phi #unit-time survival probabilities

#ESS per time (check time units)
effectiveSize(mcmc(mvSamples[-c(1:200),]))/as.numeric(time2)
