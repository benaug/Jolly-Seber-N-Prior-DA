#This is the Restricted Occupancy from Royle and Dorazio 2008 with Dorazio (2020 priors)
#Dorazio 2020 (Biometrics 76:1285) priors. Induce a flat (discrete-uniform) prior
#on N.super and equal E[B[g]] across occasions. Uniform gamma priors above do neither:
#they favor recruitment in earlier occasions and bias N.super upwards, worse with
#more occasions and lower p.
#Note: a.gam < 1 (=1/n.primary) puts a spike at gamma=0, so the chain can
# get stuck at B[g]=0 
#this version didn't do well in a simulation study. I need to look into it more.

library(nimble)
library(coda)
source("sim.JS.RO.R")
source("Nimble Model JS-RO-Dorazio.R")
source("Nimble Functions JS-RO-Dorazio.R")

n.primary <- 4 #number of primary occasions
M <- 200 #data simulator simulates from Chandler-Clark model with M as a parameter
psi <- 0.4 #expected N in primary occasion 1 is M*psi
gamma <- rep(0.2,n.primary-1) #conditional entry probabilities by primary occasion (must be between 0 and 1)
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(10,n.primary) #sampling occasions by primary occasion

M*psi #expected N[1]

set.seed(33955)
data <- sim.JS.RO(psi=psi,gamma=gamma,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.primary=n.primary,K=K,M=M)
data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size

##### Initialize z using observed data #####
z.data <- matrix(NA,M,n.primary)
z.init <- matrix(0,M,n.primary)
n.det <- nrow(data$y)
for(i in 1:n.det){
  det.idx <- which(data$y[i,]>0)
  z.data[i,min(det.idx):max(det.idx)] <- 1
  z.init[i,min(det.idx):max(det.idx)] <- 1
}

#augment data
y <- rbind(data$y,matrix(0,M-n.det,n.primary))

#phi covariate data. nimble can init for undetected inds
phi.cov.data <- rep(NA,M)
phi.cov.data[1:n.det] <- data$phi.cov

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M,
                  a.psi=1/n.primary,b.psi=2-1/n.primary,
                  a.gam=rep(1/n.primary,n.primary-1),b.gam=2-(2:n.primary)/n.primary)

#inits for Nimble
Niminits <- list(z=z.init,psi=sum(z.init[,1])/M,beta0.phi=0,beta1.phi=0,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov))

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data,z=z.data)

# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','gamma','p','phi.cov.mu','phi.cov.sd',"B")
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
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

#beta-Bernoulli Gibbs samplers for gamma
#if occasion-specific
for(g in 1:(n.primary-1)){
  target <- paste0("gamma[",g,"]")
  conf$removeSamplers(target)
  conf$addSampler(target=target,type=gammaGibbsSampler,
                  control=list(M=M,n.primary=n.primary))
}
#if gamma fixed
# conf$removeSamplers("gamma")
# conf$addSampler(target="gamma",type=gammaGibbsSampler,
#                 control=list(M=M,n.primary=n.primary))

#Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time
time2

mvSamples <- as.matrix(Cmcmc$mvSamples)
#compute N.super posterior and append it to mvSamples
mvSamples <- cbind(mvSamples,rowSums(mvSamples[,grep("B",colnames(mvSamples))])+mvSamples[,"N[1]"])
colnames(mvSamples)[ncol(mvSamples)] <- "N.super"
plot(mcmc(mvSamples[-c(1:200),]))

data$N #realized abundance
data$B #realized recruits
data$N.super #realized N.super

#ESS per time (check time units)
effectiveSize(mcmc(mvSamples[-c(1:200),]))/as.numeric(time2)
