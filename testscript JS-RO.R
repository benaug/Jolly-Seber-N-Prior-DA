#This is the Restricted Occupancy from Royle and Dorazio 2008
library(nimble)
library(coda)
source("sim.JS.RO.R")
source("Nimble Model JS-RO.R")

n.primary <- 4 #number of years
M <- 200 #data simulator simulates from Chandler-Clark model with M as a parameter
psi <- 0.4 #expected N in year 1 is M*psi
gamma <- rep(0.2,n.primary-1) #yearly conditional entry probabilities (must be between 0 and 1)
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.primary) #yearly detection probability
K <- rep(10,n.primary) #yearly sampling occasions

M*psi #expected N[1]

#data simulator will give error messages if M not set large enough
#there are two criteria to consider
#1) you max out M and run out of individuals to recruit. 
#This is very bad and will stop the data simulator and ask you to raise M
#2) you don't have enough possible recruits left for the binomial approximation of
#Poisson recruitment to be accurate. This is more of a grey area. The variance in recruits
#will artificially decrease through time. Generally, M needs to be set much higher than
#we would like for efficient MCMC to achieve this. This is a drawback of this parameterization.
#current parameter settings above do not allow for good Poisson approximation for recruits
#if gamma.prime is 0.05, poisson recruit variance underestimated 5%. if 0.1, 10%, if 0.2, 20%, etc.

data <- sim.JS.RO(psi=psi,gamma=gamma,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.primary=n.primary,K=K,M=M)
data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size

#you may need to raise M from what was used to simulate for model fitting.
#monitoring A.raw lets you see if you run out of possible recruits.
#also, gamma.prime needs to be <0.05-0.1 for Poisson recruitment variance approximation to hold within 5-10%
data$gamma.prime #can check for simulated data sets

##### Initialize z using observed data #####
z.init <- matrix(0,M,n.primary)
n.super <- nrow(data$y)
for(i in 1:n.super){
  det.idx <- which(data$y[i,]>0)
  z.init[i,min(det.idx):max(det.idx)] <- 1
}

#augment data
y <- rbind(data$y,matrix(0,M-n.super,n.primary))

#phi covariate data. nimble can init for undetected inds
phi.cov.data <- rep(NA,M)
phi.cov.data[1:n.super] <- data$phi.cov

#constants for Nimble
#Can choose between two sets of priors if appropriately commented in the model file
#1) Royle and Dorazio 2008 priors
# constants <- list(n.primary=n.primary,K=K,M=M)

#2) Dorazio 2020 (Biometrics 76:1285) priors. Induce a flat (discrete-uniform) prior
#on N.super and equal E[B[g]] across occasions. Uniform gamma priors above do neither:
#they favor recruitment in earlier occasions and bias N.super upwards, worse with
#more occasions and lower p.
#Note: a.gam < 1 (=1/n.primary) puts a spike at gamma=0, so the chain can
#occasionally stall with B[g]=0 (gamma[g] near 0 makes z flips
#unlikely, which keeps gamma[g] near 0). Usually transient.
constants <- list(n.primary=n.primary,K=K,M=M,
                  a.psi=1/n.primary,b.psi=2-1/n.primary,
                  a.gam=rep(1/n.primary,n.primary-1),b.gam=2-(2:n.primary)/n.primary)

#inits for Nimble
Niminits <- list(z=z.init,psi=sum(z.init[,1])/M,beta0.phi=0,beta1.phi=0,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov))

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data)

# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','gamma','p','phi.cov.mu','phi.cov.sd',"B")
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt)

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
#compute N.super posterior and append it to mvSamples
mvSamples <- cbind(mvSamples,rowSums(mvSamples[,grep("B",colnames(mvSamples))])+mvSamples[,"N[1]"])
colnames(mvSamples)[ncol(mvSamples)] <- "N.super"
plot(mcmc(mvSamples[-c(1:200),]))

data$N #realized abundance
data$B #realized recruits
data$N.super #realized N.super
