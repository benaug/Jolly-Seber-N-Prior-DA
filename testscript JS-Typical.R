#This is the Chandler and Clark (2014) approach with per capita recruitment
library(nimble)
library(coda)
source("sim.JS.typical.R")
source("Nimble Model JS-Typical.R")

n.year <- 4 #number of years
M <- 200 #data simulator simulates from Chandler-Clark model with M as a parameter
psi <- 0.4 #expected N in year 1 is M*psi
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.year) #yearly detection probability
K <- rep(10,n.year) #yearly sampling occasions

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

data <- sim.JS.typical(psi=psi,gamma=gamma,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.year=n.year,K=K,M=M)
data$N #yearly abundance
data$N.super #superpopulation size

#you may need to raise M from what was used to simulate for model fitting.
#monitoring A.raw lets you see if you run out of possible recruits.
#also, gamma.prime needs to be <0.05-0.1 for Poisson recruitment variance approximation to hold within 5-10%
data$gamma.prime #can check for simulated data sets

##### Initialize z using observed data #####
z.init <- matrix(0,M,n.year)
n.super <- nrow(data$y)
for(i in 1:n.super){
  det.idx <- which(data$y[i,]>0)
  z.init[i,min(det.idx):max(det.idx)] <- 1
}

#augment data
y <- rbind(data$y,matrix(0,M-n.super,n.year))

#phi covariate data. nimble can init for undetected inds
phi.cov.data <- rep(NA,M)
phi.cov.data[1:n.super] <- data$phi.cov

#constants for Nimble
constants <- list(n.year=n.year,K=K,M=M)

#inits for Nimble
Niminits <- list(z=z.init,psi=sum(z.init[,1])/M,beta0.phi=0,beta1.phi=0,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov))

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data)

# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','gamma','p','phi.cov.mu','phi.cov.sd')
parameters2 <- c('A.raw','gamma.prime.raw') #monitor these to assess whether M is large enough
nt <- 1 #thinning rate
nt2 <- 1

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,
                      monitors2=parameters2,thin2=nt2)

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(1000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:200),]))

data$N

mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
plot(mcmc(mvSamples2[-c(1:200),]))

#if A.raw hits 0, that means you maxed out M and ran out of possible recruits
#if any gamma.prime.raw > 0.05-0.1, Poisson recruitment approximation is likely not great



