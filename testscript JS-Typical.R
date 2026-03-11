#This is the Chandler and Clark (2014) approach with per capita recruitment
library(nimble)
library(coda)
source("sim.JS.typical.R")
source("Nimble Model JS-Typical.R")

n.year <- 5 #number of years
psi <- 0.33333 #expected N in year 1 is M*psi
gamma <- rep(0.2,n.year) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.year) #yearly detection probability
K <- rep(10,n.year) #yearly sampling occasions

M <- 300

data <- sim.JS.typical(psi=psi,gamma=gamma,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.year=n.year,K=K,M=M)

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
constants<-list(n.year=n.year,K=K,M=M)

#inits for Nimble
Niminits <- list(z=z.init)

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data)

# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','gamma','p')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(1000,reset=FALSE) #can extend run by rerunning this line
end.time<-Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

plot(mcmc(mvSamples[-c(1:200),]))

