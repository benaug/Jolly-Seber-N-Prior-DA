#This is the Schwarz-Arnason model from Royle and Dorazio 2008 but using direct entry probabilities
#requires custom samplers to update latent variables e and surv
library(nimble)
library(coda)
source("sim.JS.SA.R")
source("Nimble Model JS-SA-marginal.R")
source("Nimble Functions JS-SA-marginal.R") #contains required custom update for e/surv nodes

n.primary <- 4 #number of years
M <- 200 #data simulator simulates from Chandler-Clark model with M as a parameter
psi <- 0.75 #expected N.super is M*psi
pi <- c(0.5,0.2,0.2,0.1) #Probability of entry in occasion g
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.primary) #yearly detection probability
K <- rep(10,n.primary) #yearly sampling occasions

set.seed(33955)
data <- sim.JS.SA(psi=psi,pi=pi,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.primary=n.primary,K=K,M=M)
data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size


##### Initialize p, phi.cov #####
n.det <- nrow(data$y)
y <- rbind(data$y,matrix(0,M-n.det,n.primary))
aug.idx <- seq_len(M - n.det) + n.det

first.det <- last.det <- rep(0,M)
for(i in 1:n.det){
  det.idx <- which(y[i,] > 0)
  first.det[i] <- min(det.idx)
  last.det[i] <- max(det.idx)
}

#detection inits from the known-alive windows. biased high, but far better than
#letting nimble draw p from dunif(0,1)
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

#phi covariate data. latent for undetected inds, so they need inits
phi.cov.data <- rep(NA,M)
phi.cov.data[1:n.det] <- data$phi.cov
phi.cov.init <- rep(NA,M)
phi.cov.init[aug.idx] <- 0

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M)
#inits for Nimble. N, B, N.super are recovery targets, any legal value works
Niminits <- list(psi=0.5,pi=rep(1/n.primary,n.primary),
                 beta0.phi=0,beta1.phi=0,p=p.init,
                 phi.cov=phi.cov.init,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov),
                 N=rep(n.det,n.primary),B=rep(0,n.primary),N.super=n.det)
#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data)
# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','pi','p','phi.cov.mu','phi.cov.sd',"B","N.super")
nt <- 1 #thinning rate
# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('psi','pi','p','beta0.phi','beta1.phi','phi.cov.mu','phi.cov.sd','phi.cov')
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,nodes=config.nodes)

#add the sampler for N, B, N.super
conf$addSampler(target=paste0("N[1:",n.primary,"]"), type=NBSampler,
                control=list(M=M,n.primary=n.primary,y=y,K=K,
                             calcNodes=Rmodel$expandNodeNames(c("N","B","N.super"))))
# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
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
plot(mcmc(mvSamples[-c(1:200),]))

data$N #realized abundance
data$B #realized recruits
data$N.super #realized N.super

