#This is the Chandler and Clark (2014) approach with per capita recruitment. 
#Wu-type whole-history Gibbs update. So much faster
library(nimble)
library(coda)
source("sim.JS.CC.R")
source("Nimble Model JS-CC.R")
source("Nimble Functions JS-CC.R")

n.primary <- 4 #number of primary occasions
M <- 300 #data simulator simulates from Chandler-Clark model with M as a parameter
psi <- 0.25 #expected N in primary occasion 1 is M*psi
#model file currently set up with fixed gamma, so don't vary them in simulation 
#without updating model file and inits bel0ow
gamma <- rep(0.2,n.primary-1) #per-capita recruitment by primary occasion
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
tau <- rep(1,n.primary-1) #duration of each primary-occasion interval
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(10,n.primary) #sampling occasions by primary occasion

M*psi #expected N[1]

#data simulator will give error messages if M not set large enough
#there are two criteria to consider
#1) you max out M and run out of individuals to recruit.
#This is very bad and will stop the data simulator and ask you to raise M
#2) finite M makes recruitment conditionally binomial rather than Poisson.
#Conditional on N, gamma, and the number available to recruit, the recruitment
#mean is N*gamma, but the variance is N*gamma*(1-gamma.prime). Thus, relative
#to a Poisson distribution with the same conditional mean, gamma.prime=0.05
#gives 5% lower conditional variance, gamma.prime=0.1 gives 10% lower variance, etc.
#This conditional comparison does not directly imply the same reduction in the
#posterior variance of recruitment when gamma and abundance are estimated.
#Larger M reduces gamma.prime and makes the conditional recruitment distribution
#closer to Poisson, but may reduce computational efficiency.
data <- sim.JS.CC(psi=psi,gamma=gamma,tau=tau,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.primary=n.primary,K=K,M=M)

data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size

#you may need to raise M from what was used to simulate for model fitting.
#monitoring A.raw lets you see whether the augmented pool is being exhausted.
#gamma.prime measures the departure of the conditional recruitment distribution
#from its Poisson limit: conditional variance is multiplied by (1-gamma.prime).
#This does not directly measure distortion of the marginal posterior because
#uncertainty in gamma, N, and A also contributes to posterior variation.
#Sensitivity of posterior inference to M is therefore more informative than
#gamma.prime alone.
data$gamma.prime #can check conditional Poisson approximation for simulated data sets

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
constants <- list(n.primary=n.primary,K=K,M=M,tau=data$tau)

#inits for Nimble
Niminits <- list(z=z.init,psi=sum(z.init[,1])/M,beta0.phi=0,beta1.phi=0,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov))

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data,z=z.data)

# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','gamma.fixed','p','phi.cov.mu','phi.cov.sd',"B")
parameters2 <- c('A.raw','gamma.prime.raw') #monitor these to assess whether M is large enough
nt <- 1 #thinning rate
nt2 <- 1 #thinning rate for parameters 2

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt,
                      monitors2=parameters2,thin2=nt2)

#replace single-site z updates with Wu-type whole-history Gibbs sampler
z.nodes <- grep("^z\\[",Rmodel$getNodeNames(stochOnly=TRUE),value=TRUE)
conf$removeSamplers(z.nodes)
conf$addSampler(target=z.nodes,type=zSampler,
                control=list(M=M,K=K,n.primary=n.primary,tau=tau))

#if gamma varies by occasion, can remove nimble-assigned RW samplers (nimble model currently has gamma fixed)
#and replace with full conditionals. If gamma is fixed, this requires rejection sampling
# which may not be more efficient. Didn't create one.
# Note: if you add time scaling to model file, need to include that in custom update
# for(g in 1:(n.primary-1)){
#   target <- paste0("gamma[",g,"]")
#   conf$removeSamplers(target)
#   conf$addSampler(target=target,type=gammaCCSampler,control=list(M=M,n.primary=n.primary,
#                                                                  tau=data$tau,qcap=0.999))
# }

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

mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
plot(mcmc(mvSamples2[-c(1:200),]))

#Things to check for finite-M effects
#if A.raw hits 0, you maxed out M and ran out of possible recruits
#if gamma.prime.raw reaches 1, the model hits the artificial recruitment-probability cap
#larger gamma.prime.raw means a poorer conditional Poisson approximation:
#conditional recruitment variance is N*gamma*(1-gamma.prime)
#however, this does not imply an equivalent reduction in posterior recruitment variance
#because gamma, N, and A are estimated rather than fixed
#the most direct check for consequential finite-M effects is sensitivity of posterior
#inference to increasing M


