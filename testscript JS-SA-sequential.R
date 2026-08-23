#This is the Schwarz-Arnason model from Royle and Dorazio 2008 using conditional entry probabilities
library(nimble)
library(coda)
source("sim.JS.SA.R")
source("Nimble Model JS-SA-sequential.R")
source("Nimble Functions JS-SA-sequential.R")

n.primary <- 4 #number of primary occasions
M <- 200 #data simulator simulates from Chandler-Clark model with M as a parameter
psi <- 0.75 #expected N.super is M*psi
pi <- c(0.5,0.2,0.2,0.1) #Probability of entry in occasion g
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(10,n.primary) #sampling occasions by primary occasion

set.seed(33955)
data <- sim.JS.SA(psi=psi,pi=pi,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p=p,n.primary=n.primary,K=K,M=M)
data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size


##### Initialize z using observed data #####
z.data <- matrix(NA,M,n.primary)
z.init <- matrix(0,M,n.primary)
z.super.init <- rep(0,M)
n.det <- nrow(data$y)
for(i in 1:n.det){
  det.idx <- which(data$y[i,]>0)
  z.data[i,min(det.idx):max(det.idx)] <- 1
  z.init[i,min(det.idx):max(det.idx)] <- 1
  z.super.init[i] <- 1
}
#every augmented individual must enter somewhere (eta[n.primary]==1),
#so give each one an entry occasion. Spread them across occasions.
aug.idx <- seq_len(M - n.det) + n.det
for(j in seq_along(aug.idx)){
  g <- ((j - 1) %% n.primary) + 1
  z.init[aug.idx[j],g] <- 1
}

#augment data
y <- rbind(data$y,matrix(0,M-n.det,n.primary))

#phi covariate data. nimble can init for undetected inds
phi.cov.data <- rep(NA,M)
phi.cov.data[1:n.det] <- data$phi.cov

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M)

#inits for Nimble
Niminits <- list(z.super=z.super.init,z=z.init,psi=sum(z.super.init)/M,
                 pi=rep(1/n.primary,n.primary),
                 beta0.phi=0,beta1.phi=0,
                 phi.cov.mu=mean(data$phi.cov),phi.cov.sd=sd(data$phi.cov))

#data for Nimble
Nimdata <- list(y=y,phi.cov=phi.cov.data,z=z.data)

# set parameters to monitor
parameters <- c('psi','N','beta0.phi','beta1.phi','pi','p','phi.cov.mu','phi.cov.sd',"B","N.super")
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt)


#z sampler: remove sequential sampler, replace with gibbs
z.nodes <- grep("^z\\[",Rmodel$getNodeNames(stochOnly=TRUE),value=TRUE)
z.super.nodes <- grep("^z\\.super\\[",Rmodel$getNodeNames(stochOnly=TRUE),value=TRUE)
conf$removeSamplers(z.nodes)
conf$removeSamplers(z.super.nodes)
#summarize data for custom update
z.obs <- as.integer(rowSums(y)>0)
first.det <- max.col(y>0,ties.method="first")
last.det <- ncol(y)+1-max.col((y[,ncol(y):1,drop=FALSE]>0),ties.method="first")
first.det[z.obs==0] <- 0
last.det[z.obs==0] <- 0
conf$addSampler(target=c(z.nodes,z.super.nodes),type=zSampler,
                control=list(M=M,K=K,n.primary=n.primary,z.obs=z.obs,
                             first.det=first.det,last.det=last.det))

#add conjugate updates for eta that nimble does not recognize
for(g in 2:(n.primary-1)){
  target <- paste0("eta[",g,"]")
  conf$removeSamplers(target)
  conf$addSampler(target=target,type=etaGibbsSampler,control=list(g=g,M=M,n.primary=n.primary)
  )
}

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc,project=Rmodel)

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

#ESS per time (check time units)
effectiveSize(mcmc(mvSamples[-c(1:200),]))/as.numeric(time2)
