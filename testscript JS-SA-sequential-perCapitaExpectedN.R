#This is the Schwarz-Arnason model from Royle and Dorazio 2008 using conditional entry probabilities.
#Recruitment is parameterized as per-capita expected recruitment as a function of expected abundance,
#with linear scaling by interval length. The fitted model uses relative expected abundance internally.
library(nimble)
library(coda)
source("sim.JS.SA.sequentialPerCapitaExpectedN.R")
source("Nimble Model JS-SA-sequential-perCapitaExpectedN.R")
source("Nimble Functions JS-SA-sequential-perCapitaExpectedN.R")

n.primary <- 4 #number of primary occasions
M <- 200 #data augmentation size

primary.time <- c(0,1,3,4) #time of each primary occasion
tau <- diff(primary.time) #length of intervals following occasions 1:(n.primary-1)

lambda1 <- 75 #expected abundance at first primary occasion
phi <- c(0.85,0.80,0.90) #unit-time survival by interval
#model file is set up for fixed gamma per primary occasion, so keep them the same here or change model file
gamma <- c(0.2,0.2,0.2) #per-capita expected recruitment per unit time by interval
p <- rep(0.2,n.primary) #detection probability by primary occasion
K <- rep(5,n.primary) #sampling occasions by primary occasion

#simulate some data
data <- sim.JS.SA(lambda1=lambda1,phi=phi,gamma=gamma,
                  p=p,n.primary=n.primary,K=K,M=M,tau=tau)
data$EN #expected abundance
data$lambda #expected entries
data$psi #expected N.super/M
data$pi #entry probabilities
data$N #realized abundances
data$B #realized entries
data$N.super #realized superpopulation size


##### Initialize z using observed data #####
z.init <- matrix(0,M,n.primary)
z.super.init <- rep(0,M)
n.det <- nrow(data$y)
for(i in 1:n.det){
  det.idx <- which(data$y[i,] > 0)
  z.init[i,min(det.idx):max(det.idx)] <- 1
  z.super.init[i] <- 1
}

#every potential trajectory must enter somewhere because eta[n.primary]==1.
#For augmented individuals, spread initial entry occasions across primary occasions.
aug.idx <- seq_len(M-n.det)+n.det
for(j in seq_along(aug.idx)){
  g <- ((j-1)%%n.primary)+1
  z.init[aug.idx[j],g] <- 1
}

#augment data
y <- rbind(data$y,matrix(0,M-n.det,n.primary))

#constants for Nimble
constants <- list(n.primary=n.primary,K=K,M=M,tau=tau)

#inits for Nimble
Niminits <- list(z.super=z.super.init,z=z.init,
                 psi=max(sum(z.super.init)/M,0.01),
                 phi=rep(0.8,n.primary-1),
                 gamma=0.1, p=rep(0.2,n.primary))

#data for Nimble
Nimdata <- list(y=y)

#set parameters to monitor
parameters <- c('psi','N','phi','gamma',"EB",'p','B','N.super')
nt <- 1 #thinning rate

#Build the model, configure the mcmc, and compile
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

#Correlated posteriors likely if you estimate one gamma per primary occasion,
#AF_slice efficient with a few primary occasions. I'm sure this will become too slow with many occasions
# gamma.nodes <- Rmodel$expandNodeNames("gamma")
# conf$removeSampler(gamma.nodes)
# conf$addSampler(target=gamma.nodes,type='AF_slice',control=list(adaptive=TRUE),silent=TRUE)

#Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc,project=Rmodel)

#Run the model
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time #total time for compilation and fitting
time2 <- end.time-start.time2 #post-compilation run time
time2

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:200),]))

data$EN #expected abundance
data$lambda #expected recruits
data$N #realized abundance
data$B #realized recruits
data$N.super #realized N.super
