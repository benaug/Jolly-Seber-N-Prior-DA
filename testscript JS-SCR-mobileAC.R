#This is an SCR version with mobile activity centers and an individual survival covariate.
#the movement scale parameter sigma.move generally mixes poorly.
#better with higher cumulative detection prob, higher survival, more years of data, more individuals
#I am estimating a fixed per capita recruitment parameter gamma to better estimate sigma.move
#so, don't simulate different gammas in each year below

library(nimble)
library(coda)
library(truncnorm) #for data simulator
source("sim.JS.SCR.R")
source("Nimble Model JS-SCR-mobileAC.R")
source("Nimble Functions JS-SCR-mobileAC.R") #contains custom distributions and updates
source("sSampler Mobile.R")

n.year <- 5 #number of years
lambda.y1 <- 100 #expected N in year 1
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5 #phi response to individual covariate
p0 <- rep(0.1,n.year) #yearly detection probabilities at activity center
sigma <- rep(0.5,n.year) #yearly detection function scale
sigma.move <- 0.75 #movement sigma, fixed over primary periods
K <- rep(10,n.year) #yearly sampling occasions

buff <- 2 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X[[g]] <- as.matrix(expand.grid(3:11,3:11))
}

data <- sim.JS.SCR(lambda.y1=lambda.y1,gamma=gamma,n.year=n.year,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p0=p0,sigma=sigma,X=X,buff=buff,K=K,sigma.move=sigma.move)

data$N[1] + sum(data$N.recruit) #N.super

##Initialize##
M <- 300 #data augmentation level. Check N.super posterior to make sure it never hits M
N.super.init <- nrow(data$y)
X <- data$X #pull X from data (won't be in environment if not simulated directly above)
if(N.super.init > M) stop("Must augment more than number of individuals captured")
J <- unlist(lapply(X,nrow)) #traps per year
J.max <- max(J)

y.nim <- array(0,dim=c(M,n.year,J.max))
y.nim[1:N.super.init,1:n.year,1:J.max] <- data$y #all these guys must be observed
#initialize z, start with observed guys
z.init <- 1*(y.nim>0)
z.init <- matrix(0,M,n.year)
z.start.init <- z.stop.init <- rep(0,M)
y.nim2D <- apply(y.nim,c(1,2),sum)
for(i in 1:N.super.init){
  det.idx <- which(y.nim2D[i,]>0)
  if(length(det.idx)>0){ #some all 0 histories if init from truth
    z.start.init[i] <- min(det.idx)
    z.stop.init[i] <- max(det.idx)
    z.init[i,z.start.init[i]:z.stop.init[i]] <- 1
  }
}
z.super.init <- c(rep(1,N.super.init),rep(0,M-N.super.init))
z.obs <- 1*(rowSums(y.nim)>0) #indicator for "ever observed"

#initialize N structures from z.init
N.init <- colSums(z.init[z.super.init==1,])
N.survive.init <- N.recruit.init <- rep(NA,n.year-1)
for(g in 2:n.year){
  N.survive.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1)
  N.recruit.init[g-1] <- N.init[g]-N.survive.init[g-1]
}

#now individual covariate
phi.cov.data <- c(data$cov,rep(NA,M-length(data$cov)))
cov.up <- which(is.na(phi.cov.data)) #which individuals have missing cov values, used below to help nimble assign samplers

#remaining SCR stuff to initialize
#put X in ragged array
#also make K1D, year by trap operation history, as ragged array.
X.nim <- array(0,dim=c(n.year,J.max,2))
K1D <- matrix(0,n.year,J.max)
for(g in 1:n.year){
  X.nim[g,1:J[g],1:2] <- X[[g]]
  K1D[g,1:J[g]] <- rep(K[g],J[g])
}

#pull out state space with buffer around maximal trap dimensions
xlim <- data$xlim
ylim <- data$ylim
sigma.move.init <- sigma.move
#initialize s consistent with sigma.move.init
s.init <- initialize.s(sigma.move.init,z.super.init,y=y.nim,X=X.nim,xlim=xlim,ylim=ylim)

#constants for Nimble
constants <- list(n.year=n.year, M=M, J=J, xlim=xlim, ylim=ylim, K1D=K1D)
#inits for Nimble
Niminits <- list(N=N.init,lambda.y1=N.init[1], #initialize consistent with N[1] for faster convergence
                 N.survive=N.survive.init,N.recruit=N.recruit.init,
                 z=z.init,z.start=z.start.init,z.stop=z.stop.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 s=s.init,sigma.move=sigma.move.init,beta0.phi=0,beta1.phi=0,
                 phi.cov.mu=mean(data$truth$cov),phi.cov.sd=sd(data$truth$cov))

#data for Nimble
Nimdata <- list(y=y.nim,phi.cov=phi.cov.data,X=X.nim)

# set parameters to monitor
parameters <- c('N','gamma.fixed','N.recruit','N.survive','N.super',
                'lambda.y1','beta0.phi','beta1.phi',
                'phi.cov.mu','phi.cov.sd','p0','sigma','sigma.move')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#if you add/remove parameters in model file, do so in config.nodes
config.nodes <- c('beta0.phi','beta1.phi','gamma.fixed','lambda.y1',paste('phi.cov[',cov.up,']'),
               'phi.cov.mu','phi.cov.sd','p0','sigma','sigma.move')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = TRUE)

#add N/z samplers
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration? 
#205 of M seems reasonable, but optimal will depend on data set
#loop here bc potentially different numbers of traps to vectorize in each year
y.nodes <- pd.nodes <- c()
for(g in 1:n.year){
  y.nodes <- c(y.nodes,Rmodel$expandNodeNames(paste0("y[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
  pd.nodes <- c(pd.nodes,Rmodel$expandNodeNames(paste0("pd[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
}
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
s.nodes <- Rmodel$expandNodeNames(paste0("s"))
calcNodes <- c(N.nodes,ER.nodes,N.recruit.nodes,N.survive.nodes,s.nodes,z.nodes,pd.nodes,y.nodes)
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,J=J,xlim=xlim,ylim=ylim,
                                                 z.obs=z.obs,z.super.ups=z.super.ups,
                                                 y.nodes=y.nodes,pd.nodes=pd.nodes,N.nodes=N.nodes,
                                                 z.nodes=z.nodes,ER.nodes=ER.nodes,s.nodes=s.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 y2D=y.nim2D,calcNodes=calcNodes), silent = TRUE)

#activity center sampler. There are 2 samplers here for these cases
#1) z.super=1 and z=1, sSampler1 uses Metropolis-Hastings
#2) z.super=1 and z=0, sSampler2 uses Metropolis-Hastings with proposal sd tuned separately from above
#z.super=0, do nothing, s[i,1:n.year,1:2] are set to 0
#N/z sampler proposes activity centers when turning z.super on
for(i in 1:M){
  for(g in 1:n.year){
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler1',control=list(i=i,g=g,xlim=xlim,ylim=ylim,scale=1),silent = TRUE)
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler2',control=list(i=i,g=g,xlim=xlim,ylim=ylim,scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#can replace RW update with slice. nearly as fast as RW in data sets I've tried
#may be smarter when individual s trajectories move in/out of likelihood with z.super updates
#may handle funnel behavior better for smaller values
conf$removeSampler('sigma.move')
conf$addSampler(target='sigma.move',type='slice')

#optional (but recommended!) blocking 
conf$addSampler(target = c("beta0.phi","beta1.phi"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(10000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  #total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 #post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:500),]))
plot(mcmc(mvSamples[-c(1:500),"sigma.move"]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$N[1] + sum(data$N.recruit) #N.super

rem.idx <- c(grep("N",colnames(mvSamples)),
             grep("N.recruit",colnames(mvSamples)),
             grep("N.survive",colnames(mvSamples)))
tmp <- cor(mvSamples[-c(1:500),-rem.idx])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)
round(tmp,2)

# #Some sanity checks I used during debugging. Just checking that final
# #model states match between z and N objects
# 
# #check N count
# N.count <- rep(NA,n.year)
# for(g2 in 1:n.year){
#   N.count[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2]==1) 
# }
# N.count
# Cmodel$N
# 
# #check N.recruit count
# N.count <- rep(NA,n.year-1)
# for(g2 in 2:n.year){
#   N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1) 
# }
# N.count
# Cmodel$N.recruit
# 
# #check N.survive count
# N.count <- rep(NA,n.year-1)
# for(g2 in 2:n.year){
#   N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1) 
# }
# N.count
# Cmodel$N.survive
# 
# #are individual z's consistent with their z.start and z.stop?
# all(apply(Cmodel$z[Cmodel$z.super==1,],1,function(x){min(which(x==1))})==Cmodel$z.start[Cmodel$z.super==1])
# all(apply(Cmodel$z[Cmodel$z.super==1,],1,function(x){max(which(x==1))})==Cmodel$z.stop[Cmodel$z.super==1])
# 
# #zombie check
# for(i in 1:M){
#   z.tmp <- Cmodel$z[i,]
#   z.on <- which(z.tmp==1)
#   first <- z.on[1]
#   last <- z.on[length(z.on)]
#   if(length(z.on)>1){
#     if(any(z.tmp[first:last]==0)){
#       stop("rawr!")
#     }
#   }
# }
