#This is an SCR version with spatial density covariates, RSF activity center movement,
#and an individual survival covariate. Activity centers are in continuous space.
#Activity center relocation uses a proper resource selection model with an availability 
#distribution that is bivariate normal. Could use other availability distributions, bivariate t, etc.

library(nimble)
library(coda)
source("sim.JS.SCR.Dcov.mobileAC.R")
source("Nimble Model JS-SCR-Dcov-mobileAC.R")
source("Nimble Functions JS-SCR-Dcov-mobileAC.R") #contains custom distributions and updates
source("sSampler Dcov MobileAC.R") # required activity center samplers
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

n.year <- 5 #number of years
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5 #phi response to individual covariate
p0 <- rep(0.1,n.year) #yearly detection probabilities at activity center
sigma <- rep(0.75,n.year) #yearly detection function scale
sigma.move <- 2 #movement sigma, fixed over primary periods
rsf.beta <- 0.5 #selection coefficient for activity center relocation btwn primary periods
K <- rep(10,n.year) #yearly sampling occasions

buff <- 3 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X[[g]] <- as.matrix(expand.grid(3:14,3:14))
}

### Habitat Covariate stuff###
#buffer maximal trap extent
X.all <- matrix(NA,nrow=0,ncol=2)
for(g in 1:n.year){
  X.all <- rbind(X.all,X[[g]])
}

#get x and y extent by buffering state space
xlim <- range(X.all[,1])+c(-buff,buff)
ylim <- range(X.all[,2])+c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
for(g in 1:n.year){
  X[[g]][,1] <- X[[g]][,1]-x.shift
  X[[g]][,2] <- X[[g]][,2]-y.shift
}
X.all[,1] <- X.all[,1]-x.shift
X.all[,2] <- X.all[,2]-y.shift

res <- 0.50 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

#create a density covariate
#simulate a D.cov, higher cov.pars for large scale cov
set.seed(132051)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(25,25),messages=FALSE)[[2]]
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Covariate Value",col=cols1)
points(X.all,pch=4,cex=0.75,col="darkred",lwd=2)

#Additionally, maybe we want to exclude "non-habitat"
#in first year D model
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<(xlim[1]+1)&dSS.tmp[,2]<(ylim[1]+1)] <- 0
InSS[dSS.tmp[,1]<(xlim[1]+1)&dSS.tmp[,2]>=(ylim[2]-1)] <- 0
InSS[dSS.tmp[,1]>=(xlim[2]-1)&dSS.tmp[,2]<(ylim[1]+1)] <- 0
InSS[dSS.tmp[,1]>=(xlim[2]-1)&dSS.tmp[,2]>=(ylim[2]-1)] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -1.25
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

#simulate some data
data <- sim.JS.SCR.Dcov.mobileAC(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
            gamma=gamma,n.year=n.year,beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p0=p0,sigma=sigma,sigma.move=sigma.move,rsf.beta=rsf.beta,
            X=X,K=K,xlim=xlim,ylim=ylim,res=res)

#these plots are cool, you should look at these.

#visualize expected relative density and realized activity centers in each year
#year 1: cell colors depict expected relative density in year 1
#years 2 on: cell colors depict expected relative density given the realized s[g-1] and z[g-1] from previous year
#points are realized activity centers for this expectation
# par(mfrow=c(1,1),ask=FALSE)
# for(plot.year in 1:n.year){
#   image(x.vals,y.vals,matrix(data$truth$pi.cell[plot.year,],n.cells.x,n.cells.y),
#         main=paste("Expected Relative Density, Year", plot.year),
#         col=cols1)
#   points(X.all,pch=4,cex=0.75)
#   points(data$truth$s[data$truth$z[,plot.year]==1,plot.year,],pch=16)
# }

#visualize individual movement trajectories. Start year is larger circle
#will be a mess with a lot of individuals.
# ind.cols <- c("#E63946","#FF9F1C","#FFDD00","#2EC4B6","#3A86FF",
#               "#8338EC","#FB5607","#06D6A0","#FFB700","#118AB2",
#               "#EF476F","#FFC8DD","#B5E48C","#00BBF9","#9B5DE5",
#               "#F15BB5","#00F5D4","#FEE440","#FF595E","#6A4C93")
# n.colors <- length(ind.cols)
# par(mfrow=c(1,1), ask=FALSE)
# image(x.vals, y.vals, matrix(D.cov, n.cells.x, n.cells.y), col=cols1,
#       main="Movement Trajectories over D.cov")
# for(i in 1:data$truth$N.super){
#   #skip individuals alive for less than 2 years
#   if(sum(data$truth$z[i,])<2) next
#   ind.col <- ind.cols[(i-1) %% n.colors + 1]
#   alive.years <- which(data$truth$z[i,]==1) #get years alive
#   # plot points for all alive years
#   points(data$truth$s[i,alive.years,1],
#          data$truth$s[i,alive.years,2],
#          pch=16, col=ind.col, cex=0.8)
# 
#   #plot lines between consecutive alive years
#   if(length(alive.years)>1){
#     for(t in 1:(length(alive.years)-1)){
#       if(alive.years[t+1] == alive.years[t]+1){
#         lines(x=c(data$truth$s[i,alive.years[t],1],data$truth$s[i,alive.years[t+1],1]),
#               y=c(data$truth$s[i,alive.years[t],2],data$truth$s[i,alive.years[t+1],2]),
#               col=ind.col,lwd=1.2)
#       }
#     }
#   }
#   #mark entry point with larger circle
#   points(data$truth$s[i,alive.years[1],1],data$truth$s[i,alive.years[1],2],
#          pch=16,col=ind.col,cex=1.5)
# }
# points(X.all,pch=4,cex=0.75,lwd=2)

# can look at individual by year availability and use distributions
# i <- 1
# g <- 1
# par(mfrow=c(3,1))
# image(x.vals,y.vals,matrix(data$truth$avail.dist[i,g,],n.cells.x,n.cells.y),main="Availability Distribution")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="lightblue")
# image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="RSF Cov")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="lightblue")
# image(x.vals,y.vals,matrix(data$truth$use.dist[i,g,],n.cells.x,n.cells.y),main="Use Distribution")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="lightblue")
# par(mfrow=c(1,1))

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

data$truth$N.super #N.super

##Initialize##
#Hard to predict appropriate M, depends on many factors like detection prob, number of years
#level of population turnover. Maybe make sure it is at least 1.6*N.super to start
M <- 250 #data augmentation level. Check N.super posterior to make sure it never hits M
N.super.init <- nrow(data$y)
X <- data$X #pull X from data (won't be in environment if not simulated directly above)
K <- data$K #same for K
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

#pull these out of data object (won't be in environment if it wasn't simulated above, i.e. real data)
xlim <- data$xlim
ylim <- data$ylim
dSS <- data$dSS
cells <- data$cells
res <- data$res
cellArea <- res^2
D.cov <- data$D.cov
InSS <- data$InSS
x.vals <- data$x.vals
y.vals <- data$y.vals
n.cells <- data$n.cells
n.cells.x <- data$n.cells.x
n.cells.y <- data$n.cells.y

#initialize s consistent with these inits
sigma.move.init <- sigma.move
rsf.beta.init <- 0
#D.beta1.init is only used to create inits for undetected individuals
#but if we init z.super to only be 1 for detected guys as is set up above
#it will not be used
D.beta1.init <- 0
s.init <- initialize.s.hab(sigma.move.init=sigma.move.init,rsf.beta.init=rsf.beta.init,
                            z.super.init=z.super.init,D.beta1.init=D.beta1.init,
                            y=y.nim,X=X.nim,xlim=xlim,ylim=ylim,dSS=dSS,
                            cells=cells,res=res,D.cov=D.cov,InSS=InSS,x.vals=x.vals,y.vals=y.vals)

#can verify all s start inside habitat mask
image(data$x.vals,data$y.vals,matrix(data$InSS,data$n.cells.x,data$n.cells.y),
      main="Habitat Mask",xlab="X",ylab="Y",col=cols1)
points(s.init[,,1],s.init[,,2],pch=16)

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.year=n.year,M=M,J=J,K1D=K1D,D.cov=D.cov,
                  n.cells=n.cells,n.cells.x=n.cells.x,
                  n.cells.y=n.cells.y,res=res,
                  x.vals=x.vals,y.vals=y.vals,
                  xlim=xlim,ylim=ylim,
                  cellArea=cellArea,n.cells=n.cells,
                  res=res)

#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 z=z.init,z.start=z.start.init,z.stop=z.stop.init,
                 s=s.init,beta0.phi=0,beta1.phi=0,D0=N.init[1]/(sum(InSS)*res^2),
                 sigma.move=sigma.move.init,D.beta1=D.beta1.init,rsf.beta=rsf.beta.init)

#data for Nimble
Nimdata <- list(y=y.nim,phi.cov=phi.cov.data,X=X.nim,dSS=dSS,cells=cells,InSS=InSS)

# set parameters to monitor
parameters <- c('N','gamma.fixed','N.recruit','N.survive','N.super','lambda.y1',
                'beta0.phi','beta1.phi','phi.cov.mu','phi.cov.sd','p0','sigma',
                'D0','D.beta1','rsf.beta','sigma.move')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)

config.nodes <- c('beta0.phi','beta1.phi','gamma.fixed',paste('phi.cov[',cov.up,']'),
               'phi.cov.mu','phi.cov.sd','p0','sigma','rsf.beta','sigma.move')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = FALSE)

###*required* sampler replacements###
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration? 
#25% of M seems reasonable, but optimal will depend on data set
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
s.nodes <- Rmodel$expandNodeNames(paste0("s"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
calcNodes <- c(N.nodes,ER.nodes,N.recruit.nodes,N.survive.nodes,s.nodes,z.nodes,pd.nodes,y.nodes)
#need to convert cells to double to use inside custom sampler
cells.double <- matrix(as.double(cells),n.cells.y,n.cells.x)
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,J=J,cells=cells.double,
                                                 dSS=dSS,res=res,n.cells=n.cells,
                                                 xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,
                                                 n.cells.x=n.cells.x,n.cells.y=n.cells.y,
                                                 z.obs=z.obs,z.super.ups=z.super.ups,
                                                 y.nodes=y.nodes,pd.nodes=pd.nodes,N.nodes=N.nodes,
                                                 z.nodes=z.nodes,ER.nodes=ER.nodes,s.nodes=s.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 y2D=y.nim2D,calcNodes=calcNodes), silent = TRUE)

#activity center sampler. There are 2 samplers here for these cases
#1) z.super=1 and z=1, sSampler1 uses Metropolis-Hastings
#2) z.super=1 and z=0, sSampler2 uses Metropolis-Hastings with proposal sd tuned separately from above
#z.super=0, do nothing, activity centers likelihood is gated and all values are set to 0
for(i in 1:M){
  for(g in 1:n.year){
    calcNodes <- Rmodel$getDependencies(paste0("s[",i,",",g,",1:2]"))
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler1',control=list(i=i,g=g,xlim=xlim,ylim=ylim,
                                                    calcNodes=calcNodes,scale=1),silent = TRUE)
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler2',control=list(i=i,g=g,xlim=xlim,ylim=ylim,
                                                    calcNodes=calcNodes,scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#optional (but recommended!) blocking 
conf$addSampler(target = c("beta0.phi","beta1.phi"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)

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

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:500),]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$truth$N.super #N.super
sigma.move
rsf.beta


rem.idx <- c(grep("N",colnames(mvSamples)),
             grep("N.recruit",colnames(mvSamples)),
             grep("N.survive",colnames(mvSamples)))
tmp <- cor(mvSamples[-c(1:500),-rem.idx])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)
round(tmp,2)
