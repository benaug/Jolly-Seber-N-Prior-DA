#function to initialize s consistent with s.move.init
initialize.s <- function(sigma.move.init=NA,z.super.init=NA,y=NA,X=NA,xlim=NA,ylim=NA){
  M <- nrow(y)
  n.year <- dim(y)[2]
  s.init <- array(0, dim=c(M, n.year, 2))  # initialize all to 0, keep 0 if z.super[i]=0
  on.inds <- which(z.super.init == 1)
  for(i in on.inds){
    dets <- which(rowSums(matrix(y[i,,],nrow=n.year))>0)
    if(length(dets)>0){
      first.det <- min(dets)
      last.det <- max(dets)
      #Set each detected year to mean trap location in that year
      for(g in dets){
        trps <- matrix(X[g,which(y[i,g,]>0),],ncol=2)
        s.init[i,g,] <- c(mean(trps[,1]),mean(trps[,2]))
      }
      #Simulate backwards from first detection year
      if(first.det > 1){
        for(g in (first.det-1):1){
          # x coordinate
          F.a <- pnorm(xlim[1], s.init[i,g+1,1], sigma.move.init)
          F.b <- pnorm(xlim[2], s.init[i,g+1,1], sigma.move.init)
          s.init[i,g,1] <- qnorm(runif(1, F.a, F.b), s.init[i,g+1,1], sigma.move.init)
          # y coordinate
          F.a <- pnorm(ylim[1], s.init[i,g+1,2], sigma.move.init)
          F.b <- pnorm(ylim[2], s.init[i,g+1,2], sigma.move.init)
          s.init[i,g,2] <- qnorm(runif(1, F.a, F.b), s.init[i,g+1,2], sigma.move.init)
        }
      }
      #Simulate forwards from last detection year
      if(last.det < n.year){
        for(g in (last.det+1):n.year){
          # x coordinate
          F.a <- pnorm(xlim[1], s.init[i,g-1,1], sigma.move.init)
          F.b <- pnorm(xlim[2], s.init[i,g-1,1], sigma.move.init)
          s.init[i,g,1] <- qnorm(runif(1, F.a, F.b), s.init[i,g-1,1], sigma.move.init)
          # y coordinate
          F.a <- pnorm(ylim[1], s.init[i,g-1,2], sigma.move.init)
          F.b <- pnorm(ylim[2], s.init[i,g-1,2], sigma.move.init)
          s.init[i,g,2] <- qnorm(runif(1, F.a, F.b), s.init[i,g-1,2], sigma.move.init)
        }
      }
      #fill in gaps between detections with linear interpolation
      if(last.det > first.det){
        for(g in first.det:(last.det-1)){
          if(!(g+1) %in% dets){
            # find surrounding detected years
            prev.det <- max(dets[dets <= g])
            next.det <- min(dets[dets > g])
            # linearly interpolate
            frac <- (g+1 - prev.det) / (next.det - prev.det)
            s.init[i,g+1,1] <- s.init[i,prev.det,1] + frac*(s.init[i,next.det,1] - s.init[i,prev.det,1])
            s.init[i,g+1,2] <- s.init[i,prev.det,2] + frac*(s.init[i,next.det,2] - s.init[i,prev.det,2])
          }
        }
      }
    }else{
      #If you initialize z.super=0 for all undetected individuals, this is not used
      #Undetected z.super=1 individual - simulate full trajectory from random start
      s.init[i,1,] <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
      for(g in 2:n.year){
        # x coordinate
        F.a <- pnorm(xlim[1], s.init[i,g-1,1], sigma.move.init)
        F.b <- pnorm(xlim[2], s.init[i,g-1,1], sigma.move.init)
        s.init[i,g,1] <- qnorm(runif(1, F.a, F.b), s.init[i,g-1,1], sigma.move.init)
        # y coordinate
        F.a <- pnorm(ylim[1], s.init[i,g-1,2], sigma.move.init)
        F.b <- pnorm(ylim[2], s.init[i,g-1,2], sigma.move.init)
        s.init[i,g,2] <- qnorm(runif(1, F.a, F.b), s.init[i,g-1,2], sigma.move.init)
      }
    }
  }
  return(s.init)
}

#truncated normals for both x and y dimensions, used for movement
dTruncNorm <- nimbleFunction(
  run = function(x = double(1), s.prev = double(1), xlim = double(1), ylim = double(1),
                 sigma.move = double(0), z.super = double(0), log = integer(0)) {
    returnType(double(0))
    if(z.super==1){
      logProb.x <- dnorm(x[1], s.prev[1], sigma.move, log = TRUE) -
        log(pnorm(xlim[2], s.prev[1], sigma.move) - pnorm(xlim[1], s.prev[1], sigma.move))
      logProb.y <- dnorm(x[2], s.prev[2], sigma.move, log = TRUE) -
        log(pnorm(ylim[2], s.prev[2], sigma.move) - pnorm(ylim[1], s.prev[2], sigma.move))
      logProb <- logProb.x + logProb.y
    }else{
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    } 
  }
)

#function to generate Normal RVs truncated by state space
rTruncNorm <- nimbleFunction(
  run = function(n = integer(0), s.prev = double(1), xlim = double(1), ylim = double(1),
                 sigma.move = double(0),z.super = double(0)) {
    if(n!=1){
      print("rTruncNorm only accepts n=1")
    }
    returnType(double(1))
    if(z.super==1){
      #x
      F.a <- pnorm(xlim[1],s.prev[1],sigma.move)
      F.b <- pnorm(xlim[2],s.prev[1],sigma.move)
      u <- runif(n,F.a,F.b)
      s.prop.x <- qnorm(u,s.prev[1],sigma.move)
      #y
      F.a <- pnorm(ylim[1],s.prev[2],sigma.move)
      F.b <- pnorm(ylim[2],s.prev[2],sigma.move)
      u <- runif(n,F.a,F.b)
      s.prop.y <- qnorm(u,s.prev[2],sigma.move)
      s.prop <- c(s.prop.x,s.prop.y)
    }else{
      s.prop <- c(0,0)
    }
    return(s.prop)
  }
)

dunif2D <- nimbleFunction(
  run = function(x = double(1), xlim = double(1), ylim = double(1), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){#skip calculation if not is superpop
      logProb <- -log(xlim[2]-xlim[1]) - log(ylim[2]-ylim[1])
    }else{
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    } 
  }
)

#make dummy random vector generator to make nimble happy
runif2D <- nimbleFunction(
  run = function(n = integer(0), xlim = double(1), ylim = double(1), z.super = double(0)) {
    returnType(double(1))
    x <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
    return(x)
  }
)


#this is used to restrict likelihood evaluation to only the years relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0),
                 z.super = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.year <- length(phi)+1
      #extract first and last survival event years
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.year){ #if surv.start beyond last year, no survival events, logProb=0
        if(surv.stop > n.year){ #but can't survive past n.year
          surv.stop <- n.year 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event years
          logProb <- logProb + dbinom(x[g], size = 1, p = phi[g-1], log = TRUE)
        }
      }
    }
    return(logProb)
  }
)

#make dummy random vector generator to make nimble happy
rSurvival <- nimbleFunction(
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0)){
    returnType(double(1))
    n.year <- length(phi)
    return(rep(0,n.year))
  }
)

GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this year
    }
    if(z==1){ #otherwise calculate
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K = double(1), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this year
      return(0)
    }else{
      logProb <- sum(dbinom(x, size = K, p = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0), pd = double(1), K = double(1), z = double(0), z.super = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    J <- control$J
    y2D <- control$y2D
    xlim <- control$xlim
    ylim <- control$ylim
    z.super.ups <- control$z.super.ups
    n.year <- control$n.year
    z.obs <- control$z.obs
    z.nodes <- control$z.nodes
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    s.nodes <- control$s.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    entry.counts.curr <- rep(0, n.year+1)
    for(g in 1:n.year){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1)
    }
    entry.counts.curr[n.year + 1] <- sum(model$z.super==0)
    
    #1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,1]==0){ #for detected guys, skip if observed 1st year
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y2D[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y and pd nodes
        for(g in 1:first.det){ #must be recruited in year with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop

          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx]) #update pd nodes when a z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)

        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original

        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx]) #update pd nodes
          #update these logProbs
          model$calculate(y.nodes[i.idx])
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(z.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["pd",1][i,g,j] <<- model[["pd"]][i,g,j]
            }
          }
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["pd"]][i,g,j] <<- mvSaved["pd",1][i,g,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
        }
      }
    }

    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,n.year]==0){ #for detected guys, skip if observed in final year
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y2D[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y and pd nodes
        for(g in (last.det):n.year){ #can't die on or before year of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.year)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx]) #update pd nodes when z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.year)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx]) #update pd nodes when z changes
          #update these logProbs
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["pd",1][i,g,j] <<- model[["pd"]][i,g,j]
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["pd"]][i,g,j] <<- mvSaved["pd",1][i,g,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #2) undetected guy update. Only if in the superpopulation.
    # Metropolis-Hastings, Propose z vectors from priors
    #entry counts current after z.start update
    for(i in 1:M){
      if(z.obs[i]==0&model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx <- seq(i,M*n.year,M) #used to reference correct y and pd nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #get initial logProbs
        lp.initial.entry <- model$getLogProb(N.nodes[1])
        lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
        lp.initial.y <- model$getLogProb(y.nodes[i.idx])
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))

        #track proposal probs - survival is symmetric, but not recruitment and detection
        log.prop.for <- log.prop.back <- 0

        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.year)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])

        #simulate survival
        if(z.start.prop < n.year){#if you don't recruit in final year
          for(g in (z.start.prop+1):n.year){
            z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
            log.prop.for <- log.prop.for + dbinom(z.prop[g],1,model$phi[i,g-1]*z.prop[g-1],log=TRUE)
          }
        }
        z.on.prop <- which(z.prop==1)
        z.stop.prop <- max(z.on.prop)
        model$z[i,] <<- z.prop
        model$z.start[i] <<- z.start.prop
        model$z.stop[i] <<- z.stop.prop

        #update N, N.recruit, N.survive only if individual is in superpopulation
        #1) Update N
        model$N <<- model$N - z.curr + z.prop
        #2) Update N.recruit
        if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
          model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
        }
        if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
          model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g

        model$calculate(ER.nodes) #update ER when N updated
        model$calculate(pd.nodes[i.idx]) #update pd nodes when z changes
        #get proposed logProbs
        lp.proposed.entry <- model$calculate(N.nodes[1])
        lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
        lp.proposed.y <- model$calculate(y.nodes[i.idx])
        lp.proposed.surv <- model$calculate(z.nodes[i])

        # Full multinomial coefficient prior for proposed configuration
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
        entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
        log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))

        #get backwards proposal probs
        recruit.probs.back <- c(model$lambda.y1,model$ER)
        recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
        log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
        if(z.start.curr < n.year){#if you don't recruit in final year
          for(g in (z.start.curr+1):n.year){
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[i,g-1]*z.curr[g-1],log=TRUE)
          }
        }
        lp.initial.total <- lp.initial.entry + lp.initial.y + lp.initial.surv + log.prior.curr
        lp.proposed.total <- lp.proposed.entry + lp.proposed.y + lp.proposed.surv + log.prior.prop

        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        # log_MH_ratio <- (lp.proposed) - (lp.initial)
        accept <- decide(log_MH_ratio)

        if(accept) {
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["pd",1][i,g,j] <<- model[["pd"]][i,g,j]
            }
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["pd"]][i,g,j] <<- mvSaved["pd",1][i,g,j]
            }
          }
          #set these logProbs back
          model$calculate(N.recruit.nodes)
          model$calculate(N.nodes[1])
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #3) update z.super: Metropolis-Hastings
    #entry counts current coming out of undetected ind update
    #including s and z likelihoods and proposal probs for clarity, they both cancel in MH ratio
    #z survival proposals cancel with survival likelihood, but recruit prob included in z proposal probs.
    #could just use recruitment proposal prob, but including survival parts for clarity.
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z.super==1)
        non.init <- length(z.on)
        pick <- rcat(1,rep(1/non.init,non.init))
        pick <- z.on[pick]
        if(z.obs[pick]==1){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          s.curr <- model$s[pick,,]
          #p select off guy
          log.p.select.for <- log(1/non.init)
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx <- seq(pick,M*n.year,M) #used to reference correct y nodes
          
          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
          lp.initial.s <- model$getLogProb(s.nodes[pick.idx])
          lp.initial.surv <- model$getLogProb(z.nodes[pick])
          
          # propose new N.super/z.super/z.start/z.stop
          model$N.super <<-  model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.year)
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          
          #set s to all 0's
          log.prop.for.s <- 0
          for(g in 1:n.year){
            model$s[pick,g,1:2] <<- c(0,0)
          }
          model$calculate(pd.nodes[pick.idx]) #update pd nodes when z changes
          
          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1, model$ER)
          recruit.probs.back <- recruit.probs.back / sum(recruit.probs.back)
          log.prop.back.z <- log(recruit.probs.back[z.start.curr])
          if(z.start.curr < n.year){
            for(g in (z.start.curr+1):n.year){
              log.prop.back.z <- log.prop.back.z + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
            }
          }
          log.prop.back.s <- -log(xlim[2]-xlim[1]) - log(ylim[2]-ylim[1])
          #propose subsequent years from movement prior (truncated Normal here)
          for(g in 2:n.year){
            log.prop.back.s <- log.prop.back.s + dTruncNorm(s.curr[g,1:2],s.prev=s.curr[g-1,1:2],
                                                          sigma.move=model$sigma.move[1],
                                                          xlim=xlim, ylim=ylim,z.super=1,log=TRUE)
          }
          
          #get proposed logProbs for N, N.recruit, and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0
          lp.proposed.s <- model$calculate(s.nodes[pick.idx]) #will always be 0
          lp.proposed.surv <- model$calculate(z.nodes[pick]) #will always be 0
          
          lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit + lp.initial.surv + lp.initial.s
          lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.s
          
          #backwards prior and select probs
          #move from class z.start.curr in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.year + 1] <- entry.counts.prop[n.year + 1] + 1
          
          #p select on guy
          noff.back <- sum(model$z.super == 0)
          log.p.select.back <- log(1/noff.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.for.z <- log.prop.for.s <- 0
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back.z + log.prop.back.s) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for.z + log.prop.for.s)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            mvSaved["s",1][pick,1:n.year,1:2] <<- model[["s"]][pick,1:n.year,1:2]
            for(g in 1:n.year){
              for(j in 1:J[g]){
                mvSaved["pd",1][pick,g,j] <<- model[["pd"]][pick,g,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            model[["s"]][pick,1:n.year,1:2] <<- mvSaved["s",1][pick,1:n.year,1:2]
            for(g in 1:n.year){
              for(j in 1:J[g]){
                model[["pd"]][pick,g,j] <<- mvSaved["pd",1][pick,g,j]
              }
            }
            #set these logProbs back
            model$calculate(s.nodes[pick.idx])
            model$calculate(z.nodes[pick])
            model$calculate(y.nodes[pick.idx])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          z.off <- which(model$z.super==0)
          noff.init <- length(z.off)
          pick <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx <- seq(pick,M*n.year,M)

          #p select off guy
          log.p.select.for <- log(1/noff.init)
          
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          
          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx]) #will always be 0
          lp.initial.s <- model$getLogProb(s.nodes[pick.idx]) #will always be 0
          lp.initial.surv <- model$getLogProb(z.nodes[pick]) #will always be 0
          
          # Propose new z.start for the new on individual
          recruit.probs.for <- c(model$lambda.y1, model$ER)
          recruit.probs.for <- recruit.probs.for / sum(recruit.probs.for)
          z.start.prop <- rcat(1, recruit.probs.for)  # propose entry cohort
          log.prop.for.z <- log(recruit.probs.for[z.start.prop])
          model$z.start[pick] <<- z.start.prop
          
          # Simulate survival path
          model$z[pick,] <<- 0 # initialize to 0
          model$z[pick, z.start.prop] <<- 1
          if(z.start.prop < n.year){
            for(g in (z.start.prop+1):n.year){
              model$z[pick, g] <<- rbinom(1, 1, model$phi[pick, g-1] * model$z[pick, g-1])
              log.prop.for.z <- log.prop.for.z + dbinom(model$z[pick, g], 1, model$phi[pick, g-1] * model$z[pick, g-1], log=TRUE)
            }
          }
          # Update z.stop
          z.on.prop <- which(model$z[pick,] == 1)
          z.stop.prop <- max(z.on.prop)
          model$z.stop[pick] <<- z.stop.prop
          
          #propose new N/z
          model$N.super <<-  model$N.super + 1
          model$z.super[pick] <<- 1
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N + model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year] - model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          #simulate new s trajectory and record forward proposal probs (prior and likelihood cancel, but including both below for clarity)
          #propose year 1 from uniform prior
          model$s[pick,1,1:2] <<- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
          log.prop.for.s <- -log(xlim[2]-xlim[1]) - log(ylim[2]-ylim[1])
          #propose subsequent years from movement prior (truncated Normal here)
          for(g in 2:n.year){
            model$s[pick,g,1:2] <<- rTruncNorm(1,s.prev=model$s[pick,g-1,1:2],
                                               sigma.move=model$sigma.move[1],
                                               xlim=xlim, ylim=ylim,z.super=1)
            log.prop.for.s <- log.prop.for.s + dTruncNorm(model$s[pick,g,1:2],s.prev=model$s[pick,g-1,1:2],
                                                          sigma.move=model$sigma.move[1],
                                                          xlim=xlim, ylim=ylim,z.super=1,log=TRUE)
          }
          model$calculate(pd.nodes[pick.idx]) #update pd nodes when z and s change
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx])
          lp.proposed.s <- model$calculate(s.nodes[pick.idx])
          lp.proposed.surv <- model$calculate(z.nodes[pick])
          
          lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit + lp.initial.surv + lp.initial.s
          lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.s
          
          #backwards prior and select probs
          #move from class g in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.prop[n.year + 1] <- entry.counts.prop[n.year + 1] - 1
          
          #p select on guy
          non.back <- sum(model$z.super == 1)
          log.p.select.back <- log(1/non.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.back.z <- log.prop.back.s <-  0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back.z + log.prop.back.s) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for.z + log.prop.for.s)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            mvSaved["s",1][pick,1:n.year,1:2] <<- model[["s"]][pick,1:n.year,1:2]
            for(g in 1:n.year){
              for(j in 1:J[g]){
                mvSaved["pd",1][pick,g,j] <<- model[["pd"]][pick,g,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            model[["s"]][pick,1:n.year,1:2] <<- mvSaved["s",1][pick,1:n.year,1:2]
            for(g in 1:n.year){
              for(j in 1:J[g]){
                model[["pd"]][pick,g,j] <<- mvSaved["pd",1][pick,g,j]
              }
            }
            #set these logProbs back
            model$calculate(y.nodes[pick.idx])
            model$calculate(z.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            model$calculate(s.nodes[pick.idx])
          }
        }
      }
    }
    
    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)