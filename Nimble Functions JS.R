#this is used to restrict likelihood evaluation to only the primary occasions relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0),
                 z.super = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.primary <- length(phi)+1
      #extract first and last survival event primary occasions
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.primary){ #if surv.start beyond last primary occasion, no survival events, logProb=0
        if(surv.stop > n.primary){ #but can't survive past n.primary
          surv.stop <- n.primary 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event primary occasions
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
    n.primary <- length(phi)+1
    return(rep(0,n.primary))
  }
)

#custom observation model distribution that makes custom updates easier
dbinomial2 <- nimbleFunction(
  run = function(x = double(0), p = double(0), K = double(0), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
      return(0)
    }else{
      logProb <- dbinom(x, size = K, p = p, log = TRUE)
      return(logProb)
    }
  }
)

rbinomial2 <- nimbleFunction(
  run = function(n = integer(0), p = double(0), K = double(0), z = double(0), z.super = double(0)){
    returnType(double(0))
    return(0)
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    z.super.ups <- control$z.super.ups
    n.primary <- control$n.primary
    z.obs <- control$z.obs
    z.nodes <- control$z.nodes
    y.nodes <- control$y.nodes
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    entry.counts.curr <- rep(0,n.primary+1)
    for(i in 1:M){
      if(model$z.super[i]==1){
        entry.counts.curr[model$z.start[i]] <- entry.counts.curr[model$z.start[i]]+1
      }else{
        entry.counts.curr[n.primary+1] <- entry.counts.curr[n.primary+1]+1
      }
    }
    
    #1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,1]==0){ #for detected guys, skip if observed 1st primary occasion
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(model$y[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.primary)
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
        #only y nodes before first detection can change across z.start candidates
        y.idx <- i.idx[1:(first.det-1)]
        #remove focal individual from entry counts once. The candidate-specific part of the
        #multinomial coefficient is then just log(entry.counts.minus[g]+1)
        entry.counts.minus <- entry.counts.curr
        entry.counts.minus[z.start.curr] <- entry.counts.minus[z.start.curr] - 1
        #all z.start > 1 candidates have the same N[1], so only calculate that logProb once
        lp.N1.not1 <- 0
        # Here, we are looping over all valid recruit dates and storing the logProb for each
        for(g in 1:first.det){ #must be recruited in primary occasion with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #only ER nodes before first detection can change across z.start candidates
          #model$calculate(ER.nodes) #update ER when N updated
          model$calculate(ER.nodes[1:(first.det-1)])
          #get these logProbs
          #lp.N1 <- model$calculate(N.nodes[1])
          #there are only two possible N[1] values: z.start=1 and z.start>1
          if(g==1){
            lp.N1 <- model$calculate(N.nodes[1])
          }else{
            if(g==2){
              lp.N1.not1 <- model$calculate(N.nodes[1])
            }
            lp.N1 <- lp.N1.not1
          }
          #only recruitment likelihoods before first detection can change
          #lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.N.recruit <- model$calculate(N.recruit.nodes[1:(first.det-1)])
          #only observation likelihoods before first detection can change
          #lp.y <- model$calculate(y.nodes[i.idx])
          lp.y <- model$calculate(y.nodes[y.idx])
          lp.surv <- model$calculate(z.nodes[i])
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          #entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          #entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          #entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          #lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          #after removing this individual, all multinomial coefficient terms common across
          #candidates cancel, leaving log(n.g+1) for candidate entry class g
          lp.prior <- log(entry.counts.minus[g]+1)
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #model$calculate(ER.nodes) #update ER when N updated
          model$calculate(ER.nodes[1:(first.det-1)])
          #update these logProbs
          #model$calculate(y.nodes[i.idx])
          model$calculate(y.nodes[y.idx])
          model$calculate(N.nodes[1])
          #model$calculate(N.recruit.nodes)
          model$calculate(N.recruit.nodes[1:(first.det-1)])
          model$calculate(z.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
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
          #set these logProbs back
          model$calculate(N.nodes[1])
          #model$calculate(N.recruit.nodes)
          model$calculate(N.recruit.nodes[1:(first.det-1)])
          #model$calculate(y.nodes[i.idx])
          model$calculate(y.nodes[y.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    
    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,n.primary]==0){ #for detected guys, skip if observed in final primary occasion
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(model$y[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.primary)
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
        #only y nodes after last detection can change across z.stop candidates
        y.idx <- i.idx[(last.det+1):n.primary]
        #Here, we are looping over all valid z.stops and storing the logProb for each
        #entry prior does not change in z.stop update
        for(g in (last.det):n.primary){ #can't die on or before primary occasion of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          #only ER nodes after last detection can change across z.stop candidates
          #model$calculate(ER.nodes) #update ER when N updated
          if(last.det < n.primary-1){
            model$calculate(ER.nodes[(last.det+1):(n.primary-1)])
          }
          #get these logProbs
          #lp.N1 <- model$calculate(N.nodes[1])
          #N[1] cannot change in a z.stop update, so this is constant across candidates and cancels
          #lp.N.recruit <- model$calculate(N.recruit.nodes)
          if(last.det < n.primary-1){
            lp.N.recruit <- model$calculate(N.recruit.nodes[(last.det+1):(n.primary-1)])
          }else{
            lp.N.recruit <- 0
          }
          #lp.y <- model$calculate(y.nodes[i.idx])
          lp.y <- model$calculate(y.nodes[y.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #no prior term, z.stop update does not change it
          #lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv
          lp.stop[g] <- lp.N.recruit + lp.y + lp.surv
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.primary)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #model$calculate(ER.nodes) #update ER when N updated
          if(last.det < n.primary-1){
            model$calculate(ER.nodes[(last.det+1):(n.primary-1)])
          }
          #update these logProbs
          #model$calculate(N.nodes[1]) #N[1] does not change in a z.stop update
          #model$calculate(N.recruit.nodes)
          if(last.det < n.primary-1){
            model$calculate(N.recruit.nodes[(last.det+1):(n.primary-1)])
          }
          #model$calculate(y.nodes[i.idx])
          model$calculate(y.nodes[y.idx])
          model$calculate(z.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          #set these logProbs back
          #model$calculate(N.nodes[1]) #N[1] does not change in a z.stop update
          #model$calculate(N.recruit.nodes)
          if(last.det < n.primary-1){
            model$calculate(N.recruit.nodes[(last.det+1):(n.primary-1)])
          }
          #model$calculate(y.nodes[i.idx])
          model$calculate(y.nodes[y.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #2) undetected guy update. Only if in the superpopulation. 
    #Metropolis-Hastings, Propose z vectors from priors.
    #can try Gibbs below instead. not sure which is more efficient in which conditions
    #MH is cheap, but wont mix as well, Gibbs is expensive, but mixes better.
    #entry counts current after z.start update
    for(i in 1:M){
      if(z.obs[i]==0&model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        
        #track proposal probs
        #survival proposal probabilities cancel exactly with the survival likelihood because
        #the survival history is proposed from the same survival model used in the target
        log.prop.for <- log.prop.back <- 0
        
        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.primary)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])
        
        #simulate survival
        #once the individual dies, remaining z's are already 0 so no more simulation is needed
        z.stop.prop <- z.start.prop
        if(z.start.prop < n.primary){#if you don't recruit in final primary occasion
          for(g in (z.start.prop+1):n.primary){
            if(z.prop[g-1]==1){
              z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
              if(z.prop[g]==1){
                z.stop.prop <- g
              }
            }
          }
        }
        
        #if the proposed history is the current history, there is nothing to calculate or update
        if(z.start.prop!=z.start.curr|z.stop.prop!=z.stop.curr){
          #get initial logProbs
          lp.initial.entry <- model$getLogProb(N.nodes[1])
          lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
          # lp.initial.y <- model$getLogProb(y.nodes[i.idx])
          #restrict y nodes to those where z changed
          z.changed <- which(z.prop!=z.curr)
          y.idx.changed <- i.idx[z.changed]
          lp.initial.y <- model$getLogProb(y.nodes[y.idx.changed])
          #lp.initial.surv <- model$getLogProb(z.nodes[i]) #cancels exactly with backwards survival proposal probability
          #log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1))) #full multinomial coefficient calculation replaced by exact ratio below
          
          model$z[i,] <<- z.prop
          model$z.start[i] <<- z.start.prop
          model$z.stop[i] <<- z.stop.prop
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr + z.prop
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          
          model$calculate(ER.nodes) #update ER when N updated
          #get proposed logProbs
          lp.proposed.entry <- model$calculate(N.nodes[1])
          lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
          # lp.proposed.y <- model$calculate(y.nodes[i.idx])
          lp.proposed.y <- model$calculate(y.nodes[y.idx.changed])
          #lp.proposed.surv <- model$calculate(z.nodes[i]) #cancels exactly with forwards survival proposal probability
          
          # Full multinomial coefficient prior for proposed configuration
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          #log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          #In log.prior.prop-log.prior.curr, -lgamma(M+1) and all unchanged entry-class
          #lgamma terms cancel. If entry changes from class a to b, the remaining ratio is
          #log(n.b+1)-log(n.a). If entry does not change, the ratio is 0.
          if(z.start.prop!=z.start.curr){
            log.prior.ratio <- log(entry.counts.curr[z.start.prop]+1)-log(entry.counts.curr[z.start.curr])
          }else{
            log.prior.ratio <- 0
          }
          
          #get backwards proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
          #survival proposal probabilities are not calculated because they cancel exactly
          #with the survival likelihood ratio in the MH ratio
          
          lp.initial.total <- lp.initial.entry + lp.initial.y
          lp.proposed.total <- lp.proposed.entry + lp.proposed.y
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.prior.ratio + log.prop.back) - (lp.initial.total + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #update survival logProb once for accepted history; it was not needed when evaluating MH ratio
            model$calculate(z.nodes[i])
            mvSaved["z.start",1][i] <<- model[["z.start"]][i]
            mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
            mvSaved["z",1][i,] <<- model[["z"]][i,]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["ER",1] <<- model[["ER"]]
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][i] <<- mvSaved["z.start",1][i]
            model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
            model[["z"]][i,] <<- mvSaved["z",1][i,]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["ER"]] <<- mvSaved["ER",1]
            #set these logProbs back
            # model$calculate(y.nodes[i.idx])
            model$calculate(y.nodes[y.idx.changed])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            #model$calculate(z.nodes[i]) #not needed because survival logProb was never recalculated for the proposal
          }
        }
      }
    }
    # #2) undetected guy update. Only if in the superpopulation. 
    # #Gibbs updates for z.start and z.stop
    #was slower than MH in 1 comparison
    # #entry counts current after detected z.start update
    # for(i in 1:M){
    #   if(z.obs[i]==0&model$z.super[i]==1){
    #     #2a) z.start update
    #     z.curr <- model$z[i,]
    #     z.start.curr <- model$z.start[i]
    #     z.stop.curr <- model$z.stop[i]
    #     N.curr <- model$N
    #     N.recruit.curr <- model$N.recruit
    #     lp.start <- rep(-Inf,n.primary)
    #     i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
    #     #Here, we are looping over all valid recruit dates and storing the logProb for each
    #     for(g in 1:z.stop.curr){
    #       z.start.prop <- g
    #       model$z.start[i] <<- z.start.prop
    #       z.prop <- rep(0,n.primary)
    #       z.prop[g:z.stop.curr] <- 1
    #       model$z[i,] <<- z.prop
    #       
    #       #update N, N.recruit, N.survive
    #       #1) Update N
    #       model$N <<- N.curr - z.curr + z.prop
    #       #2) Update N.recruit
    #       model$N.recruit <<- N.recruit.curr #set back to original first
    #       if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
    #         model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
    #       }
    #       if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
    #         model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
    #       }
    #       #3) Update N.survive
    #       model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
    #       model$calculate(ER.nodes) #update ER when N updated
    #       
    #       #get these logProbs
    #       lp.N1 <- model$calculate(N.nodes[1])
    #       lp.N.recruit <- model$calculate(N.recruit.nodes)
    #       lp.y <- model$calculate(y.nodes[i.idx])
    #       lp.surv <- model$calculate(z.nodes[i])
    #       # Add the full multinomial coefficient prior log-prob for this proposed configuration
    #       entry.counts.prop <- entry.counts.curr
    #       entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
    #       entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
    #       lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
    #       lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.prior
    #     }
    #     maxlp <- max(lp.start)
    #     prop.probs <- exp(lp.start-maxlp)
    #     prop.probs <- prop.probs/sum(prop.probs)
    #     z.start.prop <- rcat(1,prop.probs)
    #     model$z.start[i] <<- z.start.curr
    #     
    #     if(model$z.start[i]!=z.start.prop){
    #       model$z.start[i] <<- z.start.prop
    #       z.prop <- rep(0,n.primary)
    #       z.prop[z.start.prop:z.stop.curr] <- 1
    #       model$z[i,] <<- z.prop
    #       model$N <<- N.curr - z.curr + z.prop
    #       model$N.recruit <<- N.recruit.curr #set back to original first
    #       if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
    #         model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
    #       }
    #       if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
    #         model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
    #       }
    #       model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
    #       model$calculate(ER.nodes) #update ER when N updated
    #       #update these logProbs
    #       model$calculate(y.nodes[i.idx])
    #       model$calculate(N.nodes[1])
    #       model$calculate(N.recruit.nodes)
    #       model$calculate(z.nodes[i])
    #       mvSaved["z.start",1][i] <<- model[["z.start"]][i]
    #       mvSaved["z",1][i,] <<- model[["z"]][i,]
    #       mvSaved["N",1] <<- model[["N"]]
    #       mvSaved["N.survive",1] <<- model[["N.survive"]]
    #       mvSaved["N.recruit",1] <<- model[["N.recruit"]]
    #       mvSaved["ER",1] <<- model[["ER"]]
    #       #recompute entry counts
    #       entry.counts.prop <- entry.counts.curr
    #       entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
    #       entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
    #       entry.counts.curr <- entry.counts.prop
    #     }else{
    #       model[["z.start"]][i] <<- mvSaved["z.start",1][i]
    #       model[["z"]][i,] <<- mvSaved["z",1][i,]
    #       model[["N"]] <<- mvSaved["N",1]
    #       model[["N.survive"]] <<- mvSaved["N.survive",1]
    #       model[["N.recruit"]] <<- mvSaved["N.recruit",1]
    #       model[["ER"]] <<- mvSaved["ER",1]
    #       #set these logProbs back
    #       model$calculate(N.nodes[1])
    #       model$calculate(N.recruit.nodes)
    #       model$calculate(y.nodes[i.idx])
    #       model$calculate(z.nodes[i])
    #     }
    #     
    #     #2b) z.stop update
    #     z.curr <- model$z[i,]
    #     z.start.curr <- model$z.start[i]
    #     z.stop.curr <- model$z.stop[i]
    #     N.curr <- model$N
    #     lp.stop <- rep(-Inf,n.primary)
    #     #Here, we are looping over all valid z.stops and storing the logProb for each
    #     for(g in z.start.curr:n.primary){
    #       model$z.stop[i] <<- g
    #       z.prop <- rep(0,n.primary)
    #       z.prop[z.start.curr:g] <- 1
    #       model$z[i,] <<- z.prop
    #       #update N, number of recruits does not change going backwards
    #       model$N <<- N.curr - z.curr + z.prop
    #       model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
    #       model$calculate(ER.nodes) #update ER when N updated
    #       #get these logProbs
    #       lp.N1 <- model$calculate(N.nodes[1])
    #       lp.N.recruit <- model$calculate(N.recruit.nodes)
    #       lp.y <- model$calculate(y.nodes[i.idx])
    #       lp.surv <- model$calculate(z.nodes[i])
    #       #no prior term, z.stop update does not change it
    #       lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv
    #     }
    #     maxlp <- max(lp.stop)
    #     prop.probs <- exp(lp.stop-maxlp)
    #     prop.probs <- prop.probs/sum(prop.probs)
    #     z.stop.prop <- rcat(1,prop.probs)
    #     model$z.stop[i] <<- z.stop.curr
    #     
    #     if(model$z.stop[i]!=z.stop.prop){
    #       model$z.stop[i] <<- z.stop.prop
    #       z.prop <- rep(0,n.primary)
    #       z.prop[z.start.curr:z.stop.prop] <- 1
    #       model$z[i,] <<- z.prop
    #       model$N <<- N.curr - z.curr + z.prop
    #       model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
    #       model$calculate(ER.nodes) #update ER when N updated
    #       #update these logProbs
    #       model$calculate(N.nodes[1])
    #       model$calculate(N.recruit.nodes)
    #       model$calculate(y.nodes[i.idx])
    #       model$calculate(z.nodes[i])
    #       mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
    #       mvSaved["z",1][i,] <<- model[["z"]][i,]
    #       mvSaved["N",1] <<- model[["N"]]
    #       mvSaved["N.survive",1] <<- model[["N.survive"]]
    #       mvSaved["ER",1] <<- model[["ER"]]
    #     }else{
    #       model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
    #       model[["z"]][i,] <<- mvSaved["z",1][i,]
    #       model[["N"]] <<- mvSaved["N",1]
    #       model[["N.survive"]] <<- mvSaved["N.survive",1]
    #       model[["ER"]] <<- mvSaved["ER",1]
    #       #set these logProbs back
    #       model$calculate(N.nodes[1])
    #       model$calculate(N.recruit.nodes)
    #       model$calculate(y.nodes[i.idx])
    #       model$calculate(z.nodes[i])
    #     }
    #   }
    # }
    #3) update z.super: Metropolis-Hastings
    #entry counts current coming out of undetected ind update
    #make lists of currently on/off undetected guys once, then update after accepted proposals
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(z.obs[i]==0){
        if(model$z.super[i]==1){
          non.curr <- non.curr+1
          z.on[non.curr] <- i
        }else{
          noff.curr <- noff.curr+1
          z.off[noff.curr] <- i
        }
      }
    }
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        #find all z's currently on and undetected
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          
          #p select on guy
          log.p.select.for <- log(1/non.init)
          
          #log multinomial coefficient prior
          #log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx <- seq(pick,M*n.primary,M) #used to reference correct y nodes
          
          #get initial logProbs
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          #lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
          #only currently alive primary occasions can change when this individual is removed
          lp.initial.y <- 0
          for(g in 1:n.primary){
            if(z.curr[g]==1){
              lp.initial.y <- lp.initial.y+model$getLogProb(y.nodes[pick.idx[g]])
            }
          }
          #lp.initial.surv <- model$getLogProb(z.nodes[pick]) #survival likelihood cancels exactly with reverse survival proposal probability
          
          #propose new N.super/z.super/z.start/z.stop
          model$N.super <<- model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.primary)
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr
          
          #2) Update N.recruit
          if(z.start.curr > 1){
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit
          model$calculate(ER.nodes)
          
          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          #survival proposal probability cancels exactly with current survival likelihood
          
          #get proposed logProbs
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          #lp.proposed.y <- model$calculate(y.nodes[pick.idx])
          lp.proposed.y <- 0 #all focal observation likelihood terms are zero when z.super=0
          
          #survival target/proposal terms cancel exactly
          lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit
          lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit
          
          #multinomial coefficient ratio
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.primary+1] <- entry.counts.prop[n.primary+1] + 1
          
          #backwards selection probability
          noff.back <- noff.curr+1
          log.p.select.back <- log(1/noff.back)
          
          log.z.prior.ratio <- log(entry.counts.curr[n.primary+1]+1) - log(entry.counts.curr[z.start.curr])
          log.prop.for <- 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) - (lp.initial.total + log.p.select.for + log.prop.for)
          
          accept <- decide(log_MH_ratio)
          if(accept){
            #survival term cancelled from MH ratio, but cached logProb must
            #be synchronized to the accepted z.super=0 state
            #synchronize focal observation nodes only in primary occasions that were alive before removal
            for(g in 1:n.primary){
              if(z.curr[g]==1){
                model$calculate(y.nodes[pick.idx[g]])
              }
            }
            model$calculate(z.nodes[pick])
            
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            
            entry.counts.curr <- entry.counts.prop
            
            #move guy from on list to off list
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
            
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            
            #set these logProbs back
            #model$calculate(y.nodes[pick.idx])
            #not needed: focal y logProbs were not recalculated for the rejected removal proposal
            model$calculate(z.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
          }
        }
        
      }else{#add
        if(model$N.super[1] < M){
          #find all z's currently off and undetected
          noff.init <- noff.curr
          
          if(noff.init>0){
            pick.pos <- rcat(1,rep(1/noff.init,noff.init))
            pick <- z.off[pick.pos]
            pick.idx <- seq(pick,M*n.primary,M)
            
            #p select off guy
            log.p.select.for <- log(1/noff.init)
            
            #get initial logProbs
            lp.initial.N <- model$getLogProb(N.nodes[1])
            lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
            #lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
            lp.initial.y <- 0 #all focal observation likelihood terms are zero when z.super=0
            #lp.initial.surv <- model$getLogProb(z.nodes[pick]) #cancels with forward survival proposal
            
            #propose entry cohort
            recruit.probs.for <- c(model$lambda.y1,model$ER)
            recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
            z.start.prop <- rcat(1,recruit.probs.for)
            log.prop.for <- log(recruit.probs.for[z.start.prop])
            model$z.start[pick] <<- z.start.prop
            
            #simulate survival path
            model$z[pick,] <<- 0
            model$z[pick,z.start.prop] <<- 1
            
            if(z.start.prop < n.primary){
              for(g in (z.start.prop+1):n.primary){
                model$z[pick,g] <<- rbinom(1,1,model$phi[pick,g-1]*model$z[pick,g-1])
              }
            }
            
            #update z.stop
            z.on.prop <- which(model$z[pick,]==1)
            z.stop.prop <- max(z.on.prop)
            model$z.stop[pick] <<- z.stop.prop
            
            #propose new N.super/z.super
            model$N.super <<- model$N.super + 1
            model$z.super[pick] <<- 1
            
            #update N
            model$N <<- model$N + model$z[pick,]
            
            #update N.recruit
            if(model$z.start[pick] > 1){
              model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
            }
            
            #update N.survive
            model$N.survive <<- model$N[2:n.primary]-model$N.recruit
            model$calculate(ER.nodes)
            
            #get proposed logProbs
            lp.proposed.N <- model$calculate(N.nodes[1])
            lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
            #lp.proposed.y <- model$calculate(y.nodes[pick.idx])
            lp.proposed.y <- 0
            for(g in 1:n.primary){
              if(model$z[pick,g]==1){
                lp.proposed.y <- lp.proposed.y+model$calculate(y.nodes[pick.idx[g]])
              }
            }
            
            #survival target/proposal terms cancel exactly
            lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit
            lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit
            
            #multinomial coefficient ratio
            entry.counts.prop <- entry.counts.curr
            entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
            entry.counts.prop[n.primary+1] <- entry.counts.prop[n.primary+1] - 1
            
            #backwards selection probability
            non.back <- non.curr+1
            log.p.select.back <- log(1/non.back)
            
            log.z.prior.ratio <- log(entry.counts.curr[z.start.prop]+1) - log(entry.counts.curr[n.primary+1])
            log.prop.back <- 0
            
            #MH step
            log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) - (lp.initial.total + log.p.select.for + log.prop.for)
            
            accept <- decide(log_MH_ratio)
            if(accept){
              #survival term cancelled from MH ratio, but cached logProb must
              #be synchronized to the accepted z.super=1 survival history
              model$calculate(z.nodes[pick])
              
              mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
              mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
              mvSaved["z",1][pick,] <<- model[["z"]][pick,]
              mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
              mvSaved["N",1] <<- model[["N"]]
              mvSaved["N.survive",1] <<- model[["N.survive"]]
              mvSaved["N.recruit",1] <<- model[["N.recruit"]]
              mvSaved["N.super",1][1] <<- model[["N.super"]]
              mvSaved["ER",1] <<- model[["ER"]]
              
              entry.counts.curr <- entry.counts.prop
              
              #move guy from off list to on list
              z.off[pick.pos] <- z.off[noff.curr]
              z.off[noff.curr] <- 0
              noff.curr <- noff.curr-1
              non.curr <- non.curr+1
              z.on[non.curr] <- pick
              
            }else{
              model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
              model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
              model[["z"]][pick,] <<- mvSaved["z",1][pick,]
              model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
              model[["N"]] <<- mvSaved["N",1]
              model[["N.survive"]] <<- mvSaved["N.survive",1]
              model[["N.recruit"]] <<- mvSaved["N.recruit",1]
              model[["N.super"]] <<- mvSaved["N.super",1][1]
              model[["ER"]] <<- mvSaved["ER",1]
              
              #set these logProbs back
              #model$calculate(y.nodes[pick.idx])
              #only the proposed alive primary occasions had their focal y logProbs changed
              for(g in z.start.prop:z.stop.prop){
                model$calculate(y.nodes[pick.idx[g]])
              }
              model$calculate(z.nodes[pick])
              model$calculate(N.nodes[1])
              model$calculate(N.recruit.nodes)
            }
          }
        }
      }
    }
    
    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

truncGammaPoisSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    calcNodes <- model$getDependencies(target)
    upper <- model$getBound(target,"upper")
    tau <- control$tau
    if(target=="lambda.y1"){
      is.lambda <- TRUE
      is.fixed.gamma <- FALSE
      g <- 1
      n.recruit <- 1
    }else if(target=="gamma"){
      is.lambda <- FALSE
      is.fixed.gamma <- TRUE
      g <- 1
      recruitNodes <- grep("^N.recruit\\[",model$getNodeNames(stochOnly=TRUE),value=TRUE)
      n.recruit <- length(recruitNodes)
    }else{
      is.lambda <- FALSE
      is.fixed.gamma <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
      n.recruit <- 1
    }
  },
  run = function(){
    if(is.lambda){
      count <- model$N[1]
      rate <- 1
    }else if(is.fixed.gamma){
      count <- 0
      rate <- 0
      for(j in 1:n.recruit){
        count <- count+model$N.recruit[j]
        rate <- rate+model$N[j]*tau[j]
      }
    }else{
      count <- model$N.recruit[g]
      rate <- model$N[g]*tau[g]
    }
    if(rate>0){
      shape <- count+1
      p.upper <- pgamma(upper,shape=shape,rate=rate)
      model[[target]] <<- qgamma(runif(1,0,p.upper),shape=shape,rate=rate)
    }else{
      model[[target]] <<- runif(1,0,upper)
    }
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)

pGibbsSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    calcNodes <- model$getDependencies(target)
    M <- control$M
    n.primary <- control$n.primary
    K <- control$K
    #get Beta prior parameters directly from model
    shape1 <- model$getParam(target,"shape1")
    shape2 <- model$getParam(target,"shape2")
    if(target=="p"){
      is.fixed.p <- TRUE
      g <- 1
    }else{
      is.fixed.p <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
    }
  },
  run = function(){
    success <- 0
    failure <- 0
    if(is.fixed.p){
      for(j in 1:n.primary){
        for(i in 1:M){
          if(model$z.super[i]==1&model$z[i,j]==1){
            success <- success+model$y[i,j]
            failure <- failure+K[j]-model$y[i,j]
          }
        }
      }
    }else{
      for(i in 1:M){
        if(model$z.super[i]==1&model$z[i,g]==1){
          success <- success+model$y[i,g]
          failure <- failure+K[g]-model$y[i,g]
        }
      }
    }
    model[[target]] <<- rbeta(1,shape1+success,shape2+failure)
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)