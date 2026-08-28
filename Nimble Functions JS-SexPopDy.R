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
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0)) {
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
  run = function(n = integer(0), p = double(0), K = double(0), z = double(0), z.super = double(0)) {
    returnType(double(0))
    return(0)
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.super.ups <- control$z.super.ups
    n.primary <- control$n.primary
    z.obs <- control$z.obs
    z.nodes <- control$z.nodes
    y.nodes <- control$y.nodes
    phi.nodes <- control$phi.nodes
    N.M.nodes <- control$N.M.nodes
    N.F.nodes <- control$N.F.nodes
    ER.M.nodes <- control$ER.M.nodes
    ER.F.nodes <- control$ER.F.nodes
    N.recruit.M.nodes <- control$N.recruit.M.nodes
    N.recruit.F.nodes <- control$N.recruit.F.nodes
    sex.up <- control$sex.up
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    #male slots, then female, then z.super=0 slot
    entry.counts.curr <- rep(0, 2*n.primary+1)
    for(g in 1:n.primary){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1 & model$sex==0)
      entry.counts.curr[g + n.primary] <- sum(model$z.start==g & model$z.super==1 & model$sex==1)
    }
    entry.counts.curr[2*n.primary + 1] <- sum(model$z.super==0)
    
    # 1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,1]==0){ #for detected guys, skip if observed 1st primary occasion
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        if(model$sex[i]==0){
          N.M.curr <- model$N.M
          N.recruit.M.curr <- model$N.recruit.M
        }else{
          N.F.curr <- model$N.F
          N.recruit.F.curr <- model$N.recruit.F
        }
        dets <- which(model$y[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.primary)
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
        #only y nodes before first detection can change across z.start candidates
        y.idx <- i.idx[1:(first.det-1)]
        #remove focal individual from its sex-specific entry class once. The candidate-specific
        #part of the multinomial coefficient is then just log(entry.counts.minus[cohort]+1)
        cohort.curr <- z.start.curr + model$sex[i]*n.primary
        entry.counts.minus <- entry.counts.curr
        entry.counts.minus[cohort.curr] <- entry.counts.minus[cohort.curr] - 1
        #all z.start > 1 candidates have the same focal-sex N[1], so only calculate that logProb once
        lp.N1.not1 <- 0
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
          
          #now repeat for sex
          if(model$sex[i]==0){ #male
            #1) Update N
            model$N.M <<- N.M.curr - z.curr + z.prop
            #2) Update N.recruit
            model$N.recruit.M <<- N.recruit.M.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          }else{ #female
            #1) Update N
            model$N.F <<- N.F.curr - z.curr + z.prop
            #2) Update N.recruit
            model$N.recruit.F <<- N.recruit.F.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          }
          #Updating total N changes both male and female expected recruitment.
          #Only intervals before first detection can change across z.start candidates.
          model$calculate(ER.M.nodes[1:(first.det-1)])
          model$calculate(ER.F.nodes[1:(first.det-1)])
          #there are only two possible focal-sex N[1] values: z.start=1 and z.start>1
          if(g==1){
            if(model$sex[i]==0){
              lp.N1 <- model$calculate(N.M.nodes[1])
            }else{
              lp.N1 <- model$calculate(N.F.nodes[1])
            }
          }else{
            if(g==2){
              if(model$sex[i]==0){
                lp.N1.not1 <- model$calculate(N.M.nodes[1])
              }else{
                lp.N1.not1 <- model$calculate(N.F.nodes[1])
              }
            }
            lp.N1 <- lp.N1.not1
          }
          #both sex-specific recruitment likelihoods can change because ER.M and ER.F depend on total N
          lp.N.recruit.M <- model$calculate(N.recruit.M.nodes[1:(first.det-1)])
          lp.N.recruit.F <- model$calculate(N.recruit.F.nodes[1:(first.det-1)])
          #only observation likelihoods before first detection can change
          lp.y <- model$calculate(y.nodes[y.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #after removing this individual, all unchanged multinomial terms cancel
          cohort.prop <- g + model$sex[i]*n.primary
          lp.prior <- log(entry.counts.minus[cohort.prop]+1)
          lp.start[g] <- lp.N1 + lp.N.recruit.M + lp.N.recruit.F + lp.y + lp.surv + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        
        if(model$z.start[i]!=z.start.prop){ #if proposal is same as current, no need to replace anything
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
          #now repeat for sex
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
            model$N.recruit.M <<- N.recruit.M.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- N.recruit.M.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- N.recruit.M.curr[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
            model$N.recruit.F <<- N.recruit.F.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- N.recruit.F.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- N.recruit.F.curr[z.start.prop-1] + 1
            }
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          }
          model$calculate(ER.M.nodes[1:(first.det-1)])
          model$calculate(ER.F.nodes[1:(first.det-1)])
          #update these logProbs
          model$calculate(y.nodes[y.idx])
          model$calculate(z.nodes[i])
          if(model$sex[i]==0){
            model$calculate(N.M.nodes[1])
          }else{
            model$calculate(N.F.nodes[1])
          }
          model$calculate(N.recruit.M.nodes[1:(first.det-1)])
          model$calculate(N.recruit.F.nodes[1:(first.det-1)])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          if(model$sex[i]==0){ #male
            mvSaved["N.M",1] <<- model[["N.M"]]
            mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
            mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
          }else{ #female
            mvSaved["N.F",1] <<- model[["N.F"]]
            mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
            mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
          }
          mvSaved["ER.M",1] <<- model[["ER.M"]]
          mvSaved["ER.F",1] <<- model[["ER.F"]]
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          cohort.prop <- z.start.prop + model$sex[i]*n.primary
          entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
          entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          if(model$sex[i]==0){ #male
            model[["N.M"]] <<- mvSaved["N.M",1]
            model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
            model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
          }else{ #female
            model[["N.F"]] <<- mvSaved["N.F",1]
            model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
            model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
          }
          model[["ER.M"]] <<- mvSaved["ER.M",1]
          model[["ER.F"]] <<- mvSaved["ER.F",1]
          #set these logProbs back
          model$calculate(y.nodes[y.idx])
          if(model$sex[i]==0){
            model$calculate(N.M.nodes[1])
          }else{
            model$calculate(N.F.nodes[1])
          }
          model$calculate(N.recruit.M.nodes[1:(first.det-1)])
          model$calculate(N.recruit.F.nodes[1:(first.det-1)])
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
        if(model$sex[i]==0){
          N.M.curr <- model$N.M
        }else{
          N.F.curr <- model$N.F
        }
        dets <- which(model$y[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.primary)
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
        #only y nodes after last detection can change across z.stop candidates
        y.idx <- i.idx[(last.det+1):n.primary]
        for(g in (last.det):n.primary){ #can't die on or before primary occasion of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          #Dont need to update N.survive until we select a state--does not change logProb
          model$N <<- N.curr - z.curr + z.prop
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
          }
          #Updating total N changes both ER.M and ER.F, but only intervals after last detection can change
          if(last.det < n.primary-1){
            model$calculate(ER.M.nodes[(last.det+1):(n.primary-1)])
            model$calculate(ER.F.nodes[(last.det+1):(n.primary-1)])
          }
          #N.M[1] and N.F[1] cannot change in a z.stop update, so their likelihoods cancel
          if(last.det < n.primary-1){
            lp.N.recruit.M <- model$calculate(N.recruit.M.nodes[(last.det+1):(n.primary-1)])
            lp.N.recruit.F <- model$calculate(N.recruit.F.nodes[(last.det+1):(n.primary-1)])
          }else{
            lp.N.recruit.M <- 0
            lp.N.recruit.F <- 0
          }
          lp.y <- model$calculate(y.nodes[y.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N.recruit.M + lp.N.recruit.F + lp.y + lp.surv
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
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F
          }
          if(last.det < n.primary-1){
            model$calculate(ER.M.nodes[(last.det+1):(n.primary-1)])
            model$calculate(ER.F.nodes[(last.det+1):(n.primary-1)])
          }
          #update these logProbs
          if(last.det < n.primary-1){
            model$calculate(N.recruit.M.nodes[(last.det+1):(n.primary-1)])
            model$calculate(N.recruit.F.nodes[(last.det+1):(n.primary-1)])
          }
          model$calculate(y.nodes[y.idx])
          model$calculate(z.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          if(model$sex[i]==0){ #male
            mvSaved["N.M",1] <<- model[["N.M"]]
            mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
          }else{ #female
            mvSaved["N.F",1] <<- model[["N.F"]]
            mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
          }
          mvSaved["ER.M",1] <<- model[["ER.M"]]
          mvSaved["ER.F",1] <<- model[["ER.F"]]
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          if(model$sex[i]==0){ #male
            model[["N.M"]] <<- mvSaved["N.M",1]
            model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
          }else{ #female
            model[["N.F"]] <<- mvSaved["N.F",1]
            model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
          }
          model[["ER.M"]] <<- mvSaved["ER.M",1]
          model[["ER.F"]] <<- mvSaved["ER.F",1]
          #set these logProbs back
          if(last.det < n.primary-1){
            model$calculate(N.recruit.M.nodes[(last.det+1):(n.primary-1)])
            model$calculate(N.recruit.F.nodes[(last.det+1):(n.primary-1)])
          }
          model$calculate(y.nodes[y.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    
    # 2) undetected guy update. Only if in the superpopulation. 
    # Metropolis-Hastings, simulate new sex + z vector from priors
    # entry counts current after z.start update
    for(i in 1:M){
      if(z.obs[i]==0&model$z.super[i]==1){
        sex.curr <- model$sex[i] #store this for use below
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        cohort.curr <- model$z.start[i] + model$sex[i]*n.primary
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y nodes
        i.idx2 <- seq(i,M*(n.primary-1),M) #used to reference correct phi nodes
        #get forwards recruitment probabilities
        #paste male and female recruit probs
        recruit.probs.for <- c(model$lambda.y1.M,model$ER.M,model$lambda.y1.F,model$ER.F)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #track proposal probs
        #survival proposal probabilities cancel exactly with the survival likelihood because
        #the survival history is proposed from the same sex-specific survival model used in the target
        log.prop.for <- log.prop.back <- 0
        
        #simulate recruitment, update proposed sex
        cohort.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.primary)
        if(cohort.prop<=n.primary){ #simulated male
          z.start.prop <- cohort.prop
          model$sex[i] <<- 0
        }else{ #simulated female
          z.start.prop <- cohort.prop - n.primary
          model$sex[i] <<- 1
        }
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[cohort.prop])
        
        #update phi because sex can change
        model$calculate(phi.nodes[i.idx2])
        #simulate survival with proposed-sex phi; once dead, remaining z's stay 0
        z.stop.prop <- z.start.prop
        if(z.start.prop < n.primary){
          for(g in (z.start.prop+1):n.primary){
            if(z.prop[g-1]==1){
              z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
              if(z.prop[g]==1){
                z.stop.prop <- g
              }
            }
          }
        }
        
        #cohort encodes both sex and entry occasion, so this catches sex changes as well as z changes
        if(cohort.prop!=cohort.curr|z.stop.prop!=z.stop.curr){
          #get initial stored logProbs before recalculating any sex-dependent likelihood nodes
          lp.initial.entry.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.entry.M <- lp.initial.entry.M + model$getLogProb(N.recruit.M.nodes)
          lp.initial.entry.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.entry.F <- lp.initial.entry.F + model$getLogProb(N.recruit.F.nodes)
          # lp.initial.y <- model$getLogProb(y.nodes[i.idx])
          #restrict y nodes to those where z and/or sex changed
          if(model$sex[i]!=sex.curr){
            y.idx.changed <- i.idx
          }else{
            z.changed <- which(z.prop!=z.curr)
            y.idx.changed <- i.idx[z.changed]
          }
          lp.initial.y <- model$getLogProb(y.nodes[y.idx.changed])
          #lp.initial.surv <- model$getLogProb(z.nodes[i]) #cancels exactly with backwards survival proposal probability
          #log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1))) #full multinomial coefficient replaced by exact ratio below
          
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
          #repeat for sex
          if(sex.curr==0&model$sex[i]==0){ #male to male
            model$N.M <<- model$N.M - z.curr + z.prop
            if(z.start.curr > 1){
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M
          }else if(sex.curr==1&model$sex[i]==1){ #female to female
            model$N.F <<- model$N.F - z.curr + z.prop
            if(z.start.curr > 1){
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F
          }else if(sex.curr==0&model$sex[i]==1){ #male to female
            model$N.M <<- model$N.M - z.curr
            model$N.F <<- model$N.F + z.prop
            if(z.start.curr > 1){
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F
          }else if(sex.curr==1&model$sex[i]==0){ #female to male
            model$N.F <<- model$N.F - z.curr
            model$N.M <<- model$N.M + z.prop
            if(z.start.curr > 1){
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F
          }
          model$calculate(ER.M.nodes) #update ER when N updated
          model$calculate(ER.F.nodes) #update ER when N updated
          #get proposed logProbs
          lp.proposed.entry.M <- model$calculate(N.M.nodes[1])
          lp.proposed.entry.M <- lp.proposed.entry.M + model$calculate(N.recruit.M.nodes)
          lp.proposed.entry.F <- model$calculate(N.F.nodes[1])
          lp.proposed.entry.F <- lp.proposed.entry.F + model$calculate(N.recruit.F.nodes)
          # lp.proposed.y <- model$calculate(y.nodes[i.idx])
          lp.proposed.y <- model$calculate(y.nodes[y.idx.changed])
          #lp.proposed.surv <- model$calculate(z.nodes[i]) #cancels exactly with forwards survival proposal probability
          
          #update entry counts and use the exact local multinomial coefficient ratio
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
          entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1
          if(cohort.prop!=cohort.curr){
            log.prior.ratio <- log(entry.counts.curr[cohort.prop]+1)-log(entry.counts.curr[cohort.curr])
          }else{
            log.prior.ratio <- 0
          }
          
          #get backwards proposal probability for the original sex/entry cohort
          recruit.probs.back <- c(model$lambda.y1.M,model$ER.M,model$lambda.y1.F,model$ER.F)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log.prop.back + log(recruit.probs.back[cohort.curr])
          #survival proposal probabilities are not calculated because, using the appropriate
          #current/proposed sex-specific phi in each direction, they cancel the survival likelihood ratio exactly
          
          lp.initial.total <- lp.initial.entry.M + lp.initial.entry.F + lp.initial.y
          lp.proposed.total <- lp.proposed.entry.M + lp.proposed.entry.F + lp.proposed.y
          
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.prior.ratio + log.prop.back) -
            (lp.initial.total + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #update survival logProb once for accepted history; it was not needed when evaluating MH ratio
            model$calculate(z.nodes[i])
            mvSaved["z.start",1][i] <<- model[["z.start"]][i]
            mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
            mvSaved["z",1][i,] <<- model[["z"]][i,]
            mvSaved["sex",1][i] <<- model[["sex"]][i]
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.M",1] <<- model[["N.M"]]
            mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
            mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
            mvSaved["N.F",1] <<- model[["N.F"]]
            mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
            mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
            for(g in 1:(n.primary-1)){
              mvSaved["phi",1][i,g] <<- model[["phi"]][i,g]
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][i] <<- mvSaved["z.start",1][i]
            model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
            model[["z"]][i,] <<- mvSaved["z",1][i,]
            model[["sex"]][i] <<- mvSaved["sex",1][i]
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.M"]] <<- mvSaved["N.M",1]
            model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
            model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
            model[["N.F"]] <<- mvSaved["N.F",1]
            model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
            model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
            #set these logProbs back
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
            model$calculate(phi.nodes[i.idx2])
            #model$calculate(z.nodes[i]) #not needed because survival logProb was never recalculated for the proposal
            model$calculate(y.nodes[y.idx.changed])
          }
        }
      }
    }
    
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
      #propose to add/subtract 1; direction probability is symmetric and cancels
      updown <- rbinom(1,1,0.5)
      if(updown==0){#subtract
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          sex.curr <- model$sex[pick]
          cohort.curr <- z.start.curr + sex.curr*n.primary
          pick.idx <- seq(pick,M*n.primary,M) #used to reference correct y nodes
          
          #p select on undetected guy
          log.p.select.for <- log(1/non.init)
          
          #get initial logprobs
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
          #lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
          #only currently alive primary occasions can change when this individual is removed
          lp.initial.y <- 0
          for(g in 1:n.primary){
            if(z.curr[g]==1){
              lp.initial.y <- lp.initial.y+model$getLogProb(y.nodes[pick.idx[g]])
            }
          }
          #lp.initial.surv <- model$getLogProb(z.nodes[pick]) #cancels exactly with reverse survival proposal probability
          
          # propose new N.super/z.super
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
          #repeat for sex
          if(sex.curr==0){
            model$N.M <<- model$N.M - z.curr
            if(z.start.curr > 1){
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            model$N.survive.M <<- model$N.M[2:n.primary]-model$N.recruit.M
          }else{
            model$N.F <<- model$N.F - z.curr
            if(z.start.curr > 1){
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            model$N.survive.F <<- model$N.F[2:n.primary]-model$N.recruit.F
          }
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          
          #Reverse proposal probability: only sex/entry-cohort probability is needed.
          #The reverse survival proposal is exactly the current survival likelihood and cancels.
          recruit.probs.back <- c(model$lambda.y1.M, model$ER.M, model$lambda.y1.F, model$ER.F)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[cohort.curr])
          #if(z.start.curr < n.primary){
          #  for(g in (z.start.curr+1):n.primary){
          #    log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
          #  }
          #}
          
          #get proposed logprobs for N and y
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          #lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0
          lp.proposed.y <- 0 #all focal observation likelihood terms are zero when z.super=0
          #lp.proposed.surv <- model$calculate(z.nodes[pick]) #cancels exactly with reverse survival proposal probability
          
          lp.initial.total <- lp.initial.N.M + lp.initial.N.F + lp.initial.y + lp.initial.N.recruit.M + lp.initial.N.recruit.F
          lp.proposed.total <- lp.proposed.N.M + lp.proposed.N.F + lp.proposed.y + lp.proposed.N.recruit.M + lp.proposed.N.recruit.F
          
          #local multinomial coefficient ratio for moving current sex/entry cohort to z.super=0
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
          entry.counts.prop[2*n.primary + 1] <- entry.counts.prop[2*n.primary + 1] + 1
          log.z.prior.ratio <- log(entry.counts.curr[2*n.primary+1]+1)-log(entry.counts.curr[cohort.curr])
          
          #p select this individual in the reverse add move
          noff.back <- noff.curr+1
          log.p.select.back <- log(1/noff.back)
          log.prop.for <- 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) - (lp.initial.total + log.p.select.for + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #survival logProb was omitted from the MH calculation because it cancels with the proposal;
            #calculate it once now to synchronize the accepted model state
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
            if(sex.curr==0){
              mvSaved["N.M",1] <<- model[["N.M"]]
              mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
              mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
            }else{
              mvSaved["N.F",1] <<- model[["N.F"]]
              mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
              mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
            }
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
            entry.counts.curr <- entry.counts.prop
            #move accepted individual from on list to off list using swap-delete
            z.on[pick.pos] <- z.on[non.curr]
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
            if(sex.curr==0){
              model[["N.M"]] <<- mvSaved["N.M",1]
              model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
              model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
            }else{
              model[["N.F"]] <<- mvSaved["N.F",1]
              model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
              model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
            }
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            #set these logProbs back
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
            #model$calculate(y.nodes[pick.idx])
            #not needed: focal y logProbs were not recalculated for the rejected removal proposal
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
            #model$calculate(z.nodes[pick]) #not needed because survival logProb was never recalculated for the proposal
          }
        }
      }else{#add
        noff.init <- noff.curr
        if(noff.init>0){
          pick.pos <- rcat(1,rep(1/noff.init,noff.init))
          pick <- z.off[pick.pos]
          pick.idx <- seq(pick,M*n.primary,M)
          pick.idx2 <- seq(pick,M*(n.primary-1),M)
          
          #p select off undetected guy
          log.p.select.for <- log(1/noff.init)
          
          #get initial logProbs
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
          #lp.initial.y <- model$getLogProb(y.nodes[pick.idx]) #will always be 0
          lp.initial.y <- 0 #all focal observation likelihood terms are zero when z.super=0
          #lp.initial.surv <- model$getLogProb(z.nodes[pick]) #cancels exactly with forward survival proposal probability
          
          #propose sex and entry cohort
          recruit.probs.for <- c(model$lambda.y1.M, model$ER.M, model$lambda.y1.F, model$ER.F)
          recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
          cohort.prop <- rcat(1,recruit.probs.for)
          log.prop.for <- log(recruit.probs.for[cohort.prop])
          if(cohort.prop <= n.primary){
            sex.prop <- 0
            z.start.prop <- cohort.prop
          }else{
            sex.prop <- 1
            z.start.prop <- cohort.prop - n.primary
          }
          
          #set proposed sex first, then update phi
          model$sex[pick] <<- sex.prop
          model$calculate(phi.nodes[pick.idx2])
          
          #simulate proposed z vector from the proposed-sex survival model
          #survival proposal probabilities are not calculated because they cancel the survival likelihood exactly
          z.prop <- rep(0,n.primary)
          z.prop[z.start.prop] <- 1
          z.stop.prop <- z.start.prop
          if(z.start.prop < n.primary){
            for(g in (z.start.prop+1):n.primary){
              if(z.prop[g-1]==1){
                z.prop[g] <- rbinom(1,1,model$phi[pick,g-1]*z.prop[g-1])
                if(z.prop[g]==1){
                  z.stop.prop <- g
                }
              }
            }
          }
          
          #store in model
          model$z.super[pick] <<- 1
          model$N.super <<- model$N.super + 1
          model$z.start[pick] <<- z.start.prop
          model$z.stop[pick] <<- z.stop.prop
          model$z[pick,] <<- z.prop
          
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N + model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          #repeat for sex
          if(sex.prop==0){
            model$N.M <<- model$N.M + model$z[pick,]
            if(z.start.prop > 1){
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.primary] - model$N.recruit.M
          }else{
            model$N.F <<- model$N.F + model$z[pick,]
            if(z.start.prop > 1){
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            model$N.survive.F <<- model$N.F[2:n.primary] - model$N.recruit.F
          }
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #get proposed logprobs for N and y
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          #lp.proposed.y <- model$calculate(y.nodes[pick.idx])
          lp.proposed.y <- 0
          for(g in 1:n.primary){
            if(model$z[pick,g]==1){
              lp.proposed.y <- lp.proposed.y+model$calculate(y.nodes[pick.idx[g]])
            }
          }
          #lp.proposed.surv <- model$calculate(z.nodes[pick]) #cancels exactly with forward survival proposal probability
          
          lp.initial.total <- lp.initial.N.M + lp.initial.N.F + lp.initial.y + lp.initial.N.recruit.M + lp.initial.N.recruit.F
          lp.proposed.total <- lp.proposed.N.M + lp.proposed.N.F + lp.proposed.y + lp.proposed.N.recruit.M + lp.proposed.N.recruit.F
          
          #local multinomial coefficient ratio for moving z.super=0 to proposed sex/entry cohort
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[2*n.primary + 1] <- entry.counts.prop[2*n.primary + 1] - 1
          entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1
          log.z.prior.ratio <- log(entry.counts.curr[cohort.prop]+1)-log(entry.counts.curr[2*n.primary+1])
          
          #p select this individual in the reverse subtract move
          non.back <- non.curr+1
          log.p.select.back <- log(1/non.back)
          log.prop.back <- 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) - (lp.initial.total + log.p.select.for + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #survival logProb was omitted from the MH calculation because it cancels with the proposal;
            #calculate it once now to synchronize the accepted model state
            model$calculate(z.nodes[pick])
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["sex",1][pick] <<- model[["sex"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            if(sex.prop==0){
              mvSaved["N.M",1] <<- model[["N.M"]]
              mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
              mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
            }else{
              mvSaved["N.F",1] <<- model[["N.F"]]
              mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
              mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
            }
            for(g in 1:(n.primary-1)){
              mvSaved["phi",1][pick,g] <<- model[["phi"]][pick,g]
            }
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
            entry.counts.curr <- entry.counts.prop
            #move accepted individual from off list to on list using swap-delete
            z.off[pick.pos] <- z.off[noff.curr]
            noff.curr <- noff.curr-1
            non.curr <- non.curr+1
            z.on[non.curr] <- pick
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["sex"]][pick] <<- mvSaved["sex",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            if(sex.prop==0){
              model[["N.M"]] <<- mvSaved["N.M",1]
              model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
              model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
            }else{
              model[["N.F"]] <<- mvSaved["N.F",1]
              model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
              model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
            }
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            #restore phi for original off-state sex before recalculating sex-dependent y
            model$calculate(phi.nodes[pick.idx2])
            #set these logProbs back
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
            #model$calculate(y.nodes[pick.idx])
            #only the proposed alive primary occasions had their focal y logProbs changed
            for(g in z.start.prop:z.stop.prop){
              model$calculate(y.nodes[pick.idx[g]])
            }
            #model$calculate(z.nodes[pick]) #not needed because survival logProb was never recalculated for the proposal
          }
        }
      }
    }
    #4) Finally, detected guy unobserved sex update
    for(i in 1:length(sex.up)){
      if(z.obs[sex.up[i]]==1){ #only do detected guys here
        i.idx <- seq(sex.up[i],M*n.primary,M) #used to reference correct y nodes
        i.idx2 <- seq(sex.up[i],M*(n.primary-1),M) #used to reference correct phi nodes
        z.start.curr <- model$z.start[sex.up[i]]
        z.stop.curr <- model$z.stop[sex.up[i]]
        #only observation likelihoods while the individual is alive can change with sex
        y.idx <- i.idx[z.start.curr:z.stop.curr]
        cohort.curr <- z.start.curr + model$sex[sex.up[i]]*n.primary
        
        #get initial logProbs
        #Only N.M[1]/N.F[1] change when the individual entered in primary occasion 1.
        #Otherwise, only the male/female recruitment likelihood at the entry occasion changes.
        if(z.start.curr==1){
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- 0
          lp.initial.N.recruit.F <- 0
        }else{
          lp.initial.N.M <- 0
          lp.initial.N.F <- 0
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes[z.start.curr-1])
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes[z.start.curr-1])
        }
        lp.initial.z <- model$getLogProb(z.nodes[sex.up[i]])
        lp.initial.y <- model$getLogProb(y.nodes[y.idx])
        #full multinomial coefficient calculation replaced by exact local ratio below
        #lp.initial.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))
        
        #update N variables
        if(model$sex[sex.up[i]]==0){ #initial male
          #move this guy from N.M to N.F
          model$N.M <<- model$N.M - model$z[sex.up[i],]
          model$N.F <<- model$N.F + model$z[sex.up[i],]
          #move male recruit to female recruit
          if(z.start.curr>1){ #otherwise, this is primary occasion 1, nothing to change
            model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] + 1
          }
          #update male, female survivors
          model$N.survive.M <<- model$N.M[2:n.primary] - model$N.recruit.M
          model$N.survive.F <<- model$N.F[2:n.primary] - model$N.recruit.F
        }else{ #initial female
          #move this guy from N.F to N.M
          model$N.F <<- model$N.F - model$z[sex.up[i],]
          model$N.M <<- model$N.M + model$z[sex.up[i],]
          #move female recruit to male recruit
          if(z.start.curr>1){ #otherwise, this is primary occasion 1, nothing to change
            model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] + 1
          }
          #update male, female survivors
          model$N.survive.F <<- model$N.F[2:n.primary] - model$N.recruit.F
          model$N.survive.M <<- model$N.M[2:n.primary] - model$N.recruit.M
        }
        #update sex
        model$sex[sex.up[i]] <<- 1 - model$sex[sex.up[i]]
        #update phi nodes when sex changes
        model$calculate(phi.nodes[i.idx2])
        #update prior
        cohort.prop <- z.start.curr + model$sex[sex.up[i]]*n.primary
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
        entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1
        
        #get proposed logProbs
        if(z.start.curr==1){
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- 0
          lp.proposed.N.recruit.F <- 0
        }else{
          lp.proposed.N.M <- 0
          lp.proposed.N.F <- 0
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes[z.start.curr-1])
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes[z.start.curr-1])
        }
        lp.proposed.z <- model$calculate(z.nodes[sex.up[i]])
        #lp.proposed.y <- model$calculate(y.nodes[i.idx])
        lp.proposed.y <- model$calculate(y.nodes[y.idx])
        #full multinomial coefficient calculations are unnecessary. Moving one individual from
        #cohort a to cohort b gives log prior ratio log(n.b+1)-log(n.a).
        #lp.proposed.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
        log.prior.ratio <- log(entry.counts.curr[cohort.prop]+1)-log(entry.counts.curr[cohort.curr])
        
        lp.initial.total <- lp.initial.z + lp.initial.y + lp.initial.N.recruit.M +
          lp.initial.N.recruit.F + lp.initial.N.M + lp.initial.N.F
        lp.proposed.total <- lp.proposed.z + lp.proposed.y + lp.proposed.N.recruit.M +
          lp.proposed.N.recruit.F + lp.proposed.N.M + lp.proposed.N.F
        
        #MH step
        log_MH_ratio <- lp.proposed.total - lp.initial.total + log.prior.ratio
        accept <- decide(log_MH_ratio)
        if(accept){
          mvSaved["N.M",1] <<- model[["N.M"]]
          mvSaved["N.F",1] <<- model[["N.F"]]
          mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
          mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
          mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
          mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
          mvSaved["sex",1][sex.up[i]] <<- model[["sex"]][sex.up[i]]
          for(g in 1:(n.primary-1)){
            mvSaved["phi",1][sex.up[i],g] <<- model[["phi"]][sex.up[i],g]
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["N.M"]] <<- mvSaved["N.M",1]
          model[["N.F"]] <<- mvSaved["N.F",1]
          model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
          model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
          model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
          model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
          model[["sex"]][sex.up[i]] <<- mvSaved["sex",1][sex.up[i]]
          for(g in 1:(n.primary-1)){
            model[["phi"]][sex.up[i],g] <<- mvSaved["phi",1][sex.up[i],g]
          }
          #only restore logProbs that were recalculated for the proposal
          model$calculate(y.nodes[y.idx])
          model$calculate(phi.nodes[i.idx2])
          model$calculate(z.nodes[sex.up[i]])
          if(z.start.curr==1){
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
          }else{
            model$calculate(N.recruit.M.nodes[z.start.curr-1])
            model$calculate(N.recruit.F.nodes[z.start.curr-1])
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
    n.recruit <- length(model$N.recruit.M)
    sex.ind <- 1
    g <- 1
    if(target=="lambda.y1.M"){
      sampler.type <- 1
    }else if(target=="lambda.y1.F"){
      sampler.type <- 2
    }else if(grepl("^gamma\\.sex\\[[12]\\]$",target)){
      sampler.type <- 3
      sex.ind <- as.integer(gsub("[^0-9]","",target))
    }else{
      sampler.type <- 4
      ind <- as.integer(strsplit(gsub("gamma\\.sex\\[|\\]","",target),",")[[1]])
      sex.ind <- ind[1]
      g <- ind[2]
    }
  },
  run = function(){
    if(sampler.type==1){
      count <- model$N.M[1]
      rate <- 1
    }else if(sampler.type==2){
      count <- model$N.F[1]
      rate <- 1
    }else if(sampler.type==3){
      count <- 0
      rate <- 0
      for(j in 1:n.recruit){
        if(sex.ind==1){
          count <- count+model$N.recruit.M[j]
        }else{
          count <- count+model$N.recruit.F[j]
        }
        rate <- rate+model$N[j]*tau[j]
      }
    }else{
      if(sex.ind==1){
        count <- model$N.recruit.M[g]
      }else{
        count <- model$N.recruit.F[g]
      }
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
    sex.specific <- control$sex.specific
    occasion.specific <- control$occasion.specific
    shape1 <- model$getParam(target,"shape1")
    shape2 <- model$getParam(target,"shape2")
    idx <- as.integer(strsplit(gsub("[^0-9,]","",target),",")[[1]])
    s <- idx[1]
    g <- idx[2]
  },
  run = function(){
    success <- 0
    failure <- 0
    for(j in 1:n.primary){
      if(!occasion.specific|j==g){
        for(i in 1:M){
          if(model$z.super[i]==1 & model$z[i,j]==1){
            if(!sex.specific | model$sex[i]+1==s){
              success <- success+model$y[i,j]
              failure <- failure+K[j]-model$y[i,j]
            }
          }
        }
      }
    }
    model[[target]] <<- rbeta(1,shape1+success,shape2+failure)
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)