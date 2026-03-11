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
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0)) {
    returnType(double(1))
    n.year <- length(phi)
    return(rep(0,n.year))
  }
)

#custom observation model distribution that makes custom updates easier
dbinomial2 <- nimbleFunction(
  run = function(x = double(0), p = double(0), K = double(0), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this year
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
    n.year <- control$n.year
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
    entry.counts.curr <- rep(0, 2*n.year+1)
    for(g in 1:n.year){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1 & model$sex==0)
      entry.counts.curr[g + n.year] <- sum(model$z.start==g & model$z.super==1 & model$sex==1)
    }
    entry.counts.curr[2*n.year + 1] <- sum(model$z.super==0)
    
    # 1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,1]==0){ #for detected guys, skip if observed 1st year
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
        lp.start <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y nodes
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

          #now repeat for sex
          if(model$sex[i]==0){ #male
            #1) Update N
            model$N.M <<- N.M.curr - z.curr + z.prop
            #2) Update N.recruit
            model$N.recruit.M <<- N.recruit.M.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- N.recruit.M.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- N.recruit.M.curr[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{ #female
            #1) Update N
            model$N.F <<- N.F.curr - z.curr + z.prop
            #2) Update N.recruit
            model$N.recruit.F <<- N.recruit.F.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- N.recruit.F.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- N.recruit.F.curr[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          #must account for sex-specificity, Updating N changes both ER.M and ER.F
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #get these logProbs
          lp.N1.M <- model$calculate(N.M.nodes[1])
          lp.N1.F <- model$calculate(N.F.nodes[1])
          lp.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          lp.y <- model$calculate(y.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          if(model$sex[i]==0){
            entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
            entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          }else{
            entry.counts.prop[z.start.curr + n.year] <- entry.counts.prop[z.start.curr + n.year] - 1
            entry.counts.prop[z.start.prop + n.year] <- entry.counts.prop[z.start.prop + n.year] + 1
          }
          lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          lp.start[g] <- lp.N1.M + lp.N1.F +
            lp.N.recruit.M + lp.N.recruit.F +
            lp.y + lp.surv + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)

        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original

        if(model$z.start[i]!=z.start.prop){ #if proposal is same as current, no need to replace anything
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
          #now repeat for sex
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
            model$N.recruit.M <<- N.recruit.M.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- N.recruit.M.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- N.recruit.M.curr[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
            model$N.recruit.F <<- N.recruit.F.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- N.recruit.F.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- N.recruit.F.curr[z.start.prop-1] + 1
            }
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #update these logProbs
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
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
          if(model$sex[i]==0){
            entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
            entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          }else{
            entry.counts.prop[z.start.curr + n.year] <- entry.counts.prop[z.start.curr + n.year] - 1
            entry.counts.prop[z.start.prop + n.year] <- entry.counts.prop[z.start.prop + n.year] + 1
          }
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
          model$calculate(y.nodes[i.idx])
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(z.nodes[i])
        }
      }
    }

    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,n.year]==0){ #for detected guys, skip if observed in final year
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
        lp.stop <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y nodes
        for(g in (last.det):n.year){ #can't die on or before year of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.year)
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
          model$calculate(ER.M.nodes) #update ER when N.M updated
          model$calculate(ER.F.nodes) #update ER when N.F updated
          lp.N1.M <- model$calculate(N.M.nodes[1])
          lp.N1.F <- model$calculate(N.F.nodes[1])
          lp.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          lp.y <- model$calculate(y.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N1.M + lp.N1.F +
            lp.N.recruit.M + lp.N.recruit.F +
            lp.y + lp.surv
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
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F
          }
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #update these logProbs
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(y.nodes[i.idx])
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
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #set these logProbs back
          model$calculate(y.nodes[i.idx])
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
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
        cohort.curr <- model$z.start[i] + model$sex[i]*n.year
        i.idx <- seq(i,M*n.year,M) #used to reference correct y nodes
        i.idx2 <- seq(i,M*(n.year-1),M) #used to reference correct phi nodes
        #get forwards recruitment probabilities
        #paste male and female recruit probs
        recruit.probs.for <- c(model$lambda.y1.M,model$ER.M,model$lambda.y1.F,model$ER.F)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #get initial logProbs
        lp.initial.entry.M <- model$getLogProb(N.M.nodes[1])
        lp.initial.entry.M <- lp.initial.entry.M + model$getLogProb(N.recruit.M.nodes)
        lp.initial.entry.F <- model$getLogProb(N.F.nodes[1])
        lp.initial.entry.F <- lp.initial.entry.F + model$getLogProb(N.recruit.F.nodes)
        lp.initial.y <- model$getLogProb(y.nodes[i.idx])
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))

        #track proposal probs - survival is symmetric, but not recruitment and detection
        log.prop.for <- log.prop.back <- 0
        phi.curr <- model$phi[i,] #store this for backwards proposal probs

        #simulate recruitment, update z.start and sex
        cohort.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.year)
        if(cohort.prop<=n.year){ #simulated male
          z.start.prop <- cohort.prop
          model$sex[i] <<- 0
        }else{ #simulated female
          z.start.prop <- cohort.prop - n.year
          model$sex[i] <<- 1
        }
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[cohort.prop])

        #update phi bc sex can change
        model$calculate(phi.nodes[i.idx2])
        #simulate survival with updated phi
        if(z.start.prop < n.year){ #if you don't recruit in final year
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

        #update N, N.recruit, N.survive
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
        #repeat for sex
        if(sex.curr==0&model$sex[i]==0){ #male to male
          model$N.M <<- model$N.M - z.curr + z.prop
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
          }
          model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
        }else if(sex.curr==1&model$sex[i]==1){ #female to female
          model$N.F <<- model$N.F - z.curr + z.prop
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
          }
          model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
        }else if(sex.curr==0&model$sex[i]==1){ #male to female
          #subtract current z from males, add new z to females
          model$N.M <<- model$N.M - z.curr
          model$N.F <<- model$N.F + z.prop
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
          }
          model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g

        }else if(sex.curr==1&model$sex[i]==0){ #female to male
          #subtract current z from females, add new z to males
          model$N.F <<- model$N.F - z.curr
          model$N.M <<- model$N.M + z.prop
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
          }
          model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
        }
        model$calculate(ER.M.nodes) #update ER when N updated
        model$calculate(ER.F.nodes) #update ER when N updated
        #get proposed logProbs
        lp.proposed.entry.M <- model$calculate(N.M.nodes[1])
        lp.proposed.entry.M <- lp.proposed.entry.M + model$calculate(N.recruit.M.nodes)
        lp.proposed.entry.F <- model$calculate(N.F.nodes[1])
        lp.proposed.entry.F <- lp.proposed.entry.F + model$calculate(N.recruit.F.nodes)
        lp.proposed.y <- model$calculate(y.nodes[i.idx])
        lp.proposed.surv <- model$calculate(z.nodes[i])

        # Full multinomial coefficient prior for proposed configuration
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
        entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1
        log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))

        #get backwards proposal probs
        recruit.probs.back <- c(model$lambda.y1.M,model$ER.M,model$lambda.y1.F,model$ER.F)
        recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)

        log.prop.back <- log.prop.back + log(recruit.probs.back[cohort.curr])
        if(z.start.curr < n.year){#if you don't recruit in final year
          for(g in (z.start.curr+1):n.year){
            #use original phi.curr stored above
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,phi.curr[g-1]*z.curr[g-1],log=TRUE)
          }
        }

        #Add up likelihoods and prop probs
        lp.initial.total <- lp.initial.entry.M + lp.initial.entry.F + lp.initial.y +
          lp.initial.surv + log.prior.curr
        lp.proposed.total <- lp.proposed.entry.M + lp.proposed.entry.F + lp.proposed.y +
          lp.proposed.surv + log.prior.prop

        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        accept <- decide(log_MH_ratio)
        if(accept){
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
          model$calculate(z.nodes[i])
          model$calculate(y.nodes[i.idx])
        }
      }
    }
    #3) update z.super: Metropolis-Hastings
    # entry counts current coming out of undetected ind update
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z.super==1)
        non.init <- length(z.on)
        pick <- rcat(1,rep(1/non.init,non.init)) #select one of these individuals
        pick <- z.on[pick]
        if(z.obs[pick]==1){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          sex.curr <- model$sex[pick]
          if(sex.curr==0){
            cohort.curr <- z.start.curr
          }else{
            cohort.curr <- z.start.curr + n.year
          }

          #p select off guy
          log.p.select.for <- log(1/non.init)
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx <- seq(pick,M*n.year,M) #used to reference correct y nodes

          #get initial logprobs
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
          lp.initial.surv <- model$getLogProb(z.nodes[pick])

          # propose new N.super/z.super
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
          #repeat for sex
          if(sex.curr==0){
            model$N.M <<- model$N.M - z.curr
            if(z.start.curr > 1){ #if wasn't in pop in year 1
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{
            model$N.F <<- model$N.F - z.curr
            if(z.start.curr > 1){ #if wasn't in pop in year 1
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          model$calculate(ER.M.nodes) #update ER when N updated
          model$calculate(ER.F.nodes) #update ER when N updated

          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1.M, model$ER.M,
                                  model$lambda.y1.F, model$ER.F)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[cohort.curr])
          if(z.start.curr < n.year){
            for(g in (z.start.curr + 1):n.year){
              log.prop.back <- log.prop.back +
                dbinom(z.curr[g], 1, model$phi[pick, g-1] * z.curr[g-1], log = TRUE)
            }
          }

          #get proposed logprobs for N and y
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0
          lp.proposed.surv <- model$calculate(z.nodes[pick]) #will always be 0

          lp.initial.total <- lp.initial.N.M + lp.initial.N.F + lp.initial.y +
            lp.initial.N.recruit.M + lp.initial.N.recruit.F + lp.initial.surv
          lp.proposed.total <- lp.proposed.N.M + lp.proposed.N.F + lp.proposed.y +
            lp.proposed.N.recruit.M + lp.proposed.N.recruit.F + lp.proposed.surv

          #backwards prior and select probs
          #move from current z.super=1 cohort to z.super=0 cell
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
          entry.counts.prop[2*n.year + 1] <- entry.counts.prop[2*n.year + 1] + 1
          
          #p select on guy
          noff.back <- sum(model$z.super == 0)
          log.p.select.back <- log(1/noff.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.for <- 0

          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)
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
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
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
            model$calculate(y.nodes[pick.idx])
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
            model$calculate(z.nodes[pick])
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          z.off <- which(model$z.super==0)
          noff.init <- length(z.off)
          pick <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx <- seq(pick,M*n.year,M)
          pick.idx2 <- seq(pick,M*(n.year-1),M)

          #p select off guy
          log.p.select.for <- log(1/noff.init)

          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))

          #get initial logProbs
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx]) #will always be 0
          lp.initial.surv <- model$getLogProb(z.nodes[pick]) #will always be 0

          #propose sex and entry cohort
          recruit.probs.for <- c(model$lambda.y1.M, model$ER.M,
                                 model$lambda.y1.F, model$ER.F)
          recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
          cohort.prop <- rcat(1, recruit.probs.for)
          log.prop.for <- log(recruit.probs.for[cohort.prop])
          if(cohort.prop <= n.year){
            sex.prop <- 0
            z.start.prop <- cohort.prop
          } else {
            sex.prop <- 1
            z.start.prop <- cohort.prop - n.year
          }
          
          #set proposed sex first, then update phi
          model$sex[pick] <<- sex.prop
          model$calculate(phi.nodes[pick.idx2])

          #simulate proposed z vector
          z.prop <- rep(0, n.year)
          z.prop[z.start.prop] <- 1
          if(z.start.prop < n.year){
            for(g in (z.start.prop + 1):n.year){
              z.prop[g] <- rbinom(1, 1, model$phi[pick, g-1] * z.prop[g-1])
              log.prop.for <- log.prop.for +
                dbinom(z.prop[g], 1, model$phi[pick, g-1] * z.prop[g-1], log = TRUE)
            }
          }
          z.stop.prop <- max(which(z.prop == 1))
          
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
          if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year] - model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          #repeat for sex
          if(sex.prop==0){
            model$N.M <<- model$N.M + model$z[pick,]
            if(z.start.prop > 1){ #if wasn't in pop in year 1
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.year] - model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{
            model$N.F <<- model$N.F + model$z[pick,]
            if(z.start.prop > 1){ #if wasn't in pop in year 1
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.year] - model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          #get proposed logprobs for N and y
          model$calculate(ER.M.nodes) #update ER when N updated
          model$calculate(ER.F.nodes) #update ER when N updated
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx])
          lp.proposed.surv <- model$calculate(z.nodes[pick])

          lp.initial.total <- lp.initial.N.M + lp.initial.N.F + lp.initial.y +
            lp.initial.N.recruit.M + lp.initial.N.recruit.F + lp.initial.surv
          lp.proposed.total <- lp.proposed.N.M + lp.proposed.N.F + lp.proposed.y +
            lp.proposed.N.recruit.M + lp.proposed.N.recruit.F + lp.proposed.surv

          #backwards prior and select probs
          #move from z.super==0 cell to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[2*n.year + 1] <- entry.counts.prop[2*n.year + 1] - 1
          entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1

          #p select on guy
          non.back <- sum(model$z.super == 1)
          log.p.select.back <- log(1/non.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.back <- 0

          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
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
            for(g in 1:(n.year-1)){
              mvSaved["phi",1][pick,g] <<- model[["phi"]][pick,g]
            }
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
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
            for(g in 1:(n.year-1)){
              model[["phi"]][pick,g] <<- mvSaved["phi",1][pick,g]
            }
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            #set these logProbs back
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
            model$calculate(y.nodes[pick.idx])
            model$calculate(z.nodes[pick])
          }
        }
      }
    }
    
    #4) Finally, detected guy unobserved sex update
    for(i in 1:length(sex.up)){
      if(z.obs[sex.up[i]]==1){ #only do detected guys here
        i.idx <- seq(sex.up[i],M*n.year,M) #used to reference correct y nodes
        i.idx2 <- seq(sex.up[i],M*(n.year-1),M) #used to reference correct phi nodes
        cohort.curr <- model$z.start[sex.up[i]] + model$sex[sex.up[i]]*n.year

        #get initial logProbs
        lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
        lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
        lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
        lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
        lp.initial.z <- model$getLogProb(z.nodes[sex.up[i]])
        lp.initial.y <- model$getLogProb(y.nodes[i.idx])
        lp.initial.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))
        
        #update N variables
        if(model$sex[sex.up[i]]==0){ #initial male
          #move this guy from N.M to N.F
          model$N.M <<- model$N.M - model$z[sex.up[i],]
          model$N.F <<- model$N.F + model$z[sex.up[i],]
          #move male recruit to female recruit
          if(model$z.start[sex.up[i]]>1){ #otherwise, this is year 1, nothing to change
            model$N.recruit.M[model$z.start[sex.up[i]]-1] <<- model$N.recruit.M[model$z.start[sex.up[i]]-1] - 1
            model$N.recruit.F[model$z.start[sex.up[i]]-1] <<- model$N.recruit.F[model$z.start[sex.up[i]]-1] + 1
          }
          # #update male, female survivors
          model$N.survive.M <<- model$N.M[2:n.year] - model$N.recruit.M
          model$N.survive.F <<- model$N.F[2:n.year] - model$N.recruit.F
        }else{ #initial female
          #move this guy from N.F to N.M
          model$N.F <<- model$N.F - model$z[sex.up[i],]
          model$N.M <<- model$N.M + model$z[sex.up[i],]
          #move female recruit to male recruit
          if(model$z.start[sex.up[i]]>1){ #otherwise, this is year 1, nothing to change
            model$N.recruit.F[model$z.start[sex.up[i]]-1] <<- model$N.recruit.F[model$z.start[sex.up[i]]-1] - 1
            model$N.recruit.M[model$z.start[sex.up[i]]-1] <<- model$N.recruit.M[model$z.start[sex.up[i]]-1] + 1
          }
          # #update male, female survivors
          model$N.survive.F <<- model$N.F[2:n.year] - model$N.recruit.F
          model$N.survive.M <<- model$N.M[2:n.year] - model$N.recruit.M
        }
        model$calculate(ER.M.nodes) #update ER when N.M updated
        model$calculate(ER.F.nodes) #update ER when N.F updated
        #update sex
        model$sex[sex.up[i]] <<- 1 - model$sex[sex.up[i]]
        #update phi nodes when sex changes
        model$calculate(phi.nodes[i.idx2])
        #update prior
        cohort.prop <- model$z.start[sex.up[i]] + model$sex[sex.up[i]]*n.year
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[cohort.curr] <- entry.counts.prop[cohort.curr] - 1
        entry.counts.prop[cohort.prop] <- entry.counts.prop[cohort.prop] + 1
        
        #get proposed logProbs
        lp.proposed.N.M <- model$calculate(N.M.nodes[1])
        lp.proposed.N.F <- model$calculate(N.F.nodes[1])
        lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
        lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
        lp.proposed.z <- model$calculate(z.nodes[sex.up[i]])
        lp.proposed.y <- model$calculate(y.nodes[i.idx])
        lp.proposed.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))

        lp.initial.total <- lp.initial.z + lp.initial.y + lp.initial.N.recruit.M +
          lp.initial.N.recruit.F + lp.initial.N.M + lp.initial.N.F + lp.initial.prior
        lp.proposed.total <- lp.proposed.z + lp.proposed.y + lp.proposed.N.recruit.M +
          lp.proposed.N.recruit.F + lp.proposed.N.M + lp.proposed.N.F + lp.proposed.prior

        #MH step
        log_MH_ratio <- lp.proposed.total - lp.initial.total
        accept <- decide(log_MH_ratio)
        if(accept){
          mvSaved["N.M",1] <<- model[["N.M"]]
          mvSaved["N.F",1] <<- model[["N.F"]]
          mvSaved["ER.M",1] <<- model[["ER.M"]]
          mvSaved["ER.F",1] <<- model[["ER.F"]]
          mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
          mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
          mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
          mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
          mvSaved["sex",1][sex.up[i]] <<- model[["sex"]][sex.up[i]]
          for(g in 1:(n.year-1)){
            mvSaved["phi",1][sex.up[i],g] <<- model[["phi"]][sex.up[i],g]
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["N.M"]] <<- mvSaved["N.M",1]
          model[["N.F"]] <<- mvSaved["N.F",1]
          model[["ER.M"]] <<- mvSaved["ER.M",1]
          model[["ER.F"]] <<- mvSaved["ER.F",1]
          model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
          model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
          model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
          model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
          model[["sex"]][sex.up[i]] <<- mvSaved["sex",1][sex.up[i]]
          for(g in 1:(n.year-1)){
            model[["phi"]][sex.up[i],g] <<- mvSaved["phi",1][sex.up[i],g]
          }
          model$calculate(y.nodes[i.idx])
          model$calculate(phi.nodes[i.idx2])
          model$calculate(z.nodes[sex.up[i]])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
        }
      }
    }

    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)