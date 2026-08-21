zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    M <- control$M
    K <- control$K
    n.primary <- control$n.primary
    z.obs <- control$z.obs
    first.det <- control$first.det
    last.det <- control$last.det
    calcNodes <- model$getDependencies(target)
  },
  run = function(){
    #probability of zero detections in each primary occasion if alive
    q <- rep(0,n.primary)
    for(g in 1:n.primary){
      q[g] <- (1-model$p[g])^K[g]
    }
    #unconditional probabilities of first entry in each primary occasion
    #and of never entering
    beta <- rep(0,n.primary)
    beta[1] <- model$psi[1]
    remaining.entry <- 1-model$psi[1]
    
    for(g in 2:n.primary){
      beta[g] <- remaining.entry*model$gamma[g-1]
      remaining.entry <- remaining.entry*(1-model$gamma[g-1])
    }
    beta.never <- remaining.entry
    lam <- rep(0,n.primary)
    v <- rep(0,n.primary)
    entry.probs <- rep(0,n.primary)
    trajectory.probs <- rep(0,n.primary+1)
    for(i in 1:M){
      #Wu Type I recursion:
      #lam[g] = probability of having entered, still being alive,
      #and having no detections before occasion g
      lam[1] <- beta[1]
      for(g in 1:(n.primary-1)){
        lam[g+1] <- beta[g+1] + lam[g]*q[g]*model$phi[i]
      }
      
      #Wu Type II recursion:
      #v[g] = probability of no future detections after occasion g,
      #allowing either death before g+1 or survival with no later detection
      v[n.primary] <- 1
      #nimble does not allow decreasing numbers in loops
      for(k in 1:(n.primary-1)){
        g <- n.primary-k
        v[g] <- (1-model$phi[i]) + model$phi[i]*q[g+1]*v[g+1]
      }
      #overwrite complete trajectory
      for(g in 1:n.primary){
        model$z[i,g] <<- 0
      }
      #1) detected individuals
      if(z.obs[i]==1){
        f <- first.det[i]
        l <- last.det[i]
        #1a) Wu Type I block:
        #sample entry occasion conditional on first detection at f
        for(g in 1:n.primary){
          entry.probs[g] <- 0
        }
        remaining <- 1
        for(k in 1:f){
          g <- f-k+1
          cond.entry <- beta[g]/lam[g]
          entry.probs[g] <- remaining*cond.entry
          remaining <- remaining*(1-cond.entry)
        }
        entry.probs <- entry.probs/sum(entry.probs)
        e.curr <- rcat(1,entry.probs)
        #known alive from sampled entry through last detection
        for(g in e.curr:l){
          model$z[i,g] <<- 1
        }
        #1b) Wu Type II block:
        #sample exit after last detection conditional on
        #no subsequent detections
        alive <- 1
        if(l < n.primary){
          for(g in (l+1):n.primary){
            if(alive==1){
              death.prob <- (1-model$phi[i])/v[g-1]
              if(rbinom(1,1,death.prob)==1){
                alive <- 0
              }else{
                model$z[i,g] <<- 1
              }
            }
          }
        }
      }else{
        #2) undetected individuals:
        #sample either never entering or the complete entry/survival
        #trajectory conditional on the all-zero capture history
        trajectory.probs[1] <- beta.never
        rho <- beta.never
        for(g in 1:n.primary){
          trajectory.probs[g+1] <- beta[g]*q[g]*v[g]
          rho <- rho+trajectory.probs[g+1]
        }
        for(g in 1:(n.primary+1)){
          trajectory.probs[g] <- trajectory.probs[g]/rho
        }
        trajectory.curr <- rcat(1,trajectory.probs)
        if(trajectory.curr > 1){
          e.curr <- trajectory.curr-1
          model$z[i,e.curr] <<- 1
          alive <- 1
          #Wu Type II block conditional on no detections
          if(e.curr < n.primary){
            for(g in (e.curr+1):n.primary){
              if(alive==1){
                death.prob <- (1-model$phi[i])/v[g-1]
                if(rbinom(1,1,death.prob)==1){
                  alive <- 0
                }else{
                  model$z[i,g] <<- 1
                }
              }
            }
          }
        }
      }
    }
    #recalculate z log probabilities and all downstream nodes
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)

gammaGibbsSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    calcNodes <- model$getDependencies(target)
    M <- control$M
    n.primary <- control$n.primary
    
    #get Beta prior parameters directly from model
    shape1 <- model$getParam(target,"shape1")
    shape2 <- model$getParam(target,"shape2")
    
    if(target=="gamma"){
      is.fixed.gamma <- TRUE
      g <- 1
    }else{
      is.fixed.gamma <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
    }
  },
  run = function(){
    if(is.fixed.gamma){
      count <- 0
      failures <- 0
      for(j in 1:(n.primary-1)){
        recruits <- model$B[j]
        available <- M-model$Acum[j]
        count <- count+recruits
        failures <- failures+available-recruits
      }
    }else{
      count <- model$B[g]
      available <- M-model$Acum[g]
      failures <- available-count
    }
    
    model[[target]] <<- rbeta(1,
                              shape1+count,
                              shape2+failures)
    
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)