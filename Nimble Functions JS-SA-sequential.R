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
    lam <- rep(0,n.primary)
    v <- rep(0,n.primary)
    entry.probs <- rep(0,n.primary)
    entry.zero.probs <- rep(0,n.primary)
    for(i in 1:M){
      #Wu Type I recursion:
      #lam[g] = probability of having entered, still being alive,
      #and having no detections before occasion g
      lam[1] <- model$pi[1]
      if(n.primary > 1){
        for(g in 1:(n.primary-1)){
          lam[g+1] <- model$pi[g+1] + lam[g]*q[g]*model$phi[i]
        }
      }
      #Wu Type II recursion:
      #v[g] = probability of no future detections after occasion g,
      #allowing either death before g+1 or survival with no later detection
      v[n.primary] <- 1
      if(n.primary > 1){
        #nimble does not allow decreasing numbers in loops
        for(k in 1:(n.primary-1)){
          g <- n.primary-k
          v[g] <- (1-model$phi[i]) + model$phi[i]*q[g+1]*v[g+1]
        }
      }
      #overwrite complete trajectory
      for(g in 1:n.primary){
        model$z[i,g] <<- 0
      }
      #1) detected individuals
      if(z.obs[i]==1){
        model$z.super[i] <<- 1
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
          cond.entry <- model$pi[g]/lam[g]
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
        #sample departure after last detection conditional on
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
        #2) undetected individuals
        #entry probabilities for an undetected individual that is in
        #the superpopulation:
        #P(e=g | y=0,z.super=1) proportional to pi[g]*q[g]*v[g]
        rho <- 0
        for(g in 1:n.primary){
          entry.zero.probs[g] <- model$pi[g]*q[g]*v[g]
          rho <- rho+entry.zero.probs[g]
        }
        for(g in 1:n.primary){
          entry.zero.probs[g] <- entry.zero.probs[g]/rho
        }
        #update z.super with z trajectory marginalized
        z.super.prob <- model$psi[1]*rho/((1-model$psi[1])+model$psi[1]*rho)
        model$z.super[i] <<- rbinom(1,1,z.super.prob)
        if(model$z.super[i]==0){
          #2a) not in superpopulation, draw from prior
          e.curr <- rcat(1,model$pi[1:n.primary])
          model$z[i,e.curr] <<- 1
          alive <- 1
          if(e.curr < n.primary){
            for(g in (e.curr+1):n.primary){
              if(alive==1){
                if(rbinom(1,1,model$phi[i])==1){
                  model$z[i,g] <<- 1
                }else{
                  alive <- 0
                }
              }
            }
          }
        }else{
          #2b) in superpopulation but never detected
          #Wu Type I block for an all-zero capture history
          e.curr <- rcat(1,entry.zero.probs)
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
    #recalculate z.super and z log probabilities and all downstream nodes
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)