eSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    K <- control$K
    n.primary <- control$n.primary
    z.obs <- control$z.obs
    first.det <- control$first.det
    last.det <- control$last.det
    calcNodes <- control$calcNodes
  },
  run = function(){
    #probability of zero detections in each primary occasion if alive
    q <- rep(0,n.primary)
    for(g in 1:n.primary){
      q[g] <- (1-model$p[g])^K[g]
    }
    lam <- rep(0,n.primary)
    v <- rep(0,n.primary)
    prop.probs <- rep(0,n.primary)
    for(i in 1:M){
      #Wu Type I recursion:
      #lam[g] = probability of having entered, still being alive,
      #and having no detections before occasion g
      lam[1] <- model$pi[1]
      for(g in 1:(n.primary-1)){
        lam[g+1] <- model$pi[g+1] + lam[g]*q[g]*model$phi[i]
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
      #1) detected individuals
      if(z.obs[i]==1){
        model$z.super[i] <<- 1
        f <- first.det[i]
        l <- last.det[i]
        #1a) Wu Type I block:
        #sample entry occasion conditional on first detection at f
        for(g in 1:n.primary){
          prop.probs[g] <- 0
        }
        remaining <- 1
        for(k in 1:f){
          g <- f-k+1
          cond.entry <- model$pi[g]/lam[g]
          prop.probs[g] <- remaining*cond.entry
          remaining <- remaining*(1-cond.entry)
        }
        prop.probs <- prop.probs/sum(prop.probs)
        e.curr <- rcat(1,prop.probs)
        #1b) Wu Type II block:
        #sample exit occasion conditional on no detections after last detection
        d.curr <- l
        if(l < n.primary){
          for(g in (l+1):n.primary){
            if(d.curr == g-1){
              death.prob <- (1-model$phi[i])/v[g-1]
              if(rbinom(1,1,death.prob) == 0){
                d.curr <- g
              }
            }
          }
        }
      }else{
        #2) undetected individuals:
        #update z.super with entry and survival trajectory marginalized
        rho <- 0
        for(g in 1:n.primary){
          prop.probs[g] <- model$pi[g]*q[g]*v[g]
          rho <- rho+prop.probs[g]
        }
        z.super.prob <- model$psi[1]*rho/((1-model$psi[1])+model$psi[1]*rho)
        model$z.super[i] <<- rbinom(1,1,z.super.prob)
        if(model$z.super[i] == 0){
          #2a) not in superpopulation: trajectory is independent of data,
          #so draw directly from its full conditional (prior)
          e.curr <- rcat(1,model$pi[1:n.primary])
          d.curr <- e.curr
          if(e.curr < n.primary){
            for(g in (e.curr+1):n.primary){
              if(d.curr == g-1){
                if(rbinom(1,1,model$phi[i]) == 1){
                  d.curr <- g
                }
              }
            }
          }
        }else{
          #2b) in superpopulation but undetected:
          #exact Wu Gibbs update of entry occasion and survival/exit trajectory
          for(g in 1:n.primary){
            prop.probs[g] <- prop.probs[g]/rho
          }
          e.curr <- rcat(1,prop.probs)
          d.curr <- e.curr
          #sample survival/exit conditional on the all-zero capture history
          if(e.curr < n.primary){
            for(g in (e.curr+1):n.primary){
              if(d.curr == g-1){
                death.prob <- (1-model$phi[i])/v[g-1]
                if(rbinom(1,1,death.prob) == 0){
                  d.curr <- g
                }
              }
            }
          }
        }
      }
      #set new entry occasion and survival trajectory
      model$e[i] <<- e.curr
      for(g in 2:n.primary){
        if(g > e.curr & g <= d.curr){
          model$surv[i,g] <<- 1
        }else{
          model$surv[i,g] <<- 0
        }
      }
    }
    #update deterministic nodes and log probabilities
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function () {})
)
