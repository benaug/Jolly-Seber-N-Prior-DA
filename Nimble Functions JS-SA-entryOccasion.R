eSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    K <- control$K
    n.primary <- control$n.primary
    z.obs <- control$z.obs
    first.det <- control$first.det
    last.det <- control$last.det
    e.nodes <- control$e.nodes
    surv.p.nodes <- control$surv.p.nodes
    surv.nodes <- control$surv.nodes
    z.nodes <- control$z.nodes
    y.size.nodes <- control$y.size.nodes
    y.nodes <- control$y.nodes
    N.B.nodes <- control$N.B.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    for(i in 1:M){
      y.idx <- seq(i,M*n.primary,M)
      z.idx <- y.idx #z and y.size are M x n.primary, same layout as y
      surv.idx <- seq(i,M*(n.primary-1),M) #surv and surv.p only have nodes for g=2:n.primary
      e.curr <- model$e[i] #current entry
      #current exit read from surv, not z. surv[i,g]==1 only for e < g <= d,
      #so the last 1 in the row is d. Does not depend on z being up to date.
      d.curr <- e.curr
      for(g in 2:n.primary){
        if(model$surv[i,g]==1){
          d.curr <- g
        }
      }
      #1) for detected guys, update entry and exit occasions considering when they were detected and known to be in population
      if(z.obs[i]==1){
        #1a) update entry occasion: categorical Gibbs update for values 1:first.det
        if(first.det[i] > 1){
          lp.e <- rep(-Inf,n.primary)
          for(e.prop in 1:first.det[i]){
            model$e[i] <<- e.prop
            for(g in 2:n.primary){
              if(g > e.prop & g <= d.curr){
                model$surv[i,g] <<- 1
              }else{
                model$surv[i,g] <<- 0
              }
            }
            model$calculate(z.nodes[z.idx]) #z is deterministic, let the model build it
            model$calculate(surv.p.nodes[surv.idx]) #surv.p is the parent of surv, must follow z
            model$calculate(y.size.nodes[y.idx]) #y.size is the parent of y, must follow z
            lp.e[e.prop] <- model$calculate(e.nodes[i]) +
              model$calculate(surv.nodes[surv.idx]) +
              model$calculate(y.nodes[y.idx])
          }
          maxlp <- max(lp.e)
          prop.probs <- exp(lp.e - maxlp)
          prop.probs <- prop.probs/sum(prop.probs)
          
          e.curr <- rcat(1,prop.probs)
          model$e[i] <<- e.curr
          for(g in 2:n.primary){
            if(g > e.curr & g <= d.curr){
              model$surv[i,g] <<- 1
            }else{
              model$surv[i,g] <<- 0
            }
          }
          model$calculate(z.nodes[z.idx])
          model$calculate(surv.p.nodes[surv.idx])
          model$calculate(y.size.nodes[y.idx])
          model$calculate(e.nodes[i])
          model$calculate(surv.nodes[surv.idx])
          model$calculate(y.nodes[y.idx])
        }
        
        #1b) exit occasion: categorical Gibbs update for values last.det:n.primary,
        if(last.det[i] < n.primary){
          lp.d <- rep(-Inf,n.primary)
          for(d.prop in last.det[i]:n.primary){
            for(g in 2:n.primary){
              if(g > e.curr & g <= d.prop){
                model$surv[i,g] <<- 1
              }else{
                model$surv[i,g] <<- 0
              }
            }
            model$calculate(z.nodes[z.idx])
            model$calculate(surv.p.nodes[surv.idx])
            model$calculate(y.size.nodes[y.idx])
            #entry prior unchanged, omit from the Gibbs weights
            lp.d[d.prop] <- model$calculate(surv.nodes[surv.idx]) +
              model$calculate(y.nodes[y.idx])
          }
          maxlp <- max(lp.d)
          prop.probs <- exp(lp.d - maxlp)
          prop.probs <- prop.probs/sum(prop.probs)
          d.curr <- rcat(1,prop.probs)
          
          for(g in 2:n.primary){
            if(g > e.curr & g <= d.curr){
              model$surv[i,g] <<- 1
            }else{
              model$surv[i,g] <<- 0
            }
          }
          model$calculate(z.nodes[z.idx])
          model$calculate(surv.p.nodes[surv.idx])
          model$calculate(y.size.nodes[y.idx])
          model$calculate(e.nodes[i]) #1a may have been skipped, keep this current
          model$calculate(surv.nodes[surv.idx])
          model$calculate(y.nodes[y.idx])
        }
      }else{
        #2) undetected individuals
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
          model$e[i] <<- e.curr
          for(g in 2:n.primary){
            if(g > e.curr & g <= d.curr){
              model$surv[i,g] <<- 1
            }else{
              model$surv[i,g] <<- 0
            }
          }
          model$calculate(z.nodes[z.idx])
          model$calculate(surv.p.nodes[surv.idx])
          model$calculate(e.nodes[i])
          model$calculate(surv.nodes[surv.idx])
          
        }else{
          #2b) in superpopulation but undetected:
          #2b option 1) MH proposal of full trajectory from prior
          # lp.y.curr <- model$getLogProb(y.nodes[y.idx])
          # e.prop <- rcat(1,model$pi[1:n.primary])
          # d.prop <- e.prop
          # if(e.prop < n.primary){
          #   for(g in (e.prop+1):n.primary){
          #     if(d.prop == g-1){
          #       if(rbinom(1,1,model$phi[i]) == 1){
          #         d.prop <- g
          #       }
          #     }
          #   }
          # }
          # model$e[i] <<- e.prop
          # for(g in 2:n.primary){
          #   if(g > e.prop & g <= d.prop){
          #     model$surv[i,g] <<- 1
          #   }else{
          #     model$surv[i,g] <<- 0
          #   }
          # }
          # model$calculate(z.nodes[z.idx])
          # model$calculate(surv.p.nodes[surv.idx])
          # model$calculate(y.size.nodes[y.idx])
          # model$calculate(e.nodes[i])
          # model$calculate(surv.nodes[surv.idx])
          #
          # lp.y.prop <- model$calculate(y.nodes[y.idx])
          # #proposed from prior, so prior and proposal cancel
          # log_MH_ratio <- lp.y.prop - lp.y.curr
          # if(!decide(log_MH_ratio)){
          #   model$e[i] <<- e.curr
          #   for(g in 2:n.primary){
          #     if(g > e.curr & g <= d.curr){
          #       model$surv[i,g] <<- 1
          #     }else{
          #       model$surv[i,g] <<- 0
          #     }
          #   }
          #   model$calculate(z.nodes[z.idx])
          #   model$calculate(surv.p.nodes[surv.idx])
          #   model$calculate(y.size.nodes[y.idx])
          #   model$calculate(e.nodes[i])
          #   model$calculate(surv.nodes[surv.idx])
          #   model$calculate(y.nodes[y.idx])
          # }
          
          #2b) in superpopulation but undetected:
          #2b option 2) exact Gibbs update of entry occasion and survival/exit trajectory
          #probability of zero detections in each primary occasion if alive
          #same algorithm as used in Wu et al. 2021
          q <- rep(0,n.primary)
          for(g in 1:n.primary){
            q[g] <- (1 - model$p[g])^K[g]
          }
          
          #chi[g] = probability of no detections from g:n.primary,
          #conditional on being alive in occasion g
          chi <- rep(0,n.primary)
          chi[n.primary] <- q[n.primary]
          if(n.primary > 1){
            #nimble does not allow decreasing numbers in loops
            for(k in 1:(n.primary-1)){
              g <- n.primary - k
              chi[g] <- q[g] * ((1 - model$phi[i]) +
                                  model$phi[i] * chi[g+1])
            }
          }
          
          #full conditional for entry occasion:
          #P(e=g | y=0, z.super=1, ...) proportional to pi[g] * chi[g]
          prop.probs <- rep(0,n.primary)
          for(g in 1:n.primary){
            prop.probs[g] <- model$pi[g] * chi[g]
          }
          prop.probs <- prop.probs/sum(prop.probs)
          
          e.curr <- rcat(1,prop.probs)
          d.curr <- e.curr
          
          #sample survival/exit conditional on the all-zero capture history
          if(e.curr < n.primary){
            for(g in (e.curr+1):n.primary){
              if(d.curr == g-1){
                
                #P(survive into g | alive at g-1, no detections from g onward)
                surv.prob <- model$phi[i] * chi[g] /
                  ((1 - model$phi[i]) + model$phi[i] * chi[g])
                
                if(rbinom(1,1,surv.prob) == 1){
                  d.curr <- g
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
          
          #update deterministic nodes and log probabilities
          model$calculate(z.nodes[z.idx])
          model$calculate(surv.p.nodes[surv.idx])
          model$calculate(y.size.nodes[y.idx])
          model$calculate(e.nodes[i])
          model$calculate(surv.nodes[surv.idx])
          model$calculate(y.nodes[y.idx])
        }
      }
    }
    model$calculate(N.B.nodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function () {})
)