#Wu-type whole-history Gibbs update for Chandler and Clark (2014)
#Because realized N determines recruitment probabilities, individuals are not
#conditionally independent. Conditional on all other individual histories,
#however, the focal individual's full conditional is a 3-state HMM:
#1 = not yet entered, 2 = alive, 3 = entered previously and dead.
#This sampler uses exact forward-filtering/backward-sampling in O(n.primary)
#for each individual, while accounting for the focal individual's effect on
#the recruitment probabilities of all other individuals.

zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    M <- control$M
    K <- control$K
    n.primary <- control$n.primary
    calcNodes <- model$getDependencies(target)
  },
  run = function(){
    #current abundance and cumulative-entry counts
    N.curr <- rep(0,n.primary)
    Acum.curr <- rep(0,n.primary)
    if(n.primary > 1){
      B.curr <- rep(0,n.primary-1)
    }else{
      B.curr <- rep(0,1)
    }
    for(i in 1:M){
      entered <- 0
      for(g in 1:n.primary){
        N.curr[g] <- N.curr[g] + model$z[i,g]
        if(model$z[i,g] == 1){
          entered <- 1
        }
        Acum.curr[g] <- Acum.curr[g] + entered
      }
    }
    if(n.primary > 1){
      for(g in 1:(n.primary-1)){
        B.curr[g] <- Acum.curr[g+1] - Acum.curr[g]
      }
    }

    #update complete history for each individual
    for(i in 1:M){
      z.curr <- rep(0,n.primary)
      a.curr <- rep(0,n.primary)
      z.new <- rep(0,n.primary)
      a.new <- rep(0,n.primary)
      N.minus <- rep(0,n.primary)
      Acum.minus <- rep(0,n.primary)
      state <- rep(0,n.primary)
      alpha <- rep(-Inf,3*n.primary)
      probs <- rep(0,3)
      logw <- rep(-Inf,3)

      #leave-one-out population counts
      entered <- 0
      for(g in 1:n.primary){
        z.curr[g] <- model$z[i,g]
        if(z.curr[g] == 1){
          entered <- 1
        }
        a.curr[g] <- entered
        N.minus[g] <- N.curr[g] - z.curr[g]
        Acum.minus[g] <- Acum.curr[g] - a.curr[g]
      }

      #for each interval and each possible focal state, calculate
      #1) gamma.prime implied by that state and
      #2) recruitment likelihood of all OTHER individuals
      if(n.primary > 1){
        r.state <- rep(0,3*(n.primary-1))
        logH.state <- rep(0,3*(n.primary-1))
        for(g in 1:(n.primary-1)){
          recruit.curr <- a.curr[g+1] - a.curr[g]
          R.minus <- B.curr[g] - recruit.curr
          U.minus <- (M-1) - Acum.minus[g]
          stayU.minus <- U.minus - R.minus

          for(s in 1:3){
            alive.cand <- 0
            entered.cand <- 0
            if(s == 2){
              alive.cand <- 1
              entered.cand <- 1
            }
            if(s == 3){
              entered.cand <- 1
            }

            N.cand <- N.minus[g] + alive.cand
            Acum.cand <- Acum.minus[g] + entered.cand
            A.raw.cand <- M - Acum.cand
            if(A.raw.cand < 0.01){
              A.cand <- 0.01
            }else{
              A.cand <- A.raw.cand
            }
            r.cand <- N.cand*model$gamma[g]/A.cand
            if(r.cand > 0.999){
              r.cand <- 0.999
            }

            idx <- (g-1)*3+s
            r.state[idx] <- r.cand
            lp <- 0
            if(R.minus > 0){
              if(r.cand > 0){
                lp <- lp + R.minus*log(r.cand)
              }else{
                lp <- -Inf
              }
            }
            if(stayU.minus > 0 & lp > -Inf){
              lp <- lp + stayU.minus*log(1-r.cand)
            }
            logH.state[idx] <- lp
          }
        }
      }else{
        r.state <- rep(0,3)
        logH.state <- rep(0,3)
      }

      #forward filtering
      #state 1 = not yet entered
      #state 2 = alive
      #state 3 = entered previously and dead
      if(model$y[i,1] == 0){
        log.y.abs <- 0
      }else{
        log.y.abs <- -Inf
      }
      log.y.alive <- dbinom(model$y[i,1],size=K[1],p=model$p[1],log=TRUE)
      alpha[1] <- log(1-model$psi[1]) + log.y.abs
      alpha[2] <- log(model$psi[1]) + log.y.alive
      alpha[3] <- -Inf

      if(n.primary > 1){
        for(g in 1:(n.primary-1)){
          idx.U <- (g-1)*3+1
          idx.A <- (g-1)*3+2
          idx.D <- (g-1)*3+3
          idx.next.U <- g*3+1
          idx.next.A <- g*3+2
          idx.next.D <- g*3+3

          r.U <- r.state[idx.U]
          logH.U <- logH.state[idx.U]
          logH.A <- logH.state[idx.A]
          logH.D <- logH.state[idx.D]

          if(model$y[i,g+1] == 0){
            log.y.abs <- 0
          }else{
            log.y.abs <- -Inf
          }
          log.y.alive <- dbinom(model$y[i,g+1],size=K[g+1],p=model$p[g+1],log=TRUE)

          #not yet entered -> not yet entered
          alpha[idx.next.U] <- log.y.abs + alpha[idx.U] + logH.U + log(1-r.U)

          #not yet entered -> alive OR alive -> alive
          if(r.U > 0){
            x1 <- alpha[idx.U] + logH.U + log(r.U)
          }else{
            x1 <- -Inf
          }
          x2 <- alpha[idx.A] + logH.A + log(model$phi[i])
          if(x1 > x2){
            maxlp <- x1
          }else{
            maxlp <- x2
          }
          if(maxlp == -Inf){
            alpha[idx.next.A] <- -Inf
          }else{
            alpha[idx.next.A] <- log.y.alive + maxlp +
              log(exp(x1-maxlp)+exp(x2-maxlp))
          }

          #alive -> dead OR dead -> dead
          x1 <- alpha[idx.A] + logH.A + log(1-model$phi[i])
          x2 <- alpha[idx.D] + logH.D
          if(x1 > x2){
            maxlp <- x1
          }else{
            maxlp <- x2
          }
          if(maxlp == -Inf){
            alpha[idx.next.D] <- -Inf
          }else{
            alpha[idx.next.D] <- log.y.abs + maxlp +
              log(exp(x1-maxlp)+exp(x2-maxlp))
          }
        }
      }

      #sample state in final primary occasion
      for(s in 1:3){
        logw[s] <- alpha[(n.primary-1)*3+s]
      }
      maxlp <- max(logw)
      for(s in 1:3){
        probs[s] <- exp(logw[s]-maxlp)
      }
      probs <- probs/sum(probs)
      state[n.primary] <- rcat(1,probs)

      #backward sampling
      if(n.primary > 1){
        #nimble does not allow decreasing numbers in loops
        for(k in 1:(n.primary-1)){
          g <- n.primary-k
          next.state <- state[g+1]
          for(s in 1:3){
            logw[s] <- -Inf
          }

          #previous state = not yet entered
          idx <- (g-1)*3+1
          r.U <- r.state[idx]
          if(next.state == 1){
            logw[1] <- alpha[idx] + logH.state[idx] + log(1-r.U)
          }
          if(next.state == 2 & r.U > 0){
            logw[1] <- alpha[idx] + logH.state[idx] + log(r.U)
          }

          #previous state = alive
          idx <- (g-1)*3+2
          if(next.state == 2){
            logw[2] <- alpha[idx] + logH.state[idx] + log(model$phi[i])
          }
          if(next.state == 3){
            logw[2] <- alpha[idx] + logH.state[idx] + log(1-model$phi[i])
          }

          #previous state = entered previously and dead
          idx <- (g-1)*3+3
          if(next.state == 3){
            logw[3] <- alpha[idx] + logH.state[idx]
          }

          maxlp <- max(logw)
          for(s in 1:3){
            probs[s] <- exp(logw[s]-maxlp)
          }
          probs <- probs/sum(probs)
          state[g] <- rcat(1,probs)
        }
      }

      #replace focal history and update aggregate counts for next individual
      for(g in 1:n.primary){
        if(state[g] == 2){
          z.new[g] <- 1
        }else{
          z.new[g] <- 0
        }
        if(state[g] == 1){
          a.new[g] <- 0
        }else{
          a.new[g] <- 1
        }
        model$z[i,g] <<- z.new[g]
        N.curr[g] <- N.curr[g] - z.curr[g] + z.new[g]
        Acum.curr[g] <- Acum.curr[g] - a.curr[g] + a.new[g]
      }
      if(n.primary > 1){
        for(g in 1:(n.primary-1)){
          B.curr[g] <- Acum.curr[g+1] - Acum.curr[g]
        }
      }
    }

    #update deterministic nodes and log probabilities once after the full sweep
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)
