#Wu-type entry/exit Gibbs updates for Chandler and Clark (2014)
#Because realized N determines recruitment probabilities, individuals are not
#conditionally independent. Conditional on all other individual histories,
#the focal individual's likelihood still factors by primary occasion.
#Detected individuals use separate Type I entry and Type II exit updates.
#For undetected individuals, entry/never-entry is sampled with exit
#marginalized, followed by a Type II exit draw conditional on entry.

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
    B.curr <- rep(0,max(1,n.primary-1))
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
    for(g in 1:(n.primary-1)){
      B.curr[g] <- Acum.curr[g+1] - Acum.curr[g]
    }
    
    for(i in 1:M){
      z.curr <- rep(0,n.primary)
      a.curr <- rep(0,n.primary)
      z.new <- rep(0,n.primary)
      a.new <- rep(0,n.primary)
      N.minus <- rep(0,n.primary)
      Acum.minus <- rep(0,n.primary)
      log.y.abs <- rep(0,n.primary)
      log.y.alive <- rep(0,n.primary)
      r.U <- rep(0,n.primary)
      logH.U <- rep(0,n.primary)
      logH.A <- rep(0,n.primary)
      logH.D <- rep(0,n.primary)
      logV.A <- rep(0,n.primary)
      logV.D <- rep(0,n.primary)
      logUreach <- rep(-Inf,n.primary)
      logAliveToFirst <- rep(0,n.primary)
      logw <- rep(-Inf,n.primary+1)
      probs <- rep(0,n.primary+1)
      
      #current focal history and leave-one-out population counts
      entered <- 0
      detected <- 0
      first.det <- n.primary+1
      last.det <- 0
      for(g in 1:n.primary){
        z.curr[g] <- model$z[i,g]
        if(z.curr[g] == 1){
          entered <- 1
        }
        a.curr[g] <- entered
        N.minus[g] <- N.curr[g] - z.curr[g]
        Acum.minus[g] <- Acum.curr[g] - a.curr[g]
        if(model$y[i,g] == 0){
          log.y.abs[g] <- 0
        }else{
          log.y.abs[g] <- -Inf
          detected <- 1
          if(first.det == n.primary+1){
            first.det <- g
          }
          last.det <- g
        }
        log.y.alive[g] <- dbinom(model$y[i,g],size=K[g],p=model$p[g],log=TRUE)
      }
      
      #For each interval, calculate the recruitment probability implied by the
      #focal individual being not yet entered, and the recruitment likelihood
      #of all OTHER individuals under each possible focal state at occasion g:
      #U = not yet entered, A = alive, D = entered previously and dead.
      for(g in 1:(n.primary-1)){
        recruit.curr <- a.curr[g+1] - a.curr[g]
        R.minus <- B.curr[g] - recruit.curr
        U.minus <- (M-1) - Acum.minus[g]
        stayU.minus <- U.minus - R.minus
        
        #U: focal individual has not yet entered
        A.raw.cand <- M - Acum.minus[g]
        if(A.raw.cand < 0.01){
          A.cand <- 0.01
        }else{
          A.cand <- A.raw.cand
        }
        r.cand <- N.minus[g]*model$gamma[g]/A.cand
        if(r.cand > 0.999){
          r.cand <- 0.999
        }
        r.U[g] <- r.cand
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
        logH.U[g] <- lp
        
        #A: focal individual has entered and is alive
        A.raw.cand <- M - (Acum.minus[g]+1)
        if(A.raw.cand < 0.01){
          A.cand <- 0.01
        }else{
          A.cand <- A.raw.cand
        }
        r.cand <- (N.minus[g]+1)*model$gamma[g]/A.cand
        if(r.cand > 0.999){
          r.cand <- 0.999
        }
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
        logH.A[g] <- lp
        
        #D: focal individual has entered previously and is dead
        r.cand <- N.minus[g]*model$gamma[g]/A.cand
        if(r.cand > 0.999){
          r.cand <- 0.999
        }
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
        logH.D[g] <- lp
      }
      
      #Type II backward quantities. These integrate over future exit when
      #sampling entry for an undetected individual, and are then reused to draw
      #exit conditional on the sampled entry occasion.
      logV.A[n.primary] <- 0
      logV.D[n.primary] <- 0
      for(k in 1:(n.primary-1)){
        g <- n.primary-k
        logV.D[g] <- logH.D[g] + log.y.abs[g+1] + logV.D[g+1]
        x1 <- log(model$phi[i]) + log.y.alive[g+1] + logV.A[g+1]
        x2 <- log(1-model$phi[i]) + log.y.abs[g+1] + logV.D[g+1]
        if(x1 > x2){
          maxlp <- x1
        }else{
          maxlp <- x2
        }
        if(maxlp == -Inf){
          logV.A[g] <- -Inf
        }else{
          logV.A[g] <- logH.A[g] + maxlp + log(exp(x1-maxlp)+exp(x2-maxlp))
        }
      }
      
      #probability of reaching each occasion without having entered
      logUreach[1] <- log(1-model$psi[1]) + log.y.abs[1]
      for(g in 1:(n.primary-1)){
        logUreach[g+1] <- logUreach[g] + logH.U[g] + log(1-r.U[g]) + log.y.abs[g+1]
      }
      
      if(detected == 1){
        #Type I: sample entry occasion conditional on being alive at first detection.
        logAliveToFirst[first.det] <- 0
        if(first.det > 1){
          for(k in 1:(first.det-1)){
            g <- first.det-k
            logAliveToFirst[g] <- logH.A[g] + log(model$phi[i]) +
              log.y.alive[g+1] + logAliveToFirst[g+1]
          }
        }
        logw[1] <- log(model$psi[1]) + log.y.alive[1] + logAliveToFirst[1]
        if(first.det > 1){
          for(e in 2:first.det){
            if(r.U[e-1] > 0){
              logw[e] <- logUreach[e-1] + logH.U[e-1] + log(r.U[e-1]) +
                log.y.alive[e] + logAliveToFirst[e]
            }else{
              logw[e] <- -Inf
            }
          }
        }
        maxlp <- max(logw[1:first.det])
        for(e in 1:first.det){
          probs[e] <- exp(logw[e]-maxlp)
        }
        p.denom <- sum(probs[1:first.det])
        for(e in 1:first.det){
          probs[e] <- probs[e]/p.denom
        }
        entry.new <- rcat(1,probs[1:first.det])
        #alive from entry through last detection
        for(g in 1:n.primary){
          z.new[g] <- 0
        }
        for(g in entry.new:last.det){
          z.new[g] <- 1
        }
        
        #Type II: sample exit after the last detection, conditional on entry.
        alive <- 1
        if(last.det < n.primary){
          for(g in last.det:(n.primary-1)){
            if(alive == 1){
              x1 <- log(model$phi[i]) + log.y.alive[g+1] + logV.A[g+1]
              x2 <- log(1-model$phi[i]) + log.y.abs[g+1] + logV.D[g+1]
              if(x1 > x2){
                maxlp <- x1
              }else{
                maxlp <- x2
              }
              probs[1] <- exp(x1-maxlp)
              probs[2] <- exp(x2-maxlp)
              p.denom <- probs[1]+probs[2]
              probs[1] <- probs[1]/p.denom
              probs[2] <- probs[2]/p.denom
              survive <- rcat(1,probs[1:2])
              if(survive == 1){
                z.new[g+1] <- 1
              }else{
                z.new[g+1] <- 0
                alive <- 0
              }
            }else{
              z.new[g+1] <- 0
            }
          }
        }
      }else{
        #Undetected individual: sample never-entry or entry occasion with future
        #exit marginalized, as in the Wu-type undetected-history update.
        logw[1] <- log(model$psi[1]) + log.y.alive[1] + logV.A[1]
        for(e in 2:n.primary){
          if(r.U[e-1] > 0){
            logw[e] <- logUreach[e-1] + logH.U[e-1] + log(r.U[e-1]) +
              log.y.alive[e] + logV.A[e]
          }else{
            logw[e] <- -Inf
          }
        }
        logw[n.primary+1] <- logUreach[n.primary] #never enters
      
        maxlp <- max(logw)
        for(e in 1:(n.primary+1)){
          probs[e] <- exp(logw[e]-maxlp)
        }
        probs <- probs/sum(probs)
        entry.new <- rcat(1,probs)
        
        for(g in 1:n.primary){
          z.new[g] <- 0
        }
        if(entry.new <= n.primary){
          z.new[entry.new] <- 1
          
          #Type II: sample exit conditional on the new entry occasion.
          alive <- 1
          if(entry.new < n.primary){
            for(g in entry.new:(n.primary-1)){
              if(alive == 1){
                x1 <- log(model$phi[i]) + log.y.alive[g+1] + logV.A[g+1]
                x2 <- log(1-model$phi[i]) + log.y.abs[g+1] + logV.D[g+1]
                if(x1 > x2){
                  maxlp <- x1
                }else{
                  maxlp <- x2
                }
                probs[1] <- exp(x1-maxlp)
                probs[2] <- exp(x2-maxlp)
                p.denom <- probs[1]+probs[2]
                probs[1] <- probs[1]/p.denom
                probs[2] <- probs[2]/p.denom
                survive <- rcat(1,probs[1:2])
                if(survive == 1){
                  z.new[g+1] <- 1
                }else{
                  z.new[g+1] <- 0
                  alive <- 0
                }
              }else{
                z.new[g+1] <- 0
              }
            }
          }
        }
      }
      
      #replace focal history and update aggregate counts for next individual
      entered <- 0
      for(g in 1:n.primary){
        if(z.new[g] == 1){
          entered <- 1
        }
        a.new[g] <- entered
        model$z[i,g] <<- z.new[g]
        N.curr[g] <- N.curr[g] - z.curr[g] + z.new[g]
        Acum.curr[g] <- Acum.curr[g] - a.curr[g] + a.new[g]
      }
      for(g in 1:(n.primary-1)){
        B.curr[g] <- Acum.curr[g+1] - Acum.curr[g]
      }
    }
    
    #update deterministic nodes and log probabilities once after the full sweep
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)