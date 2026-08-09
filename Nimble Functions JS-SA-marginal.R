#Marginalized Schwarz-Arnason Jolly-Seber. The discrete latent states are marginalized
#out of the likelihood with a 3-state forward algorithm (1 = not yet entered,
#2 = alive, 3 = dead), but the M augmented individuals are retained so that
#we can use an individual random effect with no numerical integration.
#N, B, and N.super are recovered each iteration by forward-filtering backward-sampling
#the trajectories from their full conditionals.


#marginal likelihood of one capture history
#forward variables are rescaled every occasion and the log of the scale factor
#accumulated to avoid underflow
dJS <- nimbleFunction(
  run = function(x = double(1), eta = double(1), p = double(1), K = double(1),
                 phi = double(0), psi = double(0), log = integer(0, default = 0)){
    returnType(double(0))
    n.primary <- length(x)
    log.scale <- 0
    obs.alive <- dbinom(x[1], size = K[1], prob = p[1])
    obs.zero <- 0
    if(x[1] == 0){
      obs.zero <- 1
    }
    a1 <- (1-eta[1])*obs.zero
    a2 <- eta[1]*obs.alive
    a3 <- 0
    scale <- a1 + a2 + a3
    if(scale > 0){
      a1 <- a1/scale
      a2 <- a2/scale
      a3 <- a3/scale
      log.scale <- log.scale + log(scale)
    }
    if(n.primary > 1){
      for(g in 2:n.primary){
        obs.alive <- dbinom(x[g], size = K[g], prob = p[g])
        obs.zero <- 0
        if(x[g] == 0){
          obs.zero <- 1
        }
        n1 <- a1*(1-eta[g])*obs.zero
        n2 <- (a1*eta[g] + a2*phi)*obs.alive
        n3 <- (a2*(1-phi) + a3)*obs.zero
        scale <- n1 + n2 + n3
        if(scale > 0){
          a1 <- n1/scale
          a2 <- n2/scale
          a3 <- n3/scale
          log.scale <- log.scale + log(scale)
        }else{
          a1 <- 0
          a2 <- 0
          a3 <- 0
          log.scale <- -Inf
        }
      }
    }
    #log likelihood of the history given membership
    log.lik.member <- -Inf
    if(log.scale > -Inf){
      log.lik.member <- log.scale + log(a1 + a2 + a3)
    }
    all.zero <- 1
    for(g in 1:n.primary){
      if(x[g] > 0){
        all.zero <- 0
      }
    }
    #z.super marginalized: member w.p. psi, non-member contributes only to all-zero
    if(all.zero == 0){
      out <- log(psi) + log.lik.member
    }else{
      lp.mem <- log(psi) + log.lik.member
      lp.non <- log(1-psi)
      mx <- lp.mem
      if(lp.non > mx){
        mx <- lp.non
      }
      if(mx == -Inf){
        out <- -Inf
      }else{
        out <- mx + log(exp(lp.mem - mx) + exp(lp.non - mx))
      }
    }
    if(log){
      return(out)
    }
    return(exp(out))
  }
)

#stub, required for registration. never used with y supplied as data
rJS <- nimbleFunction(
  run = function(n = integer(0), eta = double(1), p = double(1), K = double(1),
                 phi = double(0), psi = double(0)){
    returnType(double(1))
    return(numeric(length(eta), value = 0))
  }
)


#Sampler for N, B, N.super
#forward variables rescaled each occasion, log scale factors accumulated to avoid underflow
NBSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    n.primary <- control$n.primary
    y <- control$y
    K <- control$K
    calcNodes <- control$calcNodes
  },
  run = function(){
    N.sim <- rep(0,n.primary)
    B.sim <- rep(0,n.primary)
    N.super.sim <- 0
    alpha <- matrix(0,nrow=n.primary,ncol=3)
    probs <- rep(0,3)
    for(i in 1:M){
      state <- rep(0,n.primary)
      phi <- model$phi[i]
      log.scale <- 0
      #forward filter
      obs.alive <- dbinom(y[i,1],size=K[1],prob=model$p[1])
      obs.zero <- 0
      if(y[i,1]==0){
        obs.zero <- 1
      }
      alpha[1,1] <- (1-model$eta[1])*obs.zero
      alpha[1,2] <- model$eta[1]*obs.alive
      alpha[1,3] <- 0
      scale <- alpha[1,1] + alpha[1,2] + alpha[1,3]
      if(scale > 0){
        alpha[1,1] <- alpha[1,1]/scale
        alpha[1,2] <- alpha[1,2]/scale
        alpha[1,3] <- alpha[1,3]/scale
        log.scale <- log.scale + log(scale)
      }else{
        log.scale <- -Inf
      }
      if(n.primary > 1){
        for(g in 2:n.primary){
          obs.alive <- dbinom(y[i,g],size=K[g],prob=model$p[g])
          obs.zero <- 0
          if(y[i,g] == 0){
            obs.zero <- 1
          }
          alpha[g,1] <- alpha[g-1,1]*(1-model$eta[g])*obs.zero
          alpha[g,2] <- (alpha[g-1,1]*model$eta[g] + alpha[g-1,2]*phi)*obs.alive
          alpha[g,3] <- (alpha[g-1,2]*(1-phi) + alpha[g-1,3])*obs.zero
          scale <- alpha[g,1] + alpha[g,2] + alpha[g,3]
          if(scale > 0){
            alpha[g,1] <- alpha[g,1]/scale
            alpha[g,2] <- alpha[g,2]/scale
            alpha[g,3] <- alpha[g,3]/scale
            log.scale <- log.scale + log(scale)
          }else{
            log.scale <- -Inf
          }
        }
      }
      #log likelihood of this history given membership.
      log.lik.member <- -Inf
      if(log.scale > -Inf){
        log.lik.member <- log.scale
      }
      all.zero <- 1
      for(g in 1:n.primary){
        if(y[i,g] > 0){
          all.zero <- 0
        }
      }
      #z.super from its full conditional, on the log scale
      lp.mem <- log(model$psi[1]) + log.lik.member
      p.member <- 0
      if(all.zero == 0){
        if(lp.mem > -Inf){
          p.member <- 1
        }
      }else{
        lp.non <- log(1-model$psi[1])
        mx <- lp.mem
        if(lp.non > mx){
          mx <- lp.non
        }
        if(mx > -Inf){
          p.member <- exp(lp.mem - mx)/(exp(lp.mem - mx) + exp(lp.non - mx))
        }
      }
      if(rbinom(1,1,p.member) == 1){
        N.super.sim <- N.super.sim + 1
        #backward sample the trajectory
        for(s in 1:3){
          probs[s] <- alpha[n.primary,s]
        }
        probs <- probs/sum(probs)
        state[n.primary] <- rcat(1,probs)
        if(n.primary > 1){
          #ascending loop index. nimble does not handle descending for loops
          for(gg in 1:(n.primary-1)){
            g <- n.primary - gg
            k <- state[g+1]
            #P(s_g = j | s_{g+1} = k) prop alpha[g,j]*trans[j,k]
            if(k == 1){
              probs[1] <- alpha[g,1]*(1-model$eta[g+1])
              probs[2] <- 0
              probs[3] <- 0
            }
            if(k == 2){
              probs[1] <- alpha[g,1]*model$eta[g+1]
              probs[2] <- alpha[g,2]*phi
              probs[3] <- 0
            }
            if(k == 3){
              probs[1] <- 0
              probs[2] <- alpha[g,2]*(1-phi)
              probs[3] <- alpha[g,3]
            }
            probs <- probs/sum(probs)
            state[g] <- rcat(1,probs)
          }
        }
        #tally. entry is the first occasion in state 2
        entered <- 0
        for(g in 1:n.primary){
          if(state[g] == 2){
            N.sim[g] <- N.sim[g] + 1
            if(entered == 0){
              B.sim[g] <- B.sim[g] + 1
              entered <- 1
            }
          }
        }
      }
    }
    model$N[1:n.primary] <<- N.sim
    model$B[1:n.primary] <<- B.sim
    model$N.super[1] <<- N.super.sim
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)