
HE_horseBERGM <- function(formula, 
                       prior.sigma = NULL, 
                       prior.mean = NULL,
                       lambda2 = NULL,
                       tau2 = 100,
                       zeta2 = 1,
                       aux.iters = 100,
                       main.iters = 1000,
                       burn.in = 200, 
                       V.proposal = 0.0025,
                       gamma = 0.8, 
                       nchains = NULL,
                       startVals = NULL,
                       offset.coef = NULL,
                       thin = 1,
                       seed = NULL) {
  
  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy    <- summary(formula)
  dim   <- length(sy)
  
  tot.iters <- burn.in + main.iters
  
  # Initial settings
  
  if (is.null(seed)) {
    seed <- sample(1:99,1)
  }
  
  set.seed(seed)
  
  control <- control.ergm(MCMC.burnin = aux.iters, 
                          MCMC.interval = 1,
                          MCMC.samplesize = 1,
                          seed = seed)
  
  if (!is.null(control$init)) {
    if (length(control$init) != length(model$etamap$offsettheta)) {
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.", 
                 "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  } else { control$init <- rep(NA, length(model$etamap$offsettheta)) }
  
  if (!is.null(offset.coef)) {
    if (length(control$init[model$etamap$offsettheta]) != length(offset.coef)) {
      stop("Invalid offset parameter vector offset.coef: ",
           "wrong number of parameters: expected ",
           length(control$init[model$etamap$offsettheta]),
           " got ", length(offset.coef), ".") 
    }
    control$init[model$etamap$offsettheta] <- offset.coef   }
  
  if (is.null(lambda2)) lambda2 <- 25
  if (is.null(nchains)) nchains <- 2*dim
  if (is.null(prior.sigma)) {
    prior.sigma <- diag(100,dim)
  } else { prior.sigma <- diag(prior.sigma,dim) }
  
  lambda2 <- rep(lambda2, dim)
  zeta2 <- rep(zeta2, dim)
  
  thetaSamples <- array(NA, c(tot.iters,dim, nchains))
  sigma2Samples <- array(NA, c(tot.iters,1, nchains))
  lambda2Samples <- array(NA, c(tot.iters, dim, nchains))
  tau2Samples <- array(NA, c(tot.iters, 1, nchains))
  nuSamples <- array(NA, c(tot.iters, dim, nchains))
  xiSamples <- array(NA, c(tot.iters,1, nchains))
  
  proposal <- ergm_proposal(object = ~. , constraints = ~.,
                            arguments = control$MCMC.prop.args, nw = y)
  
  S.prop <- diag(V.proposal, dim, dim)
  Theta  <- array(NA, c(main.iters, dim, nchains))
  
  if (is.null(startVals)) {
    suppressMessages(mple <- coef(ergm(formula, estimate = "MPLE",
                                       verbose = FALSE,
                                       offset.coef = offset.coef)))
    
    theta <- matrix(mple + runif(dim*nchains, min = -0.1, max = 0.1), 
                    dim, nchains)
  } else {
    theta <- matrix(startVals + runif(dim*nchains, min = -0.1, max = 0.1),
                    dim, nchains)   
  }  
  
  theta[model$etamap$offsettheta, ] <- offset.coef
  na_theta <- which(!is.finite(theta))
  theta[na_theta] <- theta[na_theta - 1]
  theta1 <- rep(NA, dim)
  
  if (is.null(prior.mean)) { 
    if (!exists("mple")) {
      prior.mean <- rep(0,dim)
    } else {
      prior.mean <- mple
      prior.mean[which(!is.finite(prior.mean))] <- 0 
    }
  }
  
  y0 <- simulate(formula,
                 coef     = theta[,1],
                 nsim     = 1,
                 control  = control.simulate(MCMC.burnin = 1,    # !!!
                                             MCMC.interval = 1),
                 return.args = "ergm_state")$object
  
  message(" > MCMC start")
  
  clock.start <- Sys.time()
  
  for (j in 1:tot.iters) {
    
    for (h in 1:nchains) {
      
      theta1 <- theta[, h] + 
        gamma*apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff) + 
        rmvnorm(1, sigma = S.prop)[1,]
      
      delta <- ergm_MCMC_sample(y0, 
                                theta   = theta1,
                                control = control)$stats[[1]][1, ] - sy
      
      pr <- dmvnorm(rbind(theta1, theta[, h]), 
                    mean = prior.mean,
                    sigma = prior.sigma, 
                    log = TRUE)
      
      beta <- (theta[ ,h] - theta1)%*%delta + pr[1] - pr[2]
      
      if ((beta >= log(runif(1))) & (!is.na(beta))) theta[ ,h] <- theta1
      
      thetaSamples[j,,h] <- theta[ ,h]
      
      # Set-up covariance matrix
      invD <- diag(1/lambda2,dim)
      invB <- invD/tau2
      
      # sample xi
      xi_a <- 1
      xi_b <- 1 + 1/tau2
      xi <- 1/rgamma(1, xi_a, rate = xi_b)
      xiSamples[j,,h] <- xi
      
      # sample nu_j
      nu_a <- 1
      nu_b <- 1/zeta2 + 1/lambda2
      nu <- 1/rgamma(dim, nu_a, rate = nu_b)
      nuSamples[j,,h] <- nu
      
      # sample sigma^2 
      shape <- dim/2
      scale <- t(theta[, h])%*%invB%*%theta[, h]/2
      sigma2 <- 1/rgamma(1, shape, rate = scale)
      prior.sigma <- sigma2*diag(lambda2*tau2,dim)
      sigma2Samples[j,,h] <- sigma2
      
      # sample tau^2
      tau_a <- (dim+1)/2
      tau_b <- 1/xi + t(theta[, h])%*%invD%*%theta[, h]/(2*sigma2)
      tau2 <- 1/rgamma(1, tau_a, rate = tau_b)
      tau2Samples[j,,h] <- tau2
      
      # sample of lambda_j^2
      lamb_a <- 1
      lamb_b <- 1/nu + theta[ ,h]^2/(2*sigma2*tau2)
      lambda2 <- 1/rgamma(dim, lamb_a, rate = lamb_b)
      lambda2Samples[j,,h] <- lambda2
      
      # empirical sample of zeta_j^2
      if (j %% 5 == 0) {
        low <- j - 4
        high <- j
        zeta2 <- 2*colMeans(1/nuSamples[low:high,,h])
      }
    }
    
    if (j > burn.in) Theta[j - burn.in, , ] <- theta
    }

  
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  
  FF <- mcmc(na.omit(apply(Theta, 2, cbind)),thin = thin)
  lamb <- mcmc(na.omit(apply(lambda2Samples[-(1:burn.in),,], 2, cbind)),thin = thin)
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
  specs <- param_names(model)
  
  ess <- round(effectiveSize(FF), 0)
  names(ess) <- param_names(model)
  
  out <- list(Runtime = runtime, formula = formula, specs = specs,
              dim = dim, Theta = FF, AR = AR, ess = ess, lambda = lamb)
  
  # class(out) <- c("horseBERGM","bergm")
  
  return(out)
}
