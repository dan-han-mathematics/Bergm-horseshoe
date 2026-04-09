
horseBERGM <- function(formula, 
                       prior.sigma = NULL, 
                       prior.mean = NULL,
                       lambda2 = NULL,
                       tau2 = 10,
                       zeta2 = 1,
                       aux.iters = 100,
                       main.iters = 1000,
                       burn.in = 200, 
                       V.proposal = 0.0025,
                       gamma = 0.8, 
                       nchains = NULL,
                       startVals = NULL,
                       offset.coef = NULL,
                       seed = NULL) {
  
  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy    <- summary(formula)
  dim   <- length(sy)
  
  tot.iters <- burn.in + main.iters
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initial settings
  control <- control.ergm(MCMC.burnin = aux.iters, 
                          MCMC.interval = 1,
                          MCMC.samplesize = 1)
    
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
    prior.sigma <- diag(lambda2*tau2,dim)
  } else { prior.sigma <- diag(prior.sigma,dim) }
  
  if (length(zeta2) == 1) {  zeta2 <- rep(zeta2, dim)  }
  lambda2 <- rep(lambda2, dim)
  
  thetaSamples <- array(NA, c(tot.iters,dim, nchains))
  sigma2Samples <- array(NA, c(tot.iters,1, nchains))
  lambda2Samples <- array(NA, c(tot.iters, dim, nchains))
  tau2Samples <- array(NA, c(tot.iters, 1, nchains))
  # nuSamples <- array(NA, c(tot.iters, dim, nchains))
  # xiSamples <- array(NA, c(tot.iters,1, nchains))
  
  proposal <- ergm_proposal(object = ~. , constraints = ~.,
                            arguments = control$MCMC.prop.args, nw = y)
  
  S.prop <- diag(V.proposal, dim, dim)
  Theta  <- array(NA, c(main.iters, dim, nchains))
  
  if (is.null(startVals)) {
    suppressMessages(mple <- coef(ergm(formula, estimate = "MPLE",
                                       verbose = FALSE,
                                       offset.coef = offset.coef)))
    
    if (any(!is.finite(mple))) mple[which(!is.finite(mple))] <- 0
    
    theta <- matrix(mple + runif(dim*nchains, min = -0.1, max = 0.1), 
                    dim, nchains)
  } else {
    theta <- matrix(startVals + runif(dim*nchains, min = -0.1, max = 0.1),
                    dim, nchains)   
    }  
  
  theta[model$etamap$offsettheta, ] <- offset.coef
  theta1 <- rep(NA, dim)
  
  if (is.null(prior.mean)) { 
    if (!exists("mple",inherits = F)) {
      prior.mean <- rep(0,dim)
    } else   prior.mean <- mple  
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
      
      beta <- (theta[ ,h] - theta1) %*% delta + pr[1] - pr[2]
      
      if ((beta >= log(runif(1))) & (!is.na(beta))) theta[ ,h] <- theta1
      
      thetaSamples[j,,h] <- theta[ ,h]
      
      # Set-up covariance matrix
      invD <- diag(1/lambda2,dim)
      invB <- invD/tau2
      
      # sample xi
      xi_b <- 1 + 1/tau2
      xi <- 1/rgamma(1, shape = 1, rate = xi_b)
      # xiSamples[j,,h] <- xi
            
      # sample nu_j
      nu_b <- 1/zeta2 + 1/lambda2
      nu <- 1/rgamma(dim, shape = 1, rate = nu_b)
      # nuSamples[j,,h] <- nu
      
      # sample sigma^2 
      shape <- dim/2
      scale <- t(theta[, h])%*%invB%*%theta[, h]/2
      sigma2 <- 1/rgamma(1, shape, rate = scale)
      sigma2Samples[j,,h] <- sigma2 
      prior.sigma <- sigma2*diag(lambda2*tau2,dim)
      
      # sample tau^2
      tau_a <- (dim+1)/2
      tau_b <- 1/xi + t(theta[, h])%*%invD%*%theta[, h]/(2*sigma2)
      tau2 <- 1/rgamma(1, tau_a, rate = tau_b)
      tau2Samples[j,,h] <- tau2
      
     # sample of lambda_j^2
      lamb_b <- 1/nu + theta[ ,h]^2/(2*sigma2*tau2)
      lambda2 <- 1/rgamma(dim, shape = 1, rate = lamb_b)
      lambda2Samples[j,,h] <- lambda2
    }
    
    if (j > burn.in) Theta[j - burn.in, , ] <- theta
  }
  
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  
  FF <- mcmc(na.omit(apply(Theta, 2, cbind)))
  lamb <- mcmc(na.omit(apply(lambda2Samples[-(1:burn.in),,], 2, cbind)))
  TAU <- mcmc(sapply(tau2Samples[-(1:burn.in),,], rbind))
  Sigma <- mcmc(sapply(sigma2Samples[-(1:burn.in),,], rbind))
  
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
  specs <- param_names(model)
  
  ess <- round(effectiveSize(FF), 0)
  names(ess) <- param_names(model)
  
  out <- list(Runtime = runtime, formula = formula, specs = specs, 
              dim = dim, Theta = FF, AR = AR, ess = ess, 
              sigma2 = Sigma, lambda2 = lamb,
              tau2 = TAU, seed = seed)
  
  # class(out) <- c("horseBERGM","bergm")
  
  return(out)
}
