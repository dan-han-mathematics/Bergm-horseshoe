wd <- c("C:/Users/pamli/OneDrive - University of Louisville/Pamela's topic",
        "/Exponential Random Graph Models","/R codes (P)")


setwd(paste0(wd[1],wd[2],wd[3]))
source("./Horseshoe functions/HorseBERGM_v3.R")
source("./Horseshoe functions/plot.horseBERGM.R")
source("./Horseshoe functions/summary.horseBERGM.R")
source("./Horseshoe functions/bgof.R")
source("./2024.11.24_Pamela/KL_test/KL_functions.R")

library(statnet)
library(mcmc)
library(coda)
library(mvtnorm)

setwd("./2025.02.02_Pamela")

#--Network formula -------------------------------------------------------------

data(faux.magnolia.high)  ;  net <- faux.magnolia.high

form1 <- net ~ edges + nodematch(~Grade) + gwesp(0.5,fixed = T) +  kstar(2) + 
  gwdegree(0.5,fixed = T) + nodematch(~Sex) + nodematch(~Race) + triangle

(sy <- summary(form1))    ;   dims <- length(sy)
true_coef <- c(-7,1.5,rep(0,dims-2))
true_cov <- diag(1e-3,dims)

sim <- simulate(form1, nsim = 1 , basis = net, coef = true_coef, seed = 34567)

form2 <- sim ~ edges + nodematch(~Grade) + gwesp(0.5,fixed = T) +  kstar(2) + 
  gwdegree(0.5,fixed = T) + nodematch(~Sex) + nodematch(~Race) + triangle

(sy <- summary(form2))    ;   dims <- length(sy)

#-------------------------------------------------------------------------------

(mple <- coef(ergm(form2, estimate = "MPLE")))

main <- 1500  ;  burn <- 5000  ;  aux <- 1000  ;  chains <- 6
p.mean <- rep(0,dims)     ;  pads <- 0.35  ;  prop <- 0.00075

n <- 100

kld_horse <- rep(0,n)  ;  kld_bergm <- rep(0,n)
kld_bh <- rep(0,n)     ;  kld_hb <- rep(0,n)

horse_mean_est <- matrix(0, nrow = n, ncol = dims)
bergm_mean_est <- matrix(0, nrow = n, ncol = dims)

set.seed(12317)

for (i in 1:n) {
  cat("iteration",i,"\n")
  if (i == 1) print(Sys.time())
  
  horse1 <- horseBERGM(formula = form2,
                       prior.sigma = 10, 
                       lambda2 = 25, tau2 = 10,
                       zeta2 = c(1,1,rep(5e-5,dims-2)),
                       prior.mean = p.mean,
                       main.iters = main, 
                       burn.in = burn,
                       aux.iters = aux, 
                       gamma = pads,
                       V.proposal = prop,
                       nchains = chains)
  
  bergm1 <- Bergm::bergm(formula = form2,
                         prior.sigma = diag(10,dims),
                         prior.mean = p.mean,
                         main.iters = main, 
                         burn.in = burn,
                         aux.iters = aux, 
                         gamma = pads,
                         V.proposal = prop,
                         nchains = chains)
  
  horse_mean_est[i,] <- colMeans(horse1$Theta)
  bergm_mean_est[i,] <- colMeans(bergm1$Theta)
  
  kld_bergm[i] <- gaussian_kl(true_coef, true_cov, 
                              colMeans(bergm1$Theta), cov(bergm1$Theta))
  kld_horse[i] <- gaussian_kl(true_coef, true_cov, 
                              colMeans(horse1$Theta), cov(horse1$Theta))
  
  kld_bh[i] <- gaussian_kl(colMeans(bergm1$Theta), cov(bergm1$Theta),
                           colMeans(horse1$Theta), cov(horse1$Theta))
  
  kld_hb[i] <- gaussian_kl(colMeans(horse1$Theta), cov(horse1$Theta),
                           colMeans(bergm1$Theta), cov(bergm1$Theta))
}
gc()
save.image()


sink("Results_v1.txt")
options(digits = 4, scipen = 0)

cat("\n# KLD between true coef and both BHSERGM and BERGM, then between BERGM and BHSERGM #\n")
cat("# Used the faux.magnolia.high network with 8 network stats #\n")
cat("\n> formula <- sim_net ~ edges + nodematch(~Grade) + gwesp(0.5,fixed = T) +  kstar(2) +\n", 
    " ", " gwdegree(0.5,fixed = T) + nodematch(~Sex) + nodematch(~Race) + triangle\n")

cat("\n# Ran 100 times with: main.iters = 1500 ; burn.in = 5000 ; aux = 1000 ; chains = 6\n"," ",
    "   ", "   ","           p.mean = rep(0,dims) ; gamma = 0.35 ; V.prop = 0.0005\n")

cat("\n# Overall results of Horseshoe model:\n")
(horse_resuls <- data.frame(Post.Mean = colMeans(horse_mean_est),
                            Post.SD = apply(horse_mean_est, 2, sd),
                            Q_0.025 = apply(horse_mean_est,2,quantile,probs=0.025),
                            Q_0.25 = apply(horse_mean_est,2,quantile,probs=0.25),
                            Q_0.50 = apply(horse_mean_est,2,quantile,probs=0.50),
                            Q_0.75 = apply(horse_mean_est,2,quantile,probs=0.75),
                            Q_0.975 = apply(horse_mean_est,2,quantile,probs=0.975),
                            row.names = paste0("theta",1:8," (",horse1$specs,")")))

options(digits = 5, scipen = 0)

cat("\n# Overall results of Bergm model:\n")
(bergm_results <- data.frame(Post.Mean = colMeans(bergm_mean_est),
                             Post.SD = apply(bergm_mean_est, 2, sd),
                             Q_0.025 = apply(bergm_mean_est,2,quantile,probs=0.025),
                             Q_0.25 = apply(bergm_mean_est,2,quantile,probs=0.25),
                             Q_0.50 = apply(bergm_mean_est,2,quantile,probs=0.50),
                             Q_0.75 = apply(bergm_mean_est,2,quantile,probs=0.75),
                             Q_0.975 = apply(bergm_mean_est,2,quantile,probs=0.975),
                             row.names = paste0("theta",1:8," (",horse1$specs,")")))

cat("\n# KLD results:\n")
(kld_comp <- data.frame(mean = c(mean(kld_horse),mean(kld_bergm),
                                 mean(kld_bh),mean(kld_hb)),
                        sd = c(sd(kld_horse),sd(kld_bergm),
                               sd(kld_bh),sd(kld_hb)),
                        Q_0.025 = c(quantile(kld_horse, probs = 0.025),
                                    quantile(kld_bergm, probs = 0.025),
                                    quantile(kld_bh, probs = 0.025),
                                    quantile(kld_hb, probs = 0.025)),
                        Q_0.25 = c(quantile(kld_horse, probs = 0.25),
                                   quantile(kld_bergm, probs = 0.25),
                                   quantile(kld_bh, probs = 0.25),
                                   quantile(kld_hb, probs = 0.25)),
                        Q_0.50 = c(quantile(kld_horse, probs = 0.5),
                                   quantile(kld_bergm, probs = 0.5),
                                   quantile(kld_bh, probs = 0.5),
                                   quantile(kld_hb, probs = 0.5)),
                        Q_0.75 = c(quantile(kld_horse, probs = 0.75),
                                   quantile(kld_bergm, probs = 0.75),
                                   quantile(kld_bh, probs = 0.75),
                                   quantile(kld_hb, probs = 0.75)),
                        Q_0.975 = c(quantile(kld_horse, probs = 0.975),
                                    quantile(kld_bergm, probs = 0.975),
                                    quantile(kld_bh, probs = 0.975),
                                    quantile(kld_hb, probs = 0.975)),
                        row.names = c("true.vs.horse","true.vs.bergm",
                                      "bergm.vs.horse","horse.vs.bergm"))
)

sink()

