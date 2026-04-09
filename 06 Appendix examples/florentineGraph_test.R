wd <- c("C:/Users/pamli/OneDrive - University of Louisville/Pamela's topic",
        "/Exponential Random Graph Models","/R codes (P)","/2024.")

setwd(paste0(wd[1],wd[2],wd[3]))
source("./Horseshoe functions/HorseBERGM_v3.R")
source("./Horseshoe functions/HorseBERGM_w_sigma2.R")
source("./Horseshoe_v2 functions/gam_HorseBERGM_v1.R")
source("./Horseshoe functions/plot.horseBERGM.R")
source("./Horseshoe functions/summary.horseBERGM.R")
source("./Horseshoe functions/bgof.R")

setwd("./2024.05.10_Pamela/florentine_degen")

library(mvtnorm)
library(statnet)
library(mcmc)
library(coda)

#-------------------------------------------------------------------------------
# Data

data(florentine)
net <- flobusiness     ;     rm(flomarriage)

library(GGally)

pdf("flobusiness.pdf")
set.seed(1233)
ggnet2(net, size = "degree",
       color = "steelblue",
       mode = "fruchtermanreingold", 
       edge.size = 0.8)
ggsave("flobusiness.pdf")   

(sy <- summary(f1 <- formula(net ~ edges + kstar(2))))

#-------------------------------------------------------------------------------
# Models

# ERGM MPLE estimate
control <- control.ergm(MPLE.covariance.method='Godambe',
                        seed = 1711)
mple <- ergm(f1, estimate = "MPLE", control = control)

## Horseshoe w/ Cauchy prior
horse1 <- horseBERGM(formula = f1, prior.sigma = 30, 
                     lambda2 = 25, tau2 = 10,
                     startVals = c(0,0), prior.mean = c(0,0),
                     main.iters = 20000, burn.in = 20000,
                     aux.iters = 5000, gamma = 1.03,
                     V.proposal = 0.0025, nchains = 6,
                     seed = 1719)
# summary.horseBERGM(horse1)
# plot.horseBERGM(horse1,lag.max = 150)


## Generalized horseshoe with cauchy prior
horse2 <- gen_horseBERGM(formula = f1,
                         prior.sigma = 30, lambda2 = 25, tau2 = 100,
                         startVals = c(0,0), prior.mean = c(0,0),
                         main.iters = 20000, burn.in = 20000,
                         aux.iters = 5000, gamma = 1.03,
                         V.proposal = 0.0025, nchains = 6,
                         seed = 1713)
# summary.horseBERGM(horse2)
# plot.horseBERGM(horse2,lag.max = 150)


# Horseshoe w/ Gamma prior
horse3 <- gp_horseBERGM(formula = f1,
                        prior.sigma = 30, 
                        invLambda2 = 25, tau2 = 0.001,
                        nu_j = 1.5, prior.mean = 0,
                        startVals = c(0,0), 
                        main.iters = 20000, burn.in = 20000,
                        aux.iters = 5000, gamma = 1.03,
                        V.proposal = 0.0025, nchains = 6,
                        seed = 1715)
# summary.horseBERGM(horse3)
# plot.horseBERGM(horse3,lag.max = 150)


## BERGM
set.seed(1717)
bergm1 <- Bergm::bergm(formula = f1,
                       prior.sigma = diag(30,2),  
                       startVals = c(0,0), prior.mean = c(0,0),
                       main.iters = 20000, burn.in = 20000,
                       aux.iters = 5000, gamma = 1.03,
                       V.proposal = 0.0025, nchains = 6)
# summary.horseBERGM(bergm1)
# plot.horseBERGM(bergm1,lag.max = 150)

save.image()

#-------------------------------------------------------------------------------
# Simulations

library(dplyr)

# Storing seeds used during model estimation
seeds <- c('ergm_mple' = 130,'cauchy_horse' = 1719,
           'bergm' = 1717, 'sigma_horse' = 1713, 'gamma_horse' = 1715)

# Storing estimated theta coefficients to used in simulation
models_est <- list('ergm_mple' = coef(mple),
                   'cauchy_horse' = colMeans(horse1$Theta),
                   'bergm' = colMeans(bergm1$Theta),
                   'sigma_horse' = colMeans(horse2$Theta),
                   'gamma_horse' = colMeans(horse3$Theta))

# Empty list to store simulations
model_sims <- list()

# Empty list to calculate the classes of simulated graphs
mod_class <- list()

for (i in names(models_est)) {
  model_sims[[i]] <- simulate(f1, coef = models_est[[i]], 
                              nsim = 20000, output = "stats",
                              seed = seeds[i],
                              simplify = F)
  mod_class[[i]] <- as.data.frame(model_sims[[i]][[1]]) %>%
    group_by(edges,kstar2) %>%
    summarise(num = n(), .groups = "drop")
}

#--- Plots ---------------------------------------------------------------------
# Plot for main horseshoe, bergm, ergm, and general graphs on 7 nodes

pdf("flo_horse1_mcmc.pdf",width = 8, height = 7)
plot.horseBERGM(horse1, lag = 150)
dev.off()

pdf("flo_bergm1_mcmc.pdf",width = 8, height = 7)
plot.horseBERGM(bergm1, lag = 150)
dev.off()

pdf("histogram_2model_comp.pdf", width = 10, height = 7)
par(mfrow = c(2,2), mar = c(2,4,2,0.1),oma = c(1.5,0,1,0.5))
for (i in 2:3) {
  hist(model_sims[[i]][[1]][,'edges'], xlab="", ylab="", xlim = c(0,35), 
       main = ifelse(i == 2, paste("Density of edges"),""), breaks = 30)
  mtext(names(model_sims)[i], side = 2, line = 2.6, cex=0.8, font = 2)
  mtext("edges", side = 1, line = 2, cex=0.8, font = 1)
  abline(v = 15, col = "red", lty = 1, lwd = 1.3)
  abline(v = mean(model_sims[[i]][[1]][,1]),col = "blue", lty = 1, lwd = 1.3)
  hist(model_sims[[i]][[1]][,2], xlab="", ylab="", xlim = c(0,150),
       main = ifelse(i == 2, paste("Density of two-star")),breaks = 40)
  mtext("2-star", side = 1, line = 2, cex=0.8, font = 1)
  abline(v = 36, col = "red", lty = 1, lwd = 1.3)
  abline(v = mean(model_sims[[i]][[1]][,2]),col = "blue", lty = 1, lwd = 1.3)
}
dev.off()


#-------------------------------------------------------------------------------
# Saving results

sink("test1_results.txt")

cat("# Degeneracy model - Florentine graph on 16 nodes with 15 edges and 36 two-stars #\n")
cat("\n> formula <- net ~ edges + kstar(2)\n")

cat("\n# Results of Horseshoe with half-caucy prior for lambda")
cat("\n> summary.horseBERGM(horse1)\n")
summary.horseBERGM(horse1)

cat("\n# Results of Bergm")
cat("\n> summary.horseBERGM(bergm1)\n")
summary.horseBERGM(bergm1)

cat("\n# Estimated MPLE:", paste0("(",coef(mple)[1],", ",coef(mple)[2],")"))

cat("\n\n# Effective Sample Size:\n")

rbind("horseshoe"=horse1$ess,"bergm"=bergm1$ess)

cat("\n\n# Simulation results -- 20000 sims using estimates and MPLE #\n")
cat("\n# Summary of network statistics\n")
cat("\n> horseshoe_sim_stats\n")
summary(model_sims[['cauchy_horse']])

cat("\n\n> bergm_sim_stats\n")
summary(model_sims[['bergm']])

cat("\n\n> ergmMLE_sim_stats\n")
summary(model_sims[['ergm_mple']])

cat("\n\n# Number of classes, with max class by model #\n\n")

cat("> class_horse[which.max(class_horse$count), ]   (106 total classes)\n")
mod_class[['cauchy_horse']][which.max(mod_class[['cauchy_horse']]$num), 1:3] 

cat("\n\n> class_bergm[which.max(class_bergm$count), ]   (98 total classes)\n")
mod_class[['bergm']][which.max(mod_class[['bergm']]$num), 1:3] 

cat("\n\n> class_ergm[which.max(class_bergm$count), ]   (97 total classes)\n")
mod_class[['ergm_mple']][which.max(mod_class[['ergm_mple']]$num), 1:3]


cat("\n\n# main.iters = 20000, burn.in = 20000, aux.iters = 5000",
    "\n# gamma = 1.03,    V.proposal = 0.0025,    nchains = 6\n\n")

sink()

save.image(".Rdata")
