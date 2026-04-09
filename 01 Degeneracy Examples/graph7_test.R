wd <- c("C:/Users/pamli/OneDrive - University of Louisville/Pamela's topic",
        "/Exponential Random Graph Models","/R codes (P)","/2024.")

setwd(paste0(wd[1],wd[2],wd[3]))
source("./Horseshoe functions/HorseBERGM_v3.R")
source("./Horseshoe functions/HorseBERGM_w_sigma2.R")
source("./Horseshoe_v2 functions/gam_HorseBERGM_v1.R")
source("./Horseshoe functions/plot.horseBERGM.R")
source("./Horseshoe functions/summary.horseBERGM.R")
source("./Horseshoe functions/bgof.R")

setwd("./2024.05.10_Pamela/graph7_degen")

library(mvtnorm)
library(statnet)
library(mcmc)
library(coda)

#-------------------------------------------------------------------------------
#  Create network on 7 nodes

net <- network(7, directed = FALSE)

# calculate all possible vectors of statistics
stats <- ergm.allstats(net ~ edges + triangles, zeroobs = F)

# color function
summary(stats$weights)

cr <- colorRamp(c("#E8E8E8","#3A3A3A"))

stats$cols <- ifelse(stats$weight < 200, "white", 
                     rgb(cr(stats$weights/max(stats$weights)), max=255))

plot(stats$statmat, pch=21, cex=1.5, bg = stats$cols, 
     col = ifelse(stats$weights < 200, "#DFDFDF",stats$cols), 
     xlab ="Number of edges", ylab="Number of triangles",
     main = "Graphs on N = 7 nodes:\n Classified by edges and triangles")
library(rgeos)
plot(gConvexHull(sp::SpatialPoints(stats$statmat)), add = TRUE)

#-------------------------------------------------------------------------------
# Create desire network with 10 edges and 10 triangles
set.seed(1234)
net <- san(net ~ edges + triangles, target.stats = c(10,10))
summary(f1 <- formula(net ~ edges + triangles))

plot(net)

#-- Fitting models -------------------------------------------------------------

## ERGM MLE estimation
mle <- ergmito(f1, ntries = 10)

ergm.exact(coef(mle), f1, statmat = stats$statmat-10, weights = stats$weights)

## Horseshoe w/ Cauchy prior
horse1 <- horseBERGM(formula = f1,
                     prior.sigma = 30, 
                     lambda2 = 25, tau2 = 100,
                     startVals = c(0,0),
                     prior.mean = c(0,0),
                     main.iters = 20000, burn.in = 20000,
                     aux.iters = 10000, gamma = 1.03,
                     V.proposal = 0.01, nchains = 6,
                     seed = 331)
# summary.horseBERGM(horse1)
pdf("horse1_mcmc.pdf", height = 6, width = 8)
plot.horseBERGM(horse1, lag.max = 200)
dev.off()

## BERGM
set.seed(335)
bergm1 <- Bergm::bergm(formula = f1,
                       prior.sigma = diag(30,2),  
                       startVals = c(0,0), prior.mean = c(0,0),
                       main.iters = 20000, burn.in = 20000,
                       aux.iters = 10000, gamma = 1.03,
                       V.proposal = 0.01, nchains = 6)
# summary.horseBERGM(bergm1)
pdf("bergm1_mcmc.pdf", height = 6, width = 8)
plot.horseBERGM(bergm1,lag.max = 200)
dev.off()


## Generalized horseshoe with cauchy prior
horse2 <- gen_horseBERGM(formula = f1,
                         prior.sigma = 30, 
                         lambda2 = 25, tau2 = 1,
                         startVals = c(0,0),
                         prior.mean = c(0,0),
                         main.iters = 20000, burn.in = 20000,
                         aux.iters = 10000, gamma = 1.03,
                         V.proposal = 0.01, nchains = 6,
                         seed = 337)
# summary.horseBERGM(horse2)
# plot.horseBERGM(horse2,lag.max = 150)


# Horseshoe w/ Gamma prior
horse3 <- gp_horseBERGM(formula = f1,
                        prior.sigma = 30, 
                        invLambda2 = 25, tau2 = 0.001,
                        nu_j = 1.5, prior.mean = 0,
                        startVals = c(0,0),
                        main.iters = 20000, burn.in = 20000,
                        aux.iters = 10000, gamma = 1.03,
                        V.proposal = 0.01, nchains = 6,
                        seed = 339)
# summary.horseBERGM(horse3)
# plot.horseBERGM(horse3,lag.max = 150)


save.image()

#--- Simulation using model coefs ----------------------------------------------

library(dplyr)

cr <- colorRamp(c("#E8E8E8","#3A3A3A"))

# Storing seeds used during model estimation
seeds <- c('ergm_mle' = 123, 'cauchy_horse' = 331,
           'bergm' = 335, 'sigma_horse' = 337, 'gamma_horse' = 339)

# Storing estimated theta coefficients to used in simulation
models_est <- list('ergm_mle' = coef(mle),
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
    group_by(edges,triangle) %>%
    summarise(num = n(), .groups = "drop") %>%
    mutate(cols = rgb(cr(num/max(num)), max=255))
}


# Add the model classes of all graphs on 7 nodes -- 110 total classes
mod_class[['graph7']] <- stats$statmat %>%
  bind_cols(num = stats$weights,
            cols = stats$cols) %>%
  arrange(edges,triangle)


#--- Plots ---------------------------------------------------------------------
# Plot for main horseshoe, bergm, ergm, and general graphs on 7 nodes

pdf("graph7_class_vs_ergm_mle.pdf",height = 5.5,width = 10.5)
par(mfrow = c(1,2), mar = c(3.6,3.8,2.5,0.8), cex = 1, oma = c(0,0,0,0))

with(mod_class[['graph7']], {
  plot(edges, triangle, pch = 21, cex = 1.65, 
       bg = cols, 
       col = ifelse(num <= 200, "#DFDFDF",cols),  
       xlab = " ", ylab = " ", xlim = c(0,21), ylim = c(0,35),
       main = "Graphs on N = 7 nodes:\n Classified by edges and triangles")
  })
plot(gConvexHull(sp::SpatialPoints(stats$statmat)), add = TRUE, lwd = 0.85)
mtext("Number of triangles", side = 2, line = 2.45, cex = 1.1)
mtext("Number of edges", side = 1, line = 2.45, cex = 1.1)

with(mod_class[['ergm_mle']], {
  plot(edges, triangle, pch = 21, cex = 1.65, 
       bg = ifelse(num <= 1, "white",cols), 
       col = ifelse(num <= 1, "#DFDFDF",cols),   
       xlab = " ", ylab = " ", xlim = c(0,21), ylim = c(0,35),
       main = "Graph (N=7) generated from ergm mle") 
  })
mtext("Number of triangles", side = 2, line = 2.45, cex = 1.1)
mtext("Number of edges", side = 1, line = 2.45, cex = 1.1)
points(10,10,pch=19,col="red",cex=0.4)
dev.off()

pdf("graph7_horse_vs_bergm_class.pdf", height = 5.5,width = 10.5)
par(mfrow = c(1,2), mar = c(3.6,3.8,2.5,0.8), cex = 1,oma = c(0,0,0,0))

with(mod_class[['cauchy_horse']], {
  plot(edges, triangle, pch = 21, cex = 1.65, 
       bg = ifelse(num <= 1, "white",cols), 
       col = ifelse(num <= 1, "#DFDFDF",cols),   
       xlab = " ", ylab = " ", xlim = c(0,21), ylim = c(0,35),
       main = "Graphs (N=7) generated from horseBergm")
  })
mtext("Number of triangles", side = 2, line = 2.45, cex = 1.1)
mtext("Number of edges", side = 1, line = 2.45, cex = 1.1)
points(10,10,pch=19,col="red",cex=0.5,lwd=0.5)

with(mod_class[['bergm']], {
  plot(edges, triangle, pch = 21, cex = 1.65, 
       bg = ifelse(num <= 1, "white",cols), 
       col = ifelse(num <= 1, "#DFDFDF",cols),  
       xlab = " ", ylab = " ", xlim = c(0,21), ylim = c(0,35),
       main = "Graph (N=7) generated from bergm coef")
  })
mtext("Number of triangles", side = 2, line = 2.45, cex = 1.1)
mtext("Number of edges", side = 1, line = 2.45, cex = 1.1)
points(10,10,pch=19,col="red",cex=0.5,lwd=0.5)
dev.off()

pdf("histogram_3model_comp.pdf", width = 10, height = 7)
par(mfrow = c(3,2), mar = c(2,4,2.5,0.01),oma = c(1.5,0,0.1,0.6))
for (i in 1:3) {
  hist(model_sims[[i]][[1]][,'edges'], xlab="", ylab="", xlim = c(0,21), 
       main = ifelse(i == 1, paste("Density of edges"),""), breaks = 20)
  mtext("edges", side = 1, line = 2.1, cex=0.8, font = 1)
  mtext(names(model_sims)[i], side = 2, line = 2.6, cex=0.8, font = 2)
  abline(v = 10, col = "red", lty = 1, lwd = 1.25)
  abline(v = mean(model_sims[[i]][[1]][,'edges']), 
         col = "blue", lty = 1, lwd = 1.25)
  
  hist(model_sims[[i]][[1]][,2], xlab="", ylab="", xlim = c(0,35),
       main = ifelse(i == 1, paste("Density of triangle"),""))
  mtext("triangle", side = 1, line = 2, cex=0.8, font = 1)
  abline(v = 10, col = "red", lty = 1, lwd = 1.25)
  abline(v = mean(model_sims[[i]][[1]][,2]), col = "blue", lty = 1, lwd = 1.25)
}
dev.off()

#-------------------------------------------------------------------------------
## OTHER HORSESHOE PLOTS

with(mod_class[['sigma_horse']], {
  plot(edges, triangle, pch=21, cex = 2, xlim = c(0,21), ylim = c(0,35),
       bg = cols, col = ifelse(num <= 10, "#EAEAEA",cols), 
       xlab ="Number of edges", ylab="Number of triangles",
       main = paste("Graphs (N=7) generated from horseBergm",
                    "coef:\n Classified by edges and triangles"))   })
points(10,10,pch=19,col="red",cex=0.6)

with(mod_class[['gamma_horse']], {
  plot(edges, triangle, pch=21, cex = 2, xlim = c(0,21), ylim = c(0,35),
       bg = cols, col = ifelse(num <= 20, "#EAEAEA",cols), 
       xlab ="Number of edges", ylab="Number of triangles",
       main = "Graphs (N=7) generated from gen_horseBergm coef") })
points(10,10,pch=19,col="red",cex=0.6)


#-------------------------------------------------------------------------------
# Saving results

sink("g7_test1_v1_results.txt")

cat("# Degeneracy model - Graph on 7 nodes with 10 edges and 10 triangles #\n")
cat("\n> formula <- net ~ edges + triangles()\n")

cat("\n# Results of Horseshoe with half-caucy prior for lambda")
cat("\n> summary.horseBERGM(horse1)\n")
summary.horseBERGM(horse1)

cat("\n# Results of Bergm")
cat("\n> summary.horseBERGM(bergm1)\n")
summary.horseBERGM(bergm1)

cat("\n# Effective Sample Size:\n")

rbind("horseshoe"=horse1$ess,"bergm"=bergm1$ess)

cat("\n\n# Simulation results -- 20000 sims using estimates and MLE #\n")
cat("\n# Summary of network statistics\n")
cat("\n> horseshoe_sim_stats\n")
summary(model_sims[['cauchy_horse']])

cat("\n\n> bergm_sim_stats\n")
summary(model_sims[['bergm']])

cat("\n\n> MLE_sim_stats\n")
summary(model_sims[['ergm_mle']])

cat("\n\n# Number of classes, with max class by model #\n\n")

cat("> class_horse[which.max(class_horse$count), ]   (105 total classes)\n")
mod_class[['cauchy_horse']][which.max(mod_class[['cauchy_horse']]$num), 1:3] 

cat("\n\n> class_bergm[which.max(class_bergm$count), ]   (99 total classes)\n")
mod_class[['bergm']][which.max(mod_class[['bergm']]$num), 1:3] 

cat("\n\n> class_ergm[which.max(class_bergm$count), ]   (98 total classes)\n")
mod_class[['ergm_mle']][which.max(mod_class[['ergm_mle']]$num), 1:3]


cat("\n\n# main.iters = 20000, burn.in = 20000, aux.iters = 10000",
    "\n# gamma = 1.03,    V.proposal = 0.01,    nchains = 6\n\n")

sink()

save.image("g7_v1.Rdata")
