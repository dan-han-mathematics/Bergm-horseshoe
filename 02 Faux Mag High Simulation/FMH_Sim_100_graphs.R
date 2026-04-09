wd <- c("C:/Users/pamli/OneDrive - University of Louisville/Pamela's topic",
        "/Exponential Random Graph Models","/R codes (P)","/2023.")


setwd(paste0(wd[1],wd[2],wd[3]))
source("./Horseshoe functions/HorseBERGM_v3.R")
source("./Horseshoe functions/plot.horseBERGM.R")
source("./Horseshoe functions/summary.horseBERGM.R")
source("./Horseshoe functions/bgof.R")

library(mvtnorm)
library(statnet)
library(mcmc)
library(coda)

library(matrixStats)

setwd(paste0(".","/2024.11.24_Pamela"))

#--- Network Plot--- ----------------------------------------------------------

data("faux.magnolia.high")
net <- faux.magnolia.high

#--- Simulation setup ----------------------------------------------------------

# Number of simulations
n <- 100

# Vector of acceptance rate, running time, and effective sample size
AccRateBERGM <- rep(0,n)
AccRateHorse <- rep(0,n)

Timebergm <- rep(0,n)
TimeHorse <- rep(0,n)

essbergm <- matrix(0, nrow=n, ncol=3)
essHorse <- matrix(0, nrow=n, ncol=3)


# Matrix to store estimated median and mean
bergm_estimate_median <- matrix(0, nrow=n, ncol=3)
horse_estimate_median <- matrix(0, nrow=n, ncol=3)

bergm_estimate_mean <- matrix(0, nrow=n, ncol=3)
horse_estimate_mean <- matrix(0, nrow=n, ncol=3)
kld <- rep(0,n)

form1 <- net ~ edges + nodematch("Grade") + gwesp(0.25,fixed=T)
set.seed(2171)

for (i in 1:n) {
  
  sim <- simulate(form1, coef = c(-7.8,2.4,1.5))
  form2 <- sim ~ edges + nodematch("Grade") + gwesp(0.25,fixed=T)
 
  print(i)
  
  horse1 <- horseBERGM(formula = form2, 
                       prior.sigma = 100, 
                       lambda2 = 25, tau2 = 100,
                       main.iters = 2000,
                       aux.iters = 75, burn.in = 150,
                       V.proposal = 0.05,
                       gamma = 0.67, nchains = 10)
  
  bergm1 <- Bergm::bergm(formula = form2, main.iters = 2000,
                         aux.iters = 75, burn.in = 150,
                         V.proposal = 0.05,
                         gamma = 0.67, nchains = 10)
  
  bergm_estimate_median[i,] <- colMedians(bergm1$Theta)
  horse_estimate_median[i,] <- colMedians(horse1$Theta)
  
  bergm_estimate_mean[i,] <- colMeans(bergm1$Theta)
  horse_estimate_mean[i,] <- colMeans(horse1$Theta)

  AccRateBERGM[i] <- bergm1$AR
  AccRateHorse[i] <- horse1$AR
  
  Timebergm[i] <- bergm1$Time
  TimeHorse[i] <- horse1$Runtime
  
  essbergm[i,] <- as.vector(bergm1$ess)
  essHorse[i,] <- as.vector(horse1$ess)
  
  kld[i] <- gaussian_kl(colMeans(horse1$Theta), cov(horse1$Theta),
                        colMeans(bergm1$Theta), cov(bergm1$Theta))
}

save.image("v2.Rdata")

#--- Overall Mean & Median -----------------------------------------------------

sink("Sim_100_Results.txt")

cat("# Simulation of",n,"networks results #\n")
cat("\n> formula <-",paste(form2)[c(2,1,3)],"\n") 
cat("\n# Mean and Median estimate [ coef = c(-6.8,2.4,1.5) ]:")
cat("\n> colMeans(bergm_estimate_median)\n", colMeans(bergm_estimate_median))
cat("\n> colMeans(bergm_estimate_mean)\n",   colMeans(bergm_estimate_mean))
cat("\n")
cat("\n> colMeans(horse_estimate_median)\n", colMeans(horse_estimate_median))
cat("\n> colMeans(horse_estimate_mean)\n",   colMeans(horse_estimate_mean))


#--- Mean and Median Quantiles -------------------------------------------------

quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)

bergm_median_quantile <- t(apply(bergm_estimate_median,2,quantile,quantiles))
horse_median_quantile <- t(apply(horse_estimate_median,2,quantile,quantiles))

bergm_mean_quantile <- t(apply(bergm_estimate_mean,2,quantile,quantiles))
horse_mean_quantile <- t(apply(horse_estimate_mean,2,quantile,quantiles))

cat("\n")
cat("\n\n# Quantiles based on Mean and Medians:")
cat("\n> bergm_median_quantile\n")  ; bergm_median_quantile
cat("\n> horse_median_quantile\n")  ; horse_median_quantile
cat("\n")
cat("\n> bergm_mean_quantile\n")    ; bergm_mean_quantile
cat("\n> horse_mean_quantile\n")    ; horse_mean_quantile


#--- Mean and Median AccRate ---------------------------------------------------

cat("\n\n# Mean and Median Acceptace Rate:")

cat("\n> median(AccRateBERGM): ",median(AccRateBERGM))
cat("\n> median(AccRateHorse): ",median(AccRateHorse))
cat("\n")
cat("\n> mean(AccRateBERGM): ", mean(AccRateBERGM))
cat("\n> mean(AccRateHorse): ", mean(AccRateHorse), "\n")

#--- Error and MSE -------------------------------------------------------------

true_value <- matrix(rep(c(-6.8,2.4,1.5),n),nrow=n,ncol=3, byrow=T)

error_bergm_mean <- bergm_estimate_mean - true_value
error_horse_mean <- horse_estimate_mean - true_value

MSE_bergm_mean <- mean((error_bergm_mean)^2)
MSE_horse_mean <- mean((error_horse_mean)^2)

error_bergm_median <- bergm_estimate_median - true_value
error_horse_median <- horse_estimate_median - true_value

MSE_bergm_median <- mean((error_bergm_median)^2)
MSE_horse_median <- mean((error_horse_median)^2)


cat("\n\n# MSE based of true vales of c(-6.8,2.4,1.5):")
cat("\n> MSE_bergm_mean: ", MSE_bergm_mean)
cat("\n> MSE_horse_mean: ", MSE_horse_mean)
cat("\n")
cat("\n> MSE_bergm_median: ", MSE_bergm_median)
cat("\n> MSE_horse_median: ", MSE_horse_median, "\n")


#--- Avg effective sample size -------------------------------------------------
Avg_ess_bergm <- colMeans(essbergm)
Avg_ess_horse <- colMeans(essHorse)

cat("\n\n# Average effective sample size:")
cat("\n> Avg_ess_bergm\n", Avg_ess_bergm,"\n")
cat("\n> Avg_ess_horse\n", Avg_ess_horse, "\n\n")

#--- Avg running time  ---------------------------------------------------------
Avg_time_bergm <- mean(Timebergm)
Avg_time_horse <- mean(TimeHorse)

cat("\n# Average running time:")
cat("\n> Avg_time_bergm: ", Avg_time_bergm)
cat("\n> Avg_time_horse: ", Avg_time_horse)

sink()

#--- Plots and BGOF ------------------------------------------------------------
pdf("bergm1_plot.pdf")
plot.horseBERGM(bergm1,lag=500)
dev.off()

pdf("horse1_plot.pdf")
plot.horseBERGM(horse1,lag=500)
dev.off()

pdf("bergm1_bgof.pdf")
set.seed(21911)
bgof(bergm1,sample.size=50, aux.iters=75, n.deg=10, n.dist=15, n.esp=6)
dev.off()

pdf("horse1_bgof.pdf")
set.seed(21912)
bgof(horse1,sample.size=50, aux.iters=75, n.deg=10, n.dist=15, n.esp=6)
dev.off()

save.image("sim_100.Rdata")
