wd <- c("C:/Users/pamli/OneDrive/University of Louisville",
        "/My Topic for PhD Research/ERGMs","/R codes (P)","/2024")

setwd(paste0(wd[1],wd[2],wd[3]))
source("./Horseshoe functions/HorseBERGM_v3.R")
source("./Horseshoe functions/plot.horseBERGM.R")
source("./Horseshoe functions/summary.horseBERGM.R")
source("./Horseshoe functions/bgof_v2.R")
source("./bergm.R")

load("./2024/2024.10.08_Pamela/v2.Rdata")

library(mvtnorm)
library(statnet)
library(mcmc)
library(coda)
#-------------------------------------------------------------------------------
# Load and view Network 

gang_net <- readRDS("./2024/2024.10.08_Pamela/gang_network_isolates.RDS")
setwd("./2025.12.29_Pamela")

library(GGally)

btype <- get.vertex.attribute(gang_net, attrname = "birthplace")
v_color = ifelse(btype == 1, "West Africa", 
                 ifelse(btype == 2, "Jamaican",
                        ifelse(btype == 3,"UK",
                               "Somali")))

pdf("london_gang.pdf",width = 8, height = 5)
set.seed(127)
ggnet2(gang_net, color = v_color, size = 5,
       legend.size = 14, palette = c("West Africa" = "tomato",
                                     "Jamaican" = "gold", 
                                     "UK" = "steelblue", 
                                     "Somali" = "darkgreen"))
dev.off()

formula <- gang_net ~ edges + nodematch(~birthplace, diff = T) + triangle() + 
  triangle(~birthplace, diff = T, levels = 1:3) + isolates 


(sy <- summary(formula))   ;   dims <- length(sy)

main <- 15000  ;  burn <- 30000  ;  aux <- 1000  ;  pads <- 0.35  ;  chains <- 30   
p.mean <- rep(0,length(sy))     ;  prop <- 0.0005

#-------------------------------------------------------------------------------

mple <- ergm(formula,estimate = "MPLE")
summary(mple)

#-------------------------------------------------------------------------------

horse1 <- horseBERGM(formula,
                     prior.sigma = 100, 
                     lambda2 = 1, tau2 = 1,
                     main.iters = main, 
                     burn.in = burn, 
                     prior.mean = p.mean,
                     V.proposal = prop,
                     aux.iters = aux, 
                     gamma = pads,
                     nchain = chains,
                     seed = 123001)
summary.horseBERGM(horse1)
plot.horseBERGM(horse1,lag=200)

bergm1 <- bergm(formula,
                prior.sigma = diag(100,dims),
                main.iters = main, 
                burn.in = burn, 
                prior.mean = p.mean,
                V.proposal = prop,
                aux.iters = aux, 
                gamma = pads,
                nchain = chains,
                seed = 123002)

summary.horseBERGM(bergm1)
plot.horseBERGM(bergm1,lag=200)

save.image("v2.Rdata")

#-------------------------------------------------------------------------------

sink("london_results_v2.txt")
cat("# London Gang Network Model Comparison #\n")
cat("\n\n> formula <- gang_net ~ edges + nodematch(~birthplace, diff=T) + triangle() +\n",
    "  "," triangle(~birthplace, diff=T, levels=1:3) + isolates\n")

cat("\n# Summary of Horseshoe parallel sampling")
summary.horseBERGM(horse1)

cat("\n# ESS\n")
cbind('horse1'=horse1$ess, 'bergm1'=bergm1$ess)

cat("\n# Time (in minutes):\n")
horse1$Runtime
bergm1$Time


cat("\n\n# main <- 15000  ;  burn <- 30000  ;  aux <- 1000  ;  pads <- 0.35  ;  chains <- 30\n")
cat("\n# p.mean <- rep(0,length(sy))     ;  prop <- 0.0005")

sink()

#-------------------------------------------------------------------------------

plot.horseBERGM(horse1,lag=500)

pdf("horse_bgof2.pdf",width = 8, height = 6)
set.seed(123001)
bgof(horse1, sample.size = 200, aux.iters = 100, n.deg = 18, n.esp = 12, 
          n.dist = 12, n.triad = 4)

pdf("bergm_bgof1.pdf",width = 8, height = 6)
set.seed(123002)
bgof(bergm1, sample.size = 200, aux.iters = 100, n.deg = 18, n.esp = 12, 
     n.dist = 12, n.triad = 4)
dev.off()



