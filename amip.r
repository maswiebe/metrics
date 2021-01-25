library(tidyverse)
library(repr)
library(AER)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sandwich)
library(zaminfluence)
library(devtools)
library(latex2exp)
library(lfe)

# installation:
# repo_loc <- Sys.getenv("/home/michael/Documents/git_repos/zaminfluence") 
# # Or just set to the correct directory path
# devtools::document(file.path("/home/michael/Documents/git_repos/zaminfluence", "zaminfluence"))
# devtools::install_local(file.path("/home/michael/Documents/git_repos/zaminfluence", "zaminfluence"), force=TRUE)
# 
# 
# InitializePython(file.path("/home/michael/Documents/git_repos/zaminfluence/", "venv/bin/python"))

# --------------

# activate venv:
# cd Documents/git_repos/zaminfluence/
# python3 -m venv venv
# source venv/bin/activate


base_dir  <- "/home/michael/Documents/git_repos/zaminfluence"
setwd(base_dir)
py_main <- InitializePython(file.path(base_dir, "venv/bin/python"))
reticulate::py_run_string("import regsens_rgiordandev")

#------------------------------------------
N=1000
K=20
#------------------------------------------
# simulations: true effect

amipfun = function(N,K,B,vX) {
  X <- rnorm(N,0,vX)
  e <- rnorm(N)
  z <- rnorm(N)
  y <- B*X + e
  
  reg_fit <- lm(y~X,x=T,y=T)
  
  reg_infl <- ComputeModelInfluence(reg_fit)
  grad_df <- GetTargetRegressorGrads(reg_infl, "X")
  influence_dfs <- SortAndAccumulate(grad_df)
  target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
  # rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
  # select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)
  
  sign <- target_change[target_change$change=="sign","prop_removed"]
  signif <- target_change[target_change$change=="significance","prop_removed"]
  sign_signif <- target_change[target_change$change=="sign and significance","prop_removed"]
  output <- c(signif,sign,sign_signif)
  # print(output)
  return(output)
}

set.seed(1515)
grid_b <- seq(0,0.2,0.025)
amip_b_signif <- c()
amip_b_sign <- c()
amip_b_sign_signif <- c()
for (b in grid_b) {
  print(b)
  sim <- replicate(n=nsim,amipfun(N,K,b,1))
  amip_b_signif <- rbind(amip_b_signif,mean(sim[1,],na.rm=T))
  amip_b_sign <- rbind(amip_b_sign,mean(sim[2,],na.rm=T))
  amip_b_sign_signif <- rbind(amip_b_sign_signif,mean(sim[3,],na.rm=T))
}

df1 <- data.frame(cbind(grid_b,amip_b_signif,amip_b_sign,amip_b_sign_signif))

#png(file="/home/michael/Dropbox/blog/amip/output/true_b.png",res=100)
ggplot(df1,aes(x=grid_b)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V3,color="b")) +
  geom_line(aes(y=V4,color="c")) +
  labs(title='True effects are robust',x=TeX('$\\beta$'),y="Proportion dropped") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(
    labels=unname(c("Significance","Sign","Both")),
    breaks=c("a","b","c")
    )
#dev.off()

# breaks down at 0.25
# when beta is large, rerun_df implictly returns NA: there exists no proportion that achieves the goal.
# but it actually drops that row entirely.
# which makes R change the datatype from double to list


#------------------------------------------
# simulations:gamma
# false positive driven by g*cov(X,z)

amipfun2 = function(N,K,vare,g) {
  X <- matrix(rnorm(N*K),nrow=N,ncol=K)
  # B=matrix(rnorm(N*K),ncol=K,nrow=N)
  e <- rnorm(N,0,vare)
  z <- rnorm(N)
  # y <- diag(X%*%t(B)) + g*z + e
  y <- g*z + e
  
  p_val_rob <- c()
  for (k in 1:K) {
    mod <- lm(y~X[,k])
    p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
  }
  # sum(p_val_rob<0.05)/K
  # minp <- min(p_val_rob)
  small <- which.min(p_val_rob)
  x1 <- X[,small]
  reg_fit <- lm(y~x1,x=T,y=T)
  
  reg_infl <- ComputeModelInfluence(reg_fit)
  grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
  influence_dfs <- SortAndAccumulate(grad_df)
  target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
  # rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
  # select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)
  
  sign <- target_change[target_change$change=="sign","prop_removed"]
  signif <- target_change[target_change$change=="significance","prop_removed"]
  sign_signif <- target_change[target_change$change=="sign and significance","prop_removed"]
  output <- c(signif,sign,sign_signif)
  return(output)
}

set.seed(1515)
grid_g <- seq(1,10,1)
amip_g_signif <- c()
amip_g_sign <- c()
amip_g_sign_signif <- c()
for (g in grid_g) {
  sim <- replicate(n=nsim,amipfun2(N,K,1,g))
  amip_g_signif <- rbind(amip_g_signif,mean(sim[1,]))
  amip_g_sign <- rbind(amip_g_sign,mean(sim[2,]))
  amip_g_sign_signif <- rbind(amip_g_sign_signif,mean(sim[3,]))
}

df2 <- data.frame(cbind(grid_g,amip_g_signif,amip_g_sign,amip_g_sign_signif))

#png(file="/home/michael/Dropbox/blog/amip/output/falsepos_g.png",res=100)
ggplot(df2,aes(x=grid_g)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V3,color="b")) +
  geom_line(aes(y=V4,color="c")) +
  labs(title='False positives are not',x=TeX('$\\gamma$'),y="Proportion dropped") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(
    labels=unname(c("Significance","Sign","Both")),
    breaks=c("a","b","c")
  ) +
  scale_y_continuous(breaks=seq(0,0.06,0.02),limits=c(0,0.06))
#dev.off()
