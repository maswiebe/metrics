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

png(file="/home/michael/Dropbox/blog/amip/output/true_b.png",res=100)
ggplot(df1,aes(x=grid_b)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V3,color="b")) +
  geom_line(aes(y=V4,color="c")) +
  labs(title='True effects are robust',x=TeX('$\\beta$'),y="Fraction dropped") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(
    labels=unname(c("Significance","Sign","Both")),
    breaks=c("a","b","c")
    )
dev.off()

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

png(file="/home/michael/Dropbox/blog/amip/output/falsepos_g.png",res=100)
ggplot(df2,aes(x=grid_g)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V3,color="b")) +
  geom_line(aes(y=V4,color="c")) +
  labs(title='False positives are not robust',x=TeX('$\\gamma$'),y="Fraction dropped") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(
    labels=unname(c("Significance","Sign","Both")),
    breaks=c("a","b","c")
  ) +
  scale_y_continuous(breaks=seq(0,0.06,0.02),limits=c(0,0.06))
dev.off()


  ### 
# all of the factors that affected whether controlling for z matters (var(b_i), var(e), gamma) don't matter for AMIP.







#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

N <- 1000
K <- 20

X <- matrix(rnorm(N*K),nrow=N,ncol=K)
B <- matrix(rnorm(N*K),ncol=N,nrow=K)
e <- rnorm(N)
z <- rnorm(N)
y <- diag(X%*%B) + z + e

p_val_rob <- c()
for (k in 1:K) {
  mod <- lm(y~X[,k])
  # mod <- lm(y~X[,k]+z)
  p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
}
sum(p_val_rob<0.05)/K

small <- which.min(p_val_rob)
x1 <- X[,small]
reg_fit <- lm(y~x1,x=T,y=T)
summary(reg_fit)
coeftest(reg_fit, vcov = vcovHC(reg_fit,type="HC0"))


# amip
reg_infl <- ComputeModelInfluence(reg_fit)
grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)
target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
if (FALSE) {
  PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)
}
rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)
rerun_df[2:4,"prop_removed"]
# sign; sign and significance; significance

# quite sensitive: only removes 0.002 fraction to kill significance
  # only 0.02 to kill sign
  # 0.05 to change sign and get significance


#---------------------------
# simulations: B~N(0,k)

amipfun = function(N,K,varb) {
  X <- matrix(rnorm(N*K),nrow=N,ncol=K)
  # for (i in 1:K){
  #   X[,i] = rnorm(N)
  # }
  B=matrix(rnorm(K*N,0,varb),ncol=K,nrow=N)
  # for (j in 1:K){
  #   B[,j] = rnorm(N,0,varb)
  # }
  e <- rnorm(N)
  z <- rnorm(N)
  y <- diag(X%*%t(B)) + z + e
  
  p_val_rob <- c()
  for (k in 1:K) {
    mod <- lm(y~X[,k])
    p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
  }
  # sum(p_val_rob<0.05)/K
  minp <- min(p_val_rob)
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
  output <- c(minp,signif,sign,sign_signif)
  print(output)
  return(output)
}

nsim <- 100
set.seed(1515)
grid_varb <- seq(1,10,1)
amip_varb_signif <- c()
amip_varb_sign <- c()
amip_varb_sign_signif <- c()
for (b in grid_varb) {
  sim <- replicate(n=nsim,amipfun(N,K,b))
  amip_varb_signif <- cbind(amip_varb_signif,mean(sim[2,]))
  amip_varb_sign <- cbind(amip_varb_sign,mean(sim[3,]))
  amip_varb_sign_signif <- cbind(amip_varb_sign_signif,mean(sim[4,]))
}

plot(grid_varb,amip_varb_signif)

df1 <- data.frame(cbind(grid_varb,t(amip_varb_signif),t(amip_varb_sign),t(amip_varb_sign_signif)))

ggplot(df1,aes(x=grid_varb)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V3,color="b")) +
  geom_line(aes(y=V4,color="c"))

pdf(file="/home/michael/Dropbox/blog/amip/output/false_varb.pdf")
ggplot(df1, aes(grid_varb)) +
  geom_line(aes(y=signif_varb,colour="Significance")) +
  labs(x=TeX('$Var(\\beta$)'),y="Proportion dropped") +
  geom_line(aes(y=sign_varb,colour="Sign")) +
  geom_line(aes(y=sign_signif_varb,colour="Sign and significance")) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(breaks=c("Significance","Sign","Sign and significance"))
dev.off()

# var(B) doesn't matter at all for amip
  # but it did matter when controlling for Z

hist(amipsim1[1,])

# significance
hist(amipsim1[2,])
mean(amipsim1[2,]) # 0.003 to kill significance
# sign
hist(amipsim1[3,])
mean(amipsim1[3,]) # 0.02 to change sign
# sign and significance
hist(amipsim1[4,])
mean(amipsim1[4,]) # 0.05 to change sign and and get significance

# plot proportion dropped against p-value
plot(amipsim1[1,],amipsim1[2,])
plot(amipsim1[1,],amipsim1[3,])
plot(amipsim1[1,],amipsim1[4,])
# need to drop more when p-value is smaller.


#------------------------------------------
# simulations: Var(e)
  # false positive driven by cov(X,e)

amipfun2 = function(N,K,vare) {
  X <- matrix(rnorm(N*K),nrow=N,ncol=K)
  B=matrix(rnorm(N*K),ncol=K,nrow=N)
  e <- rnorm(N,0,vare)
  z <- rnorm(N)
  y <- diag(X%*%t(B)) + z + e
  
  p_val_rob <- c()
  for (k in 1:K) {
    mod <- lm(y~X[,k])
    p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
  }
  # sum(p_val_rob<0.05)/K
  minp <- min(p_val_rob)
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
  output <- c(minp,signif,sign,sign_signif)
  return(output)
}


df2 <- data.frame(cbind(signif_vare,sign_vare,sign_signif_vare,grid_vare))

pdf(file="/home/michael/Dropbox/blog/amip/output/false_vare.pdf")
ggplot(df2, aes(grid_vare)) +
  geom_line(aes(y=signif_vare,colour="Significance")) +
  labs(x=TeX('$Var(\\epsilon$)'),y="Proportion dropped") +
  geom_line(aes(y=sign_vare,colour="Sign")) +
  geom_line(aes(y=sign_signif_vare,colour="Sign and significance")) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(breaks=c("Significance","Sign","Sign and significance"))
dev.off()

nsim <- 1000
set.seed(1515)
amipsim2 <- replicate(n=nsim,amipfun2(N,K))

mean(amipsim2[1,])
hist(amipsim2[1,])

# significance
hist(amipsim2[2,])
mean(amipsim2[2,]) # 0.003 to kill significance
# sign
hist(amipsim2[3,])
mean(amipsim2[3,]) # 0.02 to change sign
# sign and significance
hist(amipsim2[4,])
mean(amipsim2[4,]) # 0.06 to change sign and and get significance
# not that different from B~N(0,1)
  # source of false positive doesn't matter?

# plot proportion dropped against p-value
plot(amipsim2[1,],amipsim2[2,])
plot(amipsim2[1,],amipsim2[3,])
plot(amipsim2[1,],amipsim2[4,])



#------------------------------------------
# simulations: B=0, B~N(0,1) and X~N(0,5)
# false positive driven by cov(X,z)

amipfun4 = function(N,K) {
  X <- matrix(rnorm(N*K,0,5),nrow=N,ncol=K)
  # B=matrix(0,ncol=K,nrow=N)
  B=matrix(5,ncol=K,nrow=N)
  # for (j in 1:K){
  # B[,j] = rnorm(N)
  # }
  e <- rnorm(N)
  z <- rnorm(N)
  y <- diag(X%*%t(B)) + z + e
  
  p_val_rob <- c()
  for (k in 1:K) {
    mod <- lm(y~X[,k])
    p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
  }
  # sum(p_val_rob<0.05)/K
  minp <- min(p_val_rob)
  small <- which.min(p_val_rob)
  x1 <- X[,small]
  reg_fit <- lm(y~x1,x=T,y=T)
  
  reg_infl <- ComputeModelInfluence(reg_fit)
  grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
  influence_dfs <- SortAndAccumulate(grad_df)
  target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
  rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
  # select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)
  
  sign <- rerun_df[rerun_df$change=="sign","prop_removed"]
  signif <- rerun_df[rerun_df$change=="significance","prop_removed"]
  sign_signif <- rerun_df[rerun_df$change=="sign and significance","prop_removed"]
  output <- c(minp,signif,sign,sign_signif)
  return(output)
}

nsim <- 1000
set.seed(1515)
amipsim4 <- replicate(n=nsim,amipfun4(N,K))

mean(amipsim4[1,])
hist(amipsim4[1,])

# significance
hist(amipsim4[2,])
mean(amipsim4[2,]) # 0.003 to kill significance
# sign
hist(amipsim4[3,])
mean(amipsim4[3,]) # 0.02 to change sign
# sign and significance
hist(amipsim4[4,])
mean(amipsim4[4,]) # 0.06 to change sign and and get significance

# plot proportion dropped against p-value
plot(amipsim4[1,],amipsim4[2,])
plot(amipsim4[1,],amipsim4[3,])
plot(amipsim4[1,],amipsim4[4,])

# taking X~N(0,5) doesn't change anything
  # if average effect is B~=0, then Var(X) doesn't matter.



#------------------------------------------
# simulations: true effect

amipfun5 = function(N,K,B,vX) {
  X <- rnorm(N,0,vX)
  e <- rnorm(N)
  z <- rnorm(N)
  y <- B*X + z + e
  
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
  return(output)
}

set.seed(1515)
grid_b <- seq(0.1,0.3,0.025)
amip_b_signif <- c()
amip_b_sign <- c()
amip_b_sign_signif <- c()
for (b in grid_b) {
  sim <- replicate(n=nsim,amipfun5(N,K,b,1))
  amip_b_signif <- rbind(amip_b_signif,mean(sim[1,],na.rm=T))
  amip_b_sign <- rbind(amip_b_sign,mean(sim[2,],na.rm=T))
  amip_b_sign_signif <- rbind(amip_b_sign_signif,mean(sim[3,],na.rm=T))
}

df5 <- data.frame(cbind(grid_b,amip_b_signif,amip_b_sign,amip_b_sign_signif))

ggplot(df5,aes(x=grid_b)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V3,color="b")) +
  geom_line(aes(y=V4,color="c"))


# breaks down at 0.25
  # when beta is large, rerun_df implictly returns NA: there exists no proportion that achieves the goal.
  # but it actually drops that row entirely.
  # which makes R change the datatype from double to list


pdf(file="/home/michael/Dropbox/blog/amip/output/true_b.pdf")
ggplot(df5, aes(grid_b)) +
  geom_line(aes(y=signif_b,colour="Significance")) +
  labs(x=TeX('$\\beta$'),y="Proportion dropped") +
  geom_line(aes(y=sign_b,colour="Sign")) +
  geom_line(aes(y=sign_signif_b,colour="Sign and significance")) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(breaks=c("Significance","Sign","Sign and significance"))
dev.off()

### vary Var(X)


pdf(file="/home/michael/Dropbox/blog/amip/output/true_varx.pdf")
ggplot(df6, aes(grid_varx)) +
  geom_line(aes(y=signif_varx,colour="Significance")) +
  labs(x="Var(X)",y="Proportion dropped") +
  geom_line(aes(y=sign_varx,colour="Sign")) +
  geom_line(aes(y=sign_signif_varx,colour="Sign and significance")) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(breaks=c("Significance","Sign","Sign and significance"))
dev.off() 


X <- rnorm(N)
e <- rnorm(N)
z <- rnorm(N)
B <- 0.4
y <- B*X + z + e

reg_fit <- lm(y~X,x=T,y=T)

reg_infl <- ComputeModelInfluence(reg_fit)
grad_df <- GetTargetRegressorGrads(reg_infl, "X")
influence_dfs <- SortAndAccumulate(grad_df)
target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
select(rerun_df, change,beta, prop_removed)

sign <- rerun_df[rerun_df$change=="sign","prop_removed"]
signif <- rerun_df[rerun_df$change=="significance","prop_removed"]
sign_signif <- rerun_df[rerun_df$change=="sign and significance","prop_removed"]
output <- c(signif,sign,sign_signif)
output

# seems pretty nonlinear: have false-positive level robustness at B=0.1, but hit 50% for sign-and-significance at B=0.34.


# post
# true effects are robust
  # robustness increasing in B; higher B reduces var(e) -> higher signal-noise ratio
  # robustness increasing in var(X) -> higher signal-noise ratio
# false positives are not; at least, are just as robust as tiny true effects.
  # none of the factors that matter for controlling for z (var(b), var(e), gamma) matter for amip.
  # why? doesn't change number of influential observations?