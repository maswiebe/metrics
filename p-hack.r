library(broom)
library(ramify)
library(lmtest)
library(sandwich)
library(tidyverse)
library(latex2exp)
library(dplyr)
library(grf)
library(ri2)
library(asbio)
library(sur)
#----------------------------------------------------------------------------
# homogeneous B=0
#----------------------------------------------------------------------------
K <- 20
N <- 1000

simfun1 = function(J,num,gamma,vare) {
  B <- c(1:J)*0
  X <- matrix(rnorm(num*J),nrow=num,ncol=J)
  e <- rnorm(num,0,vare)
  z <- rnorm(num)
  y <- X%*%B + gamma*z + e
  
  p_val <- c()
  for (k in 1:J) {
    mod <- lm(y~X[,k])
    p_val <- c(p_val,summary(mod)$coefficients[2,4])
  }
  frac <- sum(p_val<0.05)/J
  small <- which.min(p_val)
  Xmin <- X[,small]
  ratio_xe <- abs(cov(Xmin,e))/max(abs(cov(X,e)))
  ratio_xz <- abs(cov(Xmin,z))/max(abs(cov(X,z)))
  mod1 <- lm(y~Xmin)
  mod2 <- lm(y~Xmin + z)
  p1 <- summary(mod1)$coefficients[2,4]
  b1 <- summary(mod1)$coefficients[2,1]
  p2 <- summary(mod2)$coefficients[2,4]
  b2 <- summary(mod2)$coefficients[2,1]
  covxz <- cov(Xmin,z)
  covxe <- cov(Xmin,e)
  pr2 <- partial.R2(mod1,mod2)
  var_res <- var(resid(lm(Xmin~z)))
  covze <- cov(z,e)
  se2 <- summary(mod2)$coefficients[2,2]
  t2 <- summary(mod2)$coefficients[2,3]
  reg1r2 <- summary(mod1)$r.squared
  reg2r2 <- summary(mod2)$r.squared
  test <- c(p1,p2,covxz,frac,pr2,covxe,b1,b2,var_res,covze,se2,t2,ratio_xe,ratio_xz,reg2r2,reg1r2)
}

nsim <- 1000
set.seed(1515)
sims1 <- replicate(n=nsim,simfun1(K,N,1,1))

mean(sims1[4,])

sig1a <- (sims1[1,]<0.05)
sig1b <- (sims1[2,]<0.05)
sum(sig1a) # 663
sum(sig1a==1 & sig1b==1) # 245
sum(sig1a==1 & sig1b==0) # 418
sum(sig1a==1 & sig1b==1)/sum(sig1a) # 0.37 of significant b2's have a significant b1
sum(sig1a==1 & sig1b==0)/sum(sig1a) # 0.63 lose significance

#-------------------------------------------------------
fun1 <- function(K,N,gamma,vare) {
  set.seed(1515)
  sim <- replicate(n=nsim,simfun1(K,N,gamma,vare))
  sig1 <- (sim[1,]<0.05)
  sig2 <- (sim[2,]<0.05)
  still_alive <- sum(sig1==1 & sig2==1)/sum(sig1)
  still_alive2 <- sum(sig1==1 & sig2==1)/sum(sig2)
  cov_xz <- median(gamma*sim[3,]/(sim[6,] + gamma*sim[3,])) # need to use median, since denom can =0
  cov_xe <- median(sim[6,]/(sim[6,] + gamma*sim[3,]))
  pr2 <- mean(sim[5,])
  reg2r2 <- mean(sim[15,])
  reg1r2 <- mean(sim[16,])
  return(c(cov_xz,cov_xe,pr2,still_alive,reg2r2,reg1r2,still_alive2))
}

grid <- seq(0,5,0.25)
output1 <- c()
for (i in grid){
  output1 <- cbind(output1,fun1(K,N,i,1))
  print(i)
}

df1 <- data.frame(cbind(grid,t(output1))) %>%
  mutate(oms = 1-V5)

# pdf(file="/home/michael/Dropbox/blog/p-hack/output/b0_shares_signif.pdf")
#png(file="/home/michael/Dropbox/blog/p-hack/output/b0_shares_signif.png")
ggplot(df1, aes(grid)) +
  geom_line(aes(y=V2,color='a'))+
  geom_line(aes(y=V3,color='b')) +
  geom_line(aes(y=oms,color='c'))+
  geom_line(aes(y=V4,color="d"))+
  # geom_line(aes(y=V7,color="e"))+
  labs(x=TeX('$\\gamma$'),y="") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_color_discrete(
    labels=unname(TeX(c('$\\gamma Cov(X,z)$','$Cov(X,\\epsilon)$',"Lose significance",'Partial R^{2}(z)'))),
    breaks=c("a","b","c","d")
  ) +
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))
#dev.off()

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# heterogeneous B~N(0,1)
#----------------------------------------------------------------------------
set.seed(1518)
X <- matrix(rnorm(N*K), ncol=K,nrow=N)
B <- matrix(rnorm(N*K),ncol=K,nrow=N)
e <- rnorm(N)
z <- rnorm(N)
g <- 1
y <- diag(X%*%t(B)) + g*z + e
p_val_rob <- c()
for (k in 1:K) {
  mod <- lm(y~X[,k])
  p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
}
small <- which.min(p_val_rob)
xmin <- X[,small]
mod1 <- lm(y~xmin)
mod2 <- lm(y~xmin+z)
summary(mod1)
summary(mod2)

# output regression table
library(stargazer)
stargazer(mod1,mod2,
          covariate.labels = c("X","z"),
          omit = c("Constant"),
          omit.stat = c("f","ser","rsq"),
          dep.var.caption = "",
          dep.var.labels.include = F,
          no.space=T,
          out="/home/michael/Dropbox/blog/p-hack_robust/output/het_table.tex"
)

simfun2 = function(K,N,bvar,g,ve) {
  X <- matrix(rnorm(N*K), ncol=K,nrow=N)
  B <- matrix(rnorm(N*K,0,bvar),ncol=K,nrow=N)
  e <- rnorm(N,0,ve)
  z <- rnorm(N)
  y <- diag(X%*%t(B)) + g*z + e
  
  p_val_rob <- c()
  for (k in 1:K) {
    mod <- lm(y~X[,k])
    # mod <- lm(y~X[,k]+z)
    p_val_rob <- c(p_val_rob,coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4])
  }
  frac <- sum(p_val_rob<0.05)/K
  small <- which.min(p_val_rob)
  xmin <- X[,small]
  mod1 <- lm(y~xmin)
  mod2 <- lm(y~xmin+z)
  p1 <- coeftest(mod1, vcov = vcovHC(mod1,type="HC0"))[2,4]
  b1 <- summary(mod1)$coefficients[2,1]
  p2 <- coeftest(mod2, vcov = vcovHC(mod2,type="HC0"))[2,4]
  b2 <- summary(mod2)$coefficients[2,1]
  cov_xz <- cov(xmin,z)
  cov_xe <- cov(xmin,e)
  cov_xjxk <- c()
  for (j in 1:K){
    cov_xjxk <- c(cov_xjxk,cov(xmin,X[,j]*B[,j]))
  }
  sumcov_xjxk <- sum(cov_xjxk) # terms could cancel out; don't want abs
  cov_xbx <- cov(xmin,xmin*B[,small])
  ratio_xe <- abs(cov(xmin,e))/max(abs(cov(X,e)))
  ratio_xz <- abs(cov(xmin,z))/max(abs(cov(X,z)))
  pr2 <- partial.R2(mod1,mod2)
  var_res <- var(resid(lm(xmin~z)))
  test <- c(p1,p2,frac,cov_xz,cov_xe,b1,b2,ratio_xe,ratio_xz,pr2,var_res,cov_xbx,sumcov_xjxk)
  return(test)
}

set.seed(1515)
sims2 <- replicate(n=nsim,simfun2(K,N,1,1,1))

sig2a <- (sims2[1,]<0.05)
sig2b <- (sims2[2,]<0.05)
sum(sig2a) # 650
sum(sig2a==1 & sig2b==1) # 569
sum(sig2a==1 & sig2b==0) # 81
sum(sig2a==1 & sig2b==1)/sum(sig2a) # 0.88 of significant b2's have a significant b1


fun2 <- function(K,N,varb,gamma,vare) {
  set.seed(1515)
  sim <- replicate(n=nsim,simfun2(K,N,varb,gamma,vare))
  sig1 <- (sim[1,]<0.05)
  sig2 <- (sim[2,]<0.05)
  still_alive <- sum(sig1==1 & sig2==1)/sum(sig1)
  sumcov_xjxk <- median(sim[13,]/(sim[5,] + gamma*sim[4,] + sim[13,]))
  cov_xz <- median(gamma*sim[4,]/(sim[5,] + gamma*sim[4,] + sim[13,]))
  cov_xe <- median(sim[5,]/(sim[5,] + gamma*sim[4,] + sim[13,]))
  pr2 <- mean(sim[10,])
  return(c(sumcov_xjxk,cov_xz,cov_xe,still_alive,pr2))
}

output2 <- c()
for (i in grid){
  output2 <- cbind(output2,fun2(K,N,1,i,1))
  print(i)
}

# graph
df2 <- data.frame(cbind(grid,t(output2))) %>%
  mutate(oms = 1-V5)

# plot shares of hat_alpha_1: cov(XkBk,Xk), sum cov(XjBj, Xk), cov(x,z), cov(x,e)
# pdf(file="/home/michael/Dropbox/blog/p-hack/output/bhet_shares.pdf")
#png(file="/home/michael/Dropbox/blog/p-hack/output/bhet_shares.png")
ggplot(df2, aes(grid)) +
  geom_line(aes(y=V2,color="a"))+
  geom_line(aes(y=V3,color="b"))+
  geom_line(aes(y=V4,color="c"))+
  geom_line(aes(y=oms,color="d"))+
  geom_line(aes(y=V6,color="e"))+
  labs(x=TeX('$\\gamma$'),y="") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_color_discrete(
    labels=unname(TeX(c("$\\sum Cov(X_{k},\\beta_{j} X_{j})$",'$\\gamma*Cov(X_{k},z)$',"$Cov(X_{k},\\epsilon)$","Lose significance",'Partial R^{2}(z)'))),
    breaks=c("a","b","c","d","e")
  ) +
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))
#dev.off()
