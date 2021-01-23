require(broom)
require(ramify)
library("lmtest")
library("sandwich")
library(clubSandwich)
library("tidyverse")
library("latex2exp")
library(dplyr)
library(grf)
library(ri2)
library(lfe)
library(fastDummies)

#-------------------------------------------------------------------------------
# OLS, cross-section
#-------------------------------------------------------------------------------

raninf <- function(X,y){
  X_ri <- sample(X,length(X),replace=F)
  mod <- lm(y~X_ri)
  beta <- summary(mod)$coefficients[2,1]
  t <- summary(mod)$coefficients[2,3]
  t_rob <-coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,3]
  return(c(beta,t,t_rob))
} 

ri_fun_b0 = function(K,N,numsim) {
  y <- rnorm(N)
  X <- matrix(rnorm(N*K),nrow=N,ncol=K)
  p_val <- c()
  for (k in 1:K) {
    p_val <- c(p_val,summary(lm(y~X[,k]))$coefficients[2,4])
  }
  fracsig <- sum(p_val<0.05)/K
  X_min <- X[,which.min(p_val)]
  beta_min <- summary(lm(y~X_min))$coefficients[2,1]
  t_min <- summary(lm(y~X_min))$coefficients[2,3]
  p_val <- summary(lm(y~X_min))$coefficients[2,4]
  ri <- replicate(n=numsim,raninf(X_min,y))
  p_ri_b <- sum(abs(ri[1,])>abs(beta_min))/numsim
  p_ri_t <- sum(abs(ri[2,])>abs(t_min))/numsim
  
  pvals <- c(p_val,p_ri_b,p_ri_t,fracsig)
  print(pvals)
  return(pvals)
}

ri_fun_hetb = function(K,N,numsim) {
  X <- matrix(rnorm(N*K),nrow=N,ncol=K)
  B <- matrix(rnorm(N*K),ncol=K,nrow=N)
  e <- rnorm(N)
  y <- diag(X%*%t(B)) + e
  
  # p-hacking
  # using robust SE
  pvalrob <- c()
  for (k in 1:K) {
    loopmod <- lm(y~X[,k])
    pvalrob <- c(pvalrob,coeftest(loopmod, vcov = vcovHC(loopmod,type="HC0"))[2,4])
  }
  fracsig <- sum(pvalrob<0.05)/K
  X_min <- X[,which.min(pvalrob)]
  mod <- lm(y~X_min)
  beta_min <- summary(mod)$coefficients[2,1]
  p_val <- coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,4]
  t_min <- coeftest(mod, vcov = vcovHC(mod,type="HC0"))[2,3]
  
  # randomization inference
  # use robust SE for RI step as well
  ri <- replicate(n=numsim,raninf(X_min,y))
  p_ri_b <- sum(abs(ri[1,])>abs(beta_min))/numsim
  p_ri_t <- sum(abs(ri[3,])>abs(t_min))/numsim
  
  pvals <- c(p_val,p_ri_b,p_ri_t,fracsig)
  print(pvals)
  return(pvals)
}

# B=0
set.seed(1515)
pval_sims_ols <- replicate(n=nsims,ri_fun_b0(K,N,ri_sims))
sum(pval_sims_ols[1,]<0.05) # 653
sum(pval_sims_ols[1,]<0.05 & pval_sims_ols[2,]<0.05)/sum(pval_sims_ols[1,]<0.05) # 98%
sum(pval_sims_ols[1,]<0.05 & pval_sims_ols[3,]<0.05)/sum(pval_sims_ols[1,]<0.05) # 98%
sum(pval_sims_ols[1,]<0.05 & pval_sims_ols[2,]<0.05) # 638
sum(pval_sims_ols[1,]<0.05 & pval_sims_ols[3,]<0.05) # 638
# rejection rate:
mean(pval_sims_ols[4,]) # 0.051


# heterogeneous B
set.seed(1515)
output_ols_hetb <- replicate(n=nsims,ri_fun_hetb(K,N,ri_sims))

sum(output_ols_hetb[1,]<0.05) # 64.5%
sum(output_ols_hetb[1,]<0.05 & output_ols_hetb[2,]<0.05)/sum(output_ols_hetb[1,]<0.05) # 99%
sum(output_ols_hetb[1,]<0.05 & output_ols_hetb[3,]<0.05)/sum(output_ols_hetb[1,]<0.05) # 97%
sum(output_ols_hetb[1,]<0.05 & output_ols_hetb[2,]<0.05)
sum(output_ols_hetb[1,]<0.05 & output_ols_hetb[3,]<0.05)
# rejection rate
mean(output_ols_hetb[4,]) # 0.052

#-------------------------------------------------------------------------------
# DGP for diff in diff
# panel data, 50 states over 15 years; 10 treated states; 5 years pre-treatment data
# staggered rollout 
# treatments: set D=1 for all years after treatment year
# simple case: constant treatment effect (not varying over time)
#-------------------------------------------------------------------------------

N <- 50
Tn <- 10 # treatment group size
start_year <- 1995
end_year <- 2015
treatment_effect <- 0

#---
# constant B, differential trends
make_data_b0 <- function(N,Tn,start_year,end_year,treatment_effect,trend){
  timeline <- length(start_year:end_year)
  
  # state FEs
  state <- tibble(
    state = 1:N,
    state_fe = rnorm(N),
    mu=treatment_effect
  )
  
  # year FEs
  year <- tibble(
    year=start_year:end_year,
    year_fe = rnorm(timeline)
  )
  
  
  # treatment group
  treated <- tibble(
    state = sample(1:N,N,replace=F),
    treat_year = sort(c(sample(2000:2010,Tn,replace=T),rep(9999,(N-Tn)))),
    evertreat = ifelse(treat_year<9999,1,0)
  )
  
  # differential trend
  # trend <- 0.02
  # trend <- 0.2
  # trend <- 0.05
  
  # interact state X year
  expand_grid(state = 1:N, year = start_year:end_year) %>%
    left_join(state, by = "state") %>%
    left_join(year, by = "year") %>%
    left_join(treated, by = "state") %>%
    mutate(error = rnorm(N*timeline), # dimension NTx1
           treat = ifelse(year>=treat_year,1,0),
           tau = ifelse(treat==1,mu,0)) %>%
    mutate(y = state_fe + year_fe + tau + (year-start_year)*trend*evertreat + error)
}


#---
# heterogeneous B~N(0,1)
make_data_bhet <- function(N,K,Tn,start_year,end_year,het){
  timeline <- length(start_year:end_year)
  
  # state FEs
  state <- tibble(
    state = 1:N,
    state_fe = rnorm(N)
  )
  
  # year FEs
  year <- tibble(
    year=start_year:end_year,
    year_fe = rnorm(timeline)
  )
  
  paneldata <- expand_grid(state=1:N, year=start_year:end_year)
  # create K treatment variables
  for (k in 1:K) {
    treated <- tibble(
      state = sample(1:N,N,replace=F),
      treat_year = sort(c(sample(2000:2010,Tn,replace=T),rep(9999,(N-Tn))))
    )
    te <- tibble(
      state = 1:N, 
      mu = rnorm(N) # draw mu_k (Nx1) within the loop, for each k.
    )
    
    looptau <- expand_grid(state=1:N, year=start_year:end_year) %>%
      left_join(treated, by = "state") %>%
      left_join(te, by = "state") %>%
      mutate(treat = ifelse(year>=treat_year,1,0),
             tau = ifelse(treat==1,mu,0)) %>%
      # need to use !! to unquote, then := as assignment operator
      rename(!!paste0("treat_year",k) := treat_year,
             !!paste0("treat",k) := treat,
             !!paste0("tau",k) := tau) %>%
      select(-mu)
    # if creating variable names with paste, then can't use them in operations
    
    paneldata <- paneldata %>%
      left_join(looptau, by = c("state","year")) #%>%
  }
  
  paneldata %>%
    select(starts_with("tau")) %>%
    rowSums(.) -> paneldata$sumtau
  
  paneldata %>%
    left_join(state, by = c("state")) %>%
    left_join(year, by = c("year")) %>%
    mutate(error = rnorm(N*timeline), # dimension NTx1
           y = state_fe + year_fe + het*sumtau + error)
}

#------------------------------------
# randomization inference methods:
# (1) keep treatment years, randomly assign treated states
# (2) randomly assign both treatment years and states
# natural approach; and stronger

# shuffle treat_year, then generate treat_ri based on treat_year_ri
  # randomizing treatment status only, not treatment year
ri_fun <- function(datasmall,N){
  data_collapse <- datasmall %>%
    select(state,treat_year) %>%
    group_by(state) %>%
    summarize(treat_year=mean(treat_year),.groups="drop_last") %>%
    mutate(treat_year_ri = sample(treat_year,N,replace=F))  
  # this is randomizing treatment status, but keeping the same treatment years
  # this is what Gavrilova did
  data_ri <- datasmall %>%
    left_join(data_collapse,by=c("state","treat_year")) %>%
    mutate(treat_ri = ifelse(year>=treat_year_ri,1,0))
  
  mod <- felm(y~treat_ri | state + year |0| state ,data=data_ri)
  beta <- summary(mod)$coefficients[1,1]
  t <- summary(mod)$coefficients[1,3]
  return(c(beta,t))
}

# randomize treatment assignment and treatment year
ri_fun2 <- function(datasmall,N){
  treated <- tibble(
    state = sample(1:N,N,replace=F),
    treat_year_ri = sort(c(sample(2000:2010,Tn,replace=T),rep(9999,(N-Tn))))
  )
  
  data_ri <- datasmall %>%
    left_join(treated,by=c("state")) %>%
    mutate(treat_ri = ifelse(year>=treat_year_ri,1,0))
  
  mod <- felm(y~treat_ri | state + year |0| state ,data=data_ri)
  beta <- summary(mod)$coefficients[1,1]
  t <- summary(mod)$coefficients[1,3]
  return(c(beta,t))
}

#--------------------  
# simulations

ri_simfun = function(numsim,het) {
  N <- 50
  # N <- 500
  Tn <- 10
  # Tn <- 250
  start_year <- 1995
  end_year <- 2015
  
  data <- make_data_bhet(N,K,Tn,start_year,end_year,het)
  
  bhet_regressors <- data %>%
    select(starts_with("treat")) %>%
    select(-starts_with("treat_year"))
  pval <- c()
  for (i in 1:K) {
    loopmod <- felm(y~bhet_regressors[[i]] | state + year |0| state, data=data)
    pval <- c(pval,summary(loopmod)$coefficients[1,4])
  }
  psmall <- which.min(pval)
  fracsig <- sum(pval<0.05)/K
  
  data_min <- data %>%
    select(state,year,
           treat=as.name(paste0("treat",psmall)),
           treat_year=as.name(paste0("treat_year",psmall)),
           y)
  
  mod <- felm(y~treat | state + year |0| state, data=data_min)
  
  beta_min <- summary(mod)$coefficients[1,1]
  pval_min <- summary(mod)$coefficients[1,4]
  t_min <- summary(mod)$coefficients[1,3]
  
  # ri_did <- replicate(n=numsim,ri_fun(data_min,N))
  ri_did <- replicate(n=numsim,ri_fun2(data_min,N))
  
  pval_ri_b <- sum(abs(ri_did[1,])>abs(beta_min))/numsim
  pval_ri_t <- sum(abs(ri_did[2,])>abs(t_min))/numsim
  
  pvals <- c(pval_min,pval_ri_b,pval_ri_t,fracsig)
  print(pvals)
  return(pvals)
}

#------------
# beta=0

nsim <- 1000
ri_sims <- 1000

set.seed(1515)
output_dd_b0b <- replicate(n=nsim,ri_simfun(ri_sims,0))
sum(output_dd_b0b[1,]<0.05)/nsim # OLS: 82%

sum(output_dd_b0b[1,]<0.05 & output_dd_b0b[2,]<0.05)
sum(output_dd_b0b[1,]<0.05 & output_dd_b0b[2,]<0.05)/sum(output_dd_b0b[1,]<0.05) # RI-b: 65%

sum(output_dd_b0b[1,]<0.05 & output_dd_b0b[3,]<0.05)
sum(output_dd_b0b[1,]<0.05 & output_dd_b0b[3,]<0.05)/sum(output_dd_b0b[1,]<0.05) # RI-t: 80%
mean(output_dd_b0b[4,]) # rejection rate: 0.079

### is rejection rate based on small treatment group?
# try N=500, Tn=250
output_dd_b0 <- replicate(n=1000,ri_simfun(ri_sims,1))
mean(output_dd_b0[4,]) # 0.05
sum(output_dd_b0[1,]<0.05) # 64%

#-------------
### beta ~ N(0,1)

set.seed(1515)
output_dd_bhet2 <- replicate(n=nsim,ri_simfun(ri_sims,1))

sum(output_dd_bhet2[1,]<0.05)/nsim # OLS: 79%
mean(output_dd_bhet2[4,]) # rejection rate: 0.077

sum(output_dd_bhet2[1,]<0.05 & output_dd_bhet2[2,]<0.05)
sum(output_dd_bhet2[1,]<0.05 & output_dd_bhet2[2,]<0.05)/sum(output_dd_bhet2[1,]<0.05) # RI-b: 85%

sum(output_dd_bhet2[1,]<0.05 & output_dd_bhet2[3,]<0.05)
sum(output_dd_bhet2[1,]<0.05 & output_dd_bhet2[3,]<0.05)/sum(output_dd_bhet2[1,]<0.05) # RI-t: 87%



#------------------------------------------
#------------------------------------------
#------------------------------------------
# is RI a placebo/robustness test with differential trends?

ri_simfun_trend <- function(sims,trend) {
  N <- 50
  Tn <- 10
  start_year <- 1995
  end_year <- 2015
  trenddata <- make_data_b0(N,Tn,start_year,end_year,0,trend)
  tmod <- felm(y~treat | state + year |0| state, data=trenddata)
  
  btrend <- summary(tmod)$coefficients[1,1]
  ptrend <- summary(tmod)$coefficients[1,4]
  ttrend <- summary(tmod)$coefficients[1,3]
  
  ri_trend <- replicate(n=sims,ri_fun(trenddata,N))
  p_ri_b <- sum(abs(ri_trend[1,])>abs(btrend))/ri_sims
  p_ri_t <- sum(abs(ri_trend[2,])>abs(ttrend))/ri_sims
  
  output <- c(ptrend,p_ri_b,p_ri_t)
  print(output)
  return(output)
}

sims <- 100
trend <- 0.1
trend_test <- replicate(sims,ri_simfun_trend(ri_sims,trend))

# rejection rate
sum(trend_test[1,]<0.05)/sims 
sum(trend_test[2,]<0.05)/sims 
sum(trend_test[3,]<0.05)/sims 
# get 100% rejection rate when trend=0.1, for all p-values

# fraction remaining significant using RI
sum(trend_test[1,]<0.05 & trend_test[2,]<0.05)/sum(trend_test[1,]<0.05) # 1
sum(trend_test[1,]<0.05 & trend_test[3,]<0.05)/sum(trend_test[1,]<0.05) # 1

trend_fun <- function(sims,trend) {
  sim <- replicate(sims,ri_simfun_trend(ri_sims,trend))  
  rej1 <- sum(sim[1,]<0.05)/sims
  rej2 <- sum(sim[2,]<0.05)/sims
  rej3 <- sum(sim[3,]<0.05)/sims
  p1 <- mean(sim[1,])
  p2 <- mean(sim[2,])
  p3 <- mean(sim[3,])
  output <- c(rej1,rej2,rej3,p1,p2,p3)
  print(output)
  return(output)
}


trendgrid <- seq(0,0.1,0.02)
rejrate_trend <- c()
for (i in trendgrid){
  rejrate_trend <- cbind(rejrate_trend,trend_fun(sims,i))
  print(i)
}
plot(trendgrid,rejrate_trend[1,])
plot(trendgrid,rejrate_trend[2,])
plot(trendgrid,rejrate_trend[3,])

df_trend <- data.frame(cbind(trendgrid,t(rejrate_trend)))

ggplot(df_trend,aes(trendgrid)) +
  geom_line(aes(y=V2,color="a")) +
  geom_line(aes(y=V4,color="b")) +
  geom_line(aes(y=V5,color="c")) +
  geom_line(aes(y=V7,color="d")) +
  labs(title=TeX('Comparing p and p_{RI} under differential trends'),x=TeX('Trend $(\\lambda)$'),y="") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "bottom") +
  scale_color_discrete(
    labels=unname(TeX(c('Rejection $p$','Rejection $p_{RI}$','Average p','Average p_{RI}'))),
    breaks=c("a","b","c","d")
  ) +
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(breaks=seq(0,0.1,0.02),limits=c(0,0.1))

ptest <- function(n,bvar) {
  # x <- rnorm(n)
  x <- rbernoulli(n,p=0.5)
  e <- rnorm(n)
  b <- rlnorm(n,0,bvar)
  y <- b*x + e
  model <-summary(lm(y~x)) 
  p1 <- model$coefficients[2,4]
  t <- model$coefficients[2,3]
  simri <- replicate(n=ri_sims,raninf(x,y))
  p2 <- sum(abs(simri[3,])>abs(t))/ri_sims
  diff <- abs(p1-p2)
  print(diff)
  return(diff)
}

#---------------------------------
#---------------------------------
#---------------------------------
# what does it take to make get p_RI different?

ptest <- function(n,bvar) {
  # x <- rnorm(n)
  x <- rbernoulli(n,p=0.5)
  e <- rnorm(n)
  b <- rlnorm(n,0,bvar)
  y <- b*x + e
  model <-summary(lm(y~x)) 
  p1 <- model$coefficients[2,4]
  t <- model$coefficients[2,3]
  simri <- replicate(n=ri_sims,raninf(x,y))
  p2 <- sum(abs(simri[3,])>abs(t))/ri_sims
  diff <- abs(p1-p2)
  print(diff)
  return(diff)
}

# graph, X~B(p=0.5)
bvar_seq <- seq(1,10,1)
diff_fvb <- c()
for (b in bvar_seq) {
  differ <- mean(replicate(n=100,ptest(20,b)))
  diff_fvb <- c(diff_fvb,differ)
}
plot(bvar_seq,diff_fvb)

df_test <- data.frame(cbind(bvar_seq,diff_fvb))

ggplot(df_test, aes(bvar_seq)) +
  geom_line(aes(y=diff_fvb)) +
  labs(x=TeX('$Var(\\beta$)'),y="Difference") +
  scale_x_continuous(breaks=seq(0,10,2),limits=c(0,10))

