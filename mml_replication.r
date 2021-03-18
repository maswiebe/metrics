library(tidyverse)
library(tabulator)
library(bacondecomp)
library(lfe)
library(dplyr)
library(fastDummies)
library(specr)
library(fixest)

mutate <- dplyr::mutate
distinct <- dplyr::distinct
select <- dplyr::select

ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}

#----------------------------------------------------------------------------

mml_full <- haven::read_dta("/home/michael/Dropbox/replications/gavrilova-mml/zoutman_data/GKZ_MML_2017.dta")

mml_full$treat_year <- NA
mml_full$treat_year[mml_full$state_name=="Arizona"] <- 2010
mml_full$treat_year[mml_full$state_name=="California"] <- 1996
mml_full$treat_year[mml_full$state_name=="Colorado"] <- 2001
mml_full$treat_year[mml_full$state_name=="Connecticut"] <- 2012
mml_full$treat_year[mml_full$state_name=="District of Columbia"] <- 2010
mml_full$treat_year[mml_full$state_name=="Delaware"] <- 2011
mml_full$treat_year[mml_full$state_name=="Maine"] <- 1999
mml_full$treat_year[mml_full$state_name=="Michigan"] <- 2008
mml_full$treat_year[mml_full$state_name=="Montana"] <- 2004
mml_full$treat_year[mml_full$state_name=="Nevada"] <- 2001
mml_full$treat_year[mml_full$state_name=="New Jersey"] <- 2010
mml_full$treat_year[mml_full$state_name=="New Mexico"] <- 2007
mml_full$treat_year[mml_full$state_name=="Oregon"] <- 1998
mml_full$treat_year[mml_full$state_name=="Rhode Island"] <- 2006
mml_full$treat_year[mml_full$state_name=="Vermont"] <- 2004
mml_full$treat_year[mml_full$state_name=="Washington"] <- 1998
mml_full$treat_year[is.na(mml_full$treat_year)] <- 9999

mml_full <- mml_full %>%
  filter(year>=1994) %>%
  mutate(mmlXborder = mml*border,
         lrate = log(violent_rate+1),
         lhom = log(murder_reported_rate+1),
         lrob = log(robbery_reported_rate+1),
         lass = log(assault_reported_rate+1)
  )

mml <- mml_full %>%
  filter(year>=1994) %>%
  filter(black_rate<=1) %>%
  mutate(mmlXborder = mml*border,
         lrate = log(violent_rate+1),
         lhom = log(murder_reported_rate+1),
         lrob = log(robbery_reported_rate+1),
         lass = log(assault_reported_rate+1),
         ihs_rate = ihs(violent_rate),
         ihs_hom = ihs(murder_reported_rate),
         ihs_rob = ihs(robbery_reported_rate),
         ihs_ass = ihs(assault_reported_rate),
  )
# GKZ use crimes per 100,000 population; so log(rate) is well-defined

covars <- "lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec"


#----------------------------------------------------------------------------
# main regressions
feols(murder_reported_rate~mml + mmlXborder +lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)

feols(robbery_reported_rate~mml + mmlXborder +lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)

feols(assault_reported_rate~mml + mmlXborder +lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)


# heterogeneity by border state, for disaggregated outcome variables
  # paper doesn't do this
feols(murder_reported_rate~ mml_m_ar + mml_m_cali + mml_m_nm + mml_nb_c +lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)
# driven by cali, arizona; nothing for NM

feols(robbery_reported_rate~ mml_m_ar + mml_m_cali + mml_m_nm + mml_nb_c +lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)
# cali and Ariz

feols(assault_reported_rate~ mml_m_ar + mml_m_cali + mml_m_nm + mml_nb_c +lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)
# Cali and NM

#----------------------------------------------------------------------------
# regression weights

# diff in diff
mod <- as.formula(paste0("mml~",covars,"| factor(state)[year] +  factor(year)[border] + fips"))
tibble(as.data.frame(residuals(feols(mod,cluster=~state,weights=~population_weight,data=mml))))

resid = residuals(feols(mod,cluster=~state,weights=~population_weight,data=mml))
res2 = resid^2
restot = sum(res2)
regw = res2/restot

mml <- mml %>%
  mutate(
    regw=regw
  )

# mml <- mml %>%
mml %>%
  dplyr::mutate(
    # resid = residuals(feols(mml~lnmedian_income + unemploymentrate + male_rate + black_rate + hisp_rate + age2024_rate + age1019_rate + lnpop + poverty_all_percent + dec| factor(state)[year] +  factor(year)[border] + fips,cluster=~state,weights=~population_weight ,data=mml)),
    resid = residuals(feols(mod,cluster=~state,weights=~population_weight,data=mml)),
    # resid =1,
    res2 = resid^2,
    restot = sum(res2),
    regw = res2/restot
  )

mml %>% 
  group_by(state) %>%
  summarize(rw = round(sum(regw)*100)) %>%
  tab(rw)
  # View()

# 26 michigan: 16%
# 48 texas: 12%
# 8 colorado: 12%
# 35 new mexico: 10%
# 30 montana: 10%

# 6 california: 2%
# 4 arizona: 4%

# so in DD, only 16% of weight is on treated border states

#----------------------------------------------------------------------------
# multiverse analysis
# big robustness graph

ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}

#--------------------------------------------------
# models
baseline <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml", "|fips+year"))
  feols(formula,cluster=~state,data)
}
trend <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|fips + factor(state)[year]"))
  feols(formula,cluster=~state,data)
}
weights <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips+year"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}
border <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips + factor(year)[border]"))
  feols(formula,cluster=~state,data)
}
trend_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|fips + factor(state)[year]"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}
trend_border <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|fips + factor(state)[year] + factor(year)[border]"))
  feols(formula,cluster=~state,data)
}
border_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips + factor(year)[border]"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}
trend_border_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips + factor(state)[year] + factor(year)[border]"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}

pbaseline <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml", "|fips+year"))
  fepois(formula,cluster=~state,data)
}
ptrend <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|fips + factor(state)[year]"))
  fepois(formula,cluster=~state,data)
}
pweights <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips+year"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}
pborder <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips + factor(year)[border]"))
  fepois(formula,cluster=~state,data)
}
ptrend_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|fips + factor(state)[year]"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}
ptrend_border <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|fips + factor(state)[year] + factor(year)[border]"))
  fepois(formula,cluster=~state,data)
}
pborder_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips + factor(year)[border]"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}
ptrend_border_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | fips + factor(state)[year] + factor(year)[border]"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}


#-----------------------
### violent_rate

# level
results_viol <- run_specs(df = mml, 
                          y = c("violent_rate"),
                          x = c("mmlXborder"),
                          model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                          controls = c(covars) 
)
results_viol <- results_viol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )

png(file="/home/michael/Dropbox/blog/mml/output/violent_level.png",res=100)
# png(file="/home/michael/Dropbox/blog/mml/output/violent_level1.png")
plot_specs(results_viol,
           choices = c("model","controls"))
dev.off()

# log(y+1)
results_lviol <- run_specs(df = mml, 
                           y = c("lrate"),
                           x = c("mmlXborder"),
                           model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                           controls = c(covars) 
)
results_lviol <- results_lviol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/violent_log.png",res=100)
plot_specs(results_lviol,
           choices = c("model","controls"))
dev.off()

# inverse hyperbolic sine
results_ihsviol <- run_specs(df = mml, 
                             y = c("ihs_rate"),
                             x = c("mmlXborder"),
                             model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                             controls = c(covars) 
)
results_ihsviol <- results_ihsviol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/violent_ihs.png",res=100)
plot_specs(results_ihsviol,
           choices = c("model","controls"))
dev.off()

# poisson
# poisson has to use count data, instead of rate? don't have count variable here.
  # already weighting by population, and controlling for log(pop)
results_pviol <- run_specs(df = mml, 
                           y = c("violent_rate"),
                           x = c("mmlXborder"),
                           model = c("pbaseline","ptrend","pweights","pborder","ptrend_weight","ptrend_border","pborder_weight","ptrend_border_weight"),
                           controls = c(covars) 
)
results_pviol <- results_pviol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="ptrend_weight","trends + weights"),
         model = replace(model, model=="ptrend_border","trends + border"),
         model = replace(model, model=="pborder_weight","border + weights"),
         model = replace(model, model=="ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="pweights","weights"),
         model = replace(model, model=="ptrend","trends"),
         model = replace(model, model=="pborder","border"),
         model = replace(model, model=="pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/violent_pois.png",res=100)
plot_specs(results_pviol,
           choices = c("model","controls"))
dev.off()

#-----------------------
### homicides

# level
results_hom <- run_specs(df = mml, 
                         y = c("murder_reported_rate"),
                         x = c("mmlXborder"),
                         model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                         controls = c(covars) 
)
results_hom <- results_hom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/hom_level.png",res=100)
plot_specs(results_hom,
           choices = c("model","controls"))
dev.off()

# log(y+1)
results_lhom <- run_specs(df = mml, 
                          y = c("lhom"),
                          x = c("mmlXborder"),
                          model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                          controls = c(covars) 
)
results_lhom <- results_lhom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/hom_log.png",res=100)
plot_specs(results_lhom,
           choices = c("model","controls"))
dev.off()

# inverse hyperbolic sine
results_ihshom <- run_specs(df = mml, 
                            y = c("ihs_hom"),
                            x = c("mmlXborder"),
                            model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                            controls = c(covars) 
)
results_ihshom <- results_ihshom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/hom_ihs.png",res=100)
plot_specs(results_ihshom,
           choices = c("model","controls"))
dev.off()

# poisson
results_phom <- run_specs(df = mml, 
                          y = c("murder_reported"),
                          x = c("mmlXborder"),
                          model = c("pbaseline","ptrend","pweights","pborder","ptrend_weight","ptrend_border","pborder_weight","ptrend_border_weight"),
                          controls = c(covars) 
)
results_phom <- results_phom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="ptrend_weight","trends + weights"),
         model = replace(model, model=="ptrend_border","trends + border"),
         model = replace(model, model=="pborder_weight","border + weights"),
         model = replace(model, model=="ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="pweights","weights"),
         model = replace(model, model=="ptrend","trends"),
         model = replace(model, model=="pborder","border"),
         model = replace(model, model=="pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/hom_pois.png",res=100)
plot_specs(results_phom,
           choices = c("model","controls"))
dev.off()

#-----------
### robberies

# level
results_rob <- run_specs(df = mml, 
                         y = c("robbery_reported_rate"),
                         x = c("mmlXborder"),
                         model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                         controls = c(covars) 
)
results_rob <- results_rob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/rob_level.png",res=100)
plot_specs(results_rob,
           choices = c("model","controls"))
dev.off()

# log(y+1)
results_lrob <- run_specs(df = mml, 
                          y = c("lrob"),
                          x = c("mmlXborder"),
                          model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                          controls = c(covars) 
)
results_lrob <- results_lrob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/rob_log.png",res=100)
plot_specs(results_lrob,
           choices = c("model","controls"))
dev.off()

# inverse hyperbolic sine
results_ihsrob <- run_specs(df = mml, 
                            y = c("ihs_rob"),
                            x = c("mmlXborder"),
                            model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                            controls = c(covars) 
)
results_ihsrob <- results_ihsrob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/rob_ihs.png",res=100)
plot_specs(results_ihsrob,
           choices = c("model","controls"))
dev.off()

# poisson
results_prob <- run_specs(df = mml, 
                          y = c("robbery_reported"),
                          x = c("mmlXborder"),
                          model = c("pbaseline","ptrend","pweights","pborder","ptrend_weight","ptrend_border","pborder_weight","ptrend_border_weight"),
                          controls = c(covars) 
)
results_prob <- results_prob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="ptrend_weight","trends + weights"),
         model = replace(model, model=="ptrend_border","trends + border"),
         model = replace(model, model=="pborder_weight","border + weights"),
         model = replace(model, model=="ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="pweights","weights"),
         model = replace(model, model=="ptrend","trends"),
         model = replace(model, model=="pborder","border"),
         model = replace(model, model=="pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/rob_pois.png",res=100)
plot_specs(results_prob,
           choices = c("model","controls"))
dev.off()

#----------
# assaults

# level
results_ass <- run_specs(df = mml, 
                         y = c("assault_reported_rate"),
                         x = c("mmlXborder"),
                         model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                         controls = c(covars) 
)
results_ass <- results_ass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/ass_level.png",res=100)
plot_specs(results_ass,
           choices = c("model","controls"))
dev.off()

# log(y+1)
results_lass <- run_specs(df = mml, 
                          y = c("lass"),
                          x = c("mmlXborder"),
                          model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                          controls = c(covars) 
)
results_lass <- results_lass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/ass_log.png",res=100)
plot_specs(results_lass,
           choices = c("model","controls"))
dev.off()

# inverse hyperbolic sine
results_ihsass <- run_specs(df = mml, 
                            y = c("ihs_ass"),
                            x = c("mmlXborder"),
                            model = c("baseline","trend","weights","border","trend_weight","trend_border","border_weight","trend_border_weight"),
                            controls = c(covars) 
)
results_ihsass <- results_ihsass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="trend_weight","trends + weights"),
         model = replace(model, model=="trend_border","trends + border"),
         model = replace(model, model=="border_weight","border + weights"),
         model = replace(model, model=="trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/ass_ihs.png",res=100)
plot_specs(results_ihsass,
           choices = c("model","controls"))
dev.off()

# poisson
results_pass <- run_specs(df = mml, 
                          y = c("assault_reported"),
                          x = c("mmlXborder"),
                          model = c("pbaseline","ptrend","pweights","pborder","ptrend_weight","ptrend_border","pborder_weight","ptrend_border_weight"),
                          controls = c(covars) 
)
results_pass <- results_pass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == covars, "all covariates" ),
         model = replace(model, model=="ptrend_weight","trends + weights"),
         model = replace(model, model=="ptrend_border","trends + border"),
         model = replace(model, model=="pborder_weight","border + weights"),
         model = replace(model, model=="ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="pweights","weights"),
         model = replace(model, model=="ptrend","trends"),
         model = replace(model, model=="pborder","border"),
         model = replace(model, model=="pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/ass_pois.png",res=100)
plot_specs(results_pass,
           choices = c("model","controls"))
dev.off()

# overall, IHS is very similar to log; just show one

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------

# event study graph
es_mml <- function(data,depvar,pois,ddd,window){
  data <- data %>%
  # data <- mml %>%
    # make relative year indicator: for triple diff, need separate indicators for border states
    mutate(rel_year = ifelse(treat_year==9999,-1, year-treat_year),
           border_rel_year = ifelse(border==1,ifelse(treat_year==9999,-1, year-treat_year),-1)) 
    # assign -1 to nevertreated group; omitted time category
  
  # get the minimum relative year - we need this to reindex
  min_year <- min(data$rel_year) # -18
  min_byear <- min(data$border_rel_year) # -16
  
  # reindex the relative years
  data <- data %>%
    mutate(rel_year = rel_year - min_year) %>%
    dummy_cols(select_columns = "rel_year")
  data <- data %>%
    mutate(border_rel_year = border_rel_year - min_byear) %>%
    dummy_cols(select_columns = "border_rel_year")
    
  table(data$border_rel_year) # omitted year=15
  table(data$rel_year)  # omitted year = 17
  
  # make regression formula 
    # include all time indicators, no binning or dropping

  indics1 <- paste("rel_year", (0:(max(data$rel_year)))[-(-min_year)], sep = "_", collapse = " + ")
  indics2 <- paste("border_rel_year", (0:(max(data$border_rel_year)))[-(-min_byear)], sep = "_", collapse = " + ")
      # 15,16 years after treatment is California in 2011,2012: no variation in Border, hence collinear with DD indicator

  indics <- paste0(indics1," + ",indics2)
  
  # formula <- as.formula(paste("violent_rate~", indics,"+",covars,"| fips + factor(state)[year] + factor(year)[border]"))
  formula <- as.formula(paste(depvar,"~", indics,"+",covars,"| fips + factor(state)[year] + factor(year)[border]"))
  
  # run mod
  if (pois==0) {
    mod <- feols(formula, cluster=~state,weights=~population_weight,data = data)
  } else if (pois==1) {
    mod <- fepois(formula, cluster=~state,weights=~population_weight,data = data)  
  }
  
  if (window==1){
    win <- c(-15:-2, 0:13)
    breaks <- seq(-15,13,2)
  } else if (window==0) {
    win <- c(-5:-2,0:5)
    # win <- c(-5:5)
    breaks <- -5:5
  }
  
  if (ddd==1) {
    keepvars <- paste("border_rel_year", win - min_byear, sep = "_")
  } else if (ddd==0) {
    keepvars <- paste("rel_year", win - min_year, sep = "_")  
    # might need a different 'win' for DD (ddd=0)
  }

  # grab the obs we need
  es <- broom::tidy(mod, conf.int = TRUE) %>% 
    filter(term %in% keepvars) %>%
    mutate(t = win) %>%
    select(t, estimate, conf.low, conf.high)
  
  es %>%
    bind_rows(tibble(t = -1, estimate = 0, conf.low = 0, conf.high = 0)) %>% 
    ggplot(aes(x = t, y = estimate)) + 
    geom_linerange(aes(ymin = conf.low, ymax = conf.high), color = 'darkgrey', size = 2) + 
    geom_point(color = 'blue', size = 3) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = -1, linetype = "dashed",alpha=0.5) + 
    scale_x_continuous(breaks = breaks) +
    labs(x = "Relative Time", y = "Estimate") + 
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
}

es_mml(mml,"violent_rate",pois=0,ddd=1,window=1)


png(file="/home/michael/Dropbox/blog/mml/output/es_violent.png",res=85)
es_mml(mml,"violent_rate",0,1,0)
dev.off()

es_mml(mml,"murder_reported_rate",0,1,1)

png(file="/home/michael/Dropbox/blog/mml/output/es_hom.png",res=85)
es_mml(mml,"murder_reported_rate",0,1,0)
dev.off()

es_mml(mml,"robbery_reported_rate",0,1,1)

png(file="/home/michael/Dropbox/blog/mml/output/es_rob.png",res=85)
es_mml(mml,"robbery_reported_rate",0,1,0)
dev.off()

es_mml(mml,"assault_reported_rate",0,1,1)

png(file="/home/michael/Dropbox/blog/mml/output/es_ass.png",res=85)
es_mml(mml,"assault_reported_rate",0,1,0)
dev.off()

# poisson
es_mml(mml,"murder_reported",1,1,1)
es_mml(mml,"robbery_reported",1,1,1) 
es_mml(mml,"assault_reported",1,1,1)


#----------------------------
# run GKZ binned event study
es_mml_bin <- function(data,depvar,pois,ddd,unw){
  data <- data %>%
  # data <- mml %>%
  # data <- mml_full %>%
    mutate(rel_year = ifelse(treat_year==9999,-1, year-treat_year),
           border_rel_year = ifelse(border==1,ifelse(treat_year==9999,-1, year-treat_year),-1))    
  min_year <- min(data$rel_year) # -18
  min_byear <- min(data$border_rel_year) # -16   
  data <- data %>%
    mutate(post = ifelse(rel_year >= 5,1,0)) %>% 
    mutate(rel_year = rel_year - min_year) %>%
    dummy_cols(select_columns = "rel_year")
  data <- data %>%
    mutate(bpost = ifelse(border_rel_year >= 5,1,0)) %>% 
    mutate(border_rel_year = border_rel_year - min_byear) %>%
    dummy_cols(select_columns = "border_rel_year")
  # treatment year is 16 (borderrelyear) and 18 (relyear)
    # want to include: 
      # border_rel_year: 11-15,16,17-20,bpost
      # rel_year: 13-17,18,19-22,post
  indics2 <- paste0(paste("border_rel_year", c(((-1-min_byear)-4):((-1-min_byear)+5)), sep = "_", collapse = " + ")," + bpost") # extended ddd
  indics1 <- paste0(paste("rel_year", c(((-1-min_year)-4):((-1-min_year)+5)), sep = "_", collapse = " + ")," + post")
  # indics2 <- paste0(paste("border_rel_year", c(((-1-min_byear)-1):((-1-min_byear)+5)), sep = "_", collapse = " + ")," + bpost") # extended ddd
  # indics1 <- paste0(paste("rel_year", c(((-1-min_year)-1):((-1-min_year)+5)), sep = "_", collapse = " + ")," + post")
      # to replicate GKZ; remember to add b4+b7 to get their coefficient
      # use mml_full, don't drop any obs
  indics <- paste0(indics1," + ",indics2)
  
  win <- c(-5:4)
  breaks <- -5:5
  # win <- c(-2:4)
  # breaks <- -2:5
    # to replicate GKZ
  
  if (ddd==1) {
    keepvars <- c(paste("border_rel_year", win -min_byear, sep = "_"), "bpost")
  } else if (ddd==0) {
    keepvars <- c(paste("rel_year", win -min_year, sep = "_") , "post")
  }
  
  # formula <- as.formula(paste("violent_rate~", indics,"+",covars,"| fips + factor(state)[year] + factor(year)[border]"))
  formula <- as.formula(paste(depvar,"~", indics,"+",covars,"| fips + factor(state)[year] + factor(year)[border]"))
  
  # run mod
  if (pois==0) {
    if (unw==0){
    mod <- feols(formula, cluster=~state,weights=~population_weight,data = data)
    } else if (unw==1) {
      mod <- feols(formula, cluster=~state,data = data)
    }
  } else if (pois==1) {
    if (unw==0) {
    mod <- fepois(formula, cluster=~state,weights=~population_weight,data = data)
    } else if (unw==1) {
    mod <- fepois(formula, cluster=~state,data = data)
    }
  }

  # grab the obs we need
  es <- broom::tidy(mod, conf.int = TRUE) %>% 
    filter(term %in% keepvars) %>%
    mutate(t = c(-5:5)) %>%
    select(t, estimate, conf.low, conf.high)
  
  es %>%
    ggplot(aes(x = t, y = estimate)) + 
    geom_linerange(aes(ymin = conf.low, ymax = conf.high), color = 'darkgrey', size = 2) + 
    geom_point(color = 'blue', size = 3) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_x_continuous(breaks = breaks) +
    # geom_vline(xintercept = 0, linetype = "dashed",alpha=0.5) + 
    labs(x = "Relative Time", y = "Estimate") + 
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
}

# ddd
png(file="/home/michael/Dropbox/blog/mml/output/es_violent_bin.png",res=85)
es_mml_bin(mml,"violent_rate",0,1,0)
dev.off()

png(file="/home/michael/Dropbox/blog/mml/output/es_hom_bin.png",res=85)
es_mml_bin(mml,"murder_reported_rate",0,1,0)
dev.off()

png(file="/home/michael/Dropbox/blog/mml/output/es_rob_bin.png",res=85)
es_mml_bin(mml,"robbery_reported_rate",0,1,0)
dev.off()
# unweighted robbery results
png(file="/home/michael/Dropbox/blog/mml/output/es_rob_bin_unw.png",res=85)
es_mml_bin(mml,"robbery_reported_rate",0,1,1) 
dev.off()

png(file="/home/michael/Dropbox/blog/mml/output/es_ass_bin.png",res=85)
es_mml_bin(mml,"assault_reported_rate",0,1,0)
dev.off()

# poisson
png(file="/home/michael/Dropbox/blog/mml/output/es_phom_bin.png",res=85)
es_mml_bin(mml,"murder_reported",1,1,0)
dev.off()
png(file="/home/michael/Dropbox/blog/mml/output/es_prob_bin.png",res=85)
es_mml_bin(mml,"robbery_reported",1,1,1) # doesn't converge; need to drop blackrate>1? YES.
dev.off()
png(file="/home/michael/Dropbox/blog/mml/output/es_pass_bin.png",res=85)
es_mml_bin(mml,"assault_reported",1,1,0)
dev.off()

# log-level
png(file="/home/michael/Dropbox/blog/mml/output/es_lhom_bin.png",res=85)
es_mml_bin(mml,"lhom",1,1,0)
dev.off()
png(file="/home/michael/Dropbox/blog/mml/output/es_lrob_bin.png",res=85)
es_mml_bin(mml,"lrob",1,1,1) 
dev.off()
png(file="/home/michael/Dropbox/blog/mml/output/es_lass_bin.png",res=85)
es_mml_bin(mml,"lass",1,1,0)
dev.off()

# dd
es_mml_bin(mml,"murder_reported_rate",0,0)
es_mml_bin(mml,"robbery_reported_rate",0,0)
es_mml_bin(mml,"assault_reported_rate",0,0)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### Synthetic control

# aggregate to state level

mml3 <- mml %>%
  group_by(state,year) %>%
  mutate(pop_state=sum(pop)) %>%
  ungroup %>%
  mutate(pop_w = pop/pop_state)

# state_mml <- mml %>%
state_mml <- mml3 %>%
  group_by(state,state_name,year) %>%
  summarize(
    border = mean(border),
    mml = mean(mml),
    treat_year = mean(treat_year),
    pop_state = sum(pop),
    hom_state = sum(murder_reported),
    rob_state = sum(robbery_reported),
    ass_state = sum(assault_reported),
    male_state = sum(male),
    black_state = sum(black),
    hisp_state = sum(hisp),
    age1019_state = sum(age1019),
    age2024_state = sum(age2024),
    dec_state = mean(dec),
    medinc_state = weighted.mean(median_income,pop_w),
    unempl_state = weighted.mean(unemploymentrate,pop_w),
    pov_state = weighted.mean(poverty_all_percent,pop_w)
  )

state_mml <- state_mml %>%
  mutate(
    hom_rate = (hom_state/pop_state)*100000,
    rob_rate = (rob_state/pop_state)*100000,
    ass_rate = (ass_state/pop_state)*100000,
    male_rate = (male_state/pop_state),
    black_rate =(black_state/pop_state),
    hisp_rate = (hisp_state/pop_state),
    age1019_rate = (age1019_state/pop_state),
    age2024_rate = (age2024_state/pop_state),
    mmlXborder = mml*border,
    lnmedian_income = log(medinc_state),
    lnpop = log(pop_state)
  )

#---------------------------------
# california: 6, treated in 1996
cali_mml <- state_mml %>%
  # filter(state==6 | treat_year==9999 | year < treat_year)
  # can't use treated states in pre-treatment years; need balanced panel
  filter(state==6 | treat_year==9999)

# hom
dataprep.cali_hom <- 
  dataprep(
    foo = as.data.frame(cali_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "hom_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 6,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:1995,
    time.optimize.ssr = 1994:1995,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

cali.synth_hom <- synth(dataprep.cali_hom)
# 99.98% of weight is on texas
png(file="/home/michael/Dropbox/blog/mml/output/sc_cali_hom.png",res=100)
path.plot(dataprep.res = dataprep.cali_hom,synth.res = cali.synth_hom,Ylab="Homicide rate",Xlab='',Main='Synthetic control: California, homicides',tr.intake=1996)
dev.off()
gaps.plot(dataprep.res = dataprep.cali_hom,synth.res = cali.synth_hom)

# rob
dataprep.cali_rob <- 
  dataprep(
    foo = as.data.frame(cali_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "rob_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 6,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:1995,
    time.optimize.ssr = 1994:1995,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

cali.synth_rob <- synth(dataprep.cali_rob)
# 68% on 36 New York, 28% on 27 Minnesota

png(file="/home/michael/Dropbox/blog/mml/output/sc_cali_rob.png",res=100)
path.plot(dataprep.res = dataprep.cali_rob,synth.res = cali.synth_rob,Ylab="Robbery rate",Xlab='',Main='Synthetic control: California, robberies',tr.intake=1996)
dev.off()

gaps.plot(dataprep.res = dataprep.cali_rob,synth.res = cali.synth_rob)

# ass
dataprep.cali_ass <- 
  dataprep(
    foo = as.data.frame(cali_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "ass_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 6,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:1995,
    time.optimize.ssr = 1994:1995,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

cali.synth_ass <- synth(dataprep.cali_ass)
# 67% on texas, 32% on 12 Florida

png(file="/home/michael/Dropbox/blog/mml/output/sc_cali_ass.png",res=100)
path.plot(dataprep.res = dataprep.cali_ass,synth.res = cali.synth_ass,Ylab="Assault rate",Xlab='',Main='Synthetic control: California, assaults',tr.intake=1996)
dev.off() 

gaps.plot(dataprep.res = dataprep.cali_ass,synth.res = cali.synth_ass)


#------------------------------------
# arizona: 4, treated in 2010
ariz_mml <- state_mml %>%
  filter(state_name=="Arizona"|treat_year==9999)

# hom
dataprep.ariz_hom <- 
  dataprep(
    foo = as.data.frame(ariz_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "hom_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 4,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:2009,
    time.optimize.ssr = 1994:2009,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

ariz.synth_hom <- synth(dataprep.ariz_hom)
# homicide: 63% on Texas (48), 25% on Florida, 6% on 51 Virginia

png(file="/home/michael/Dropbox/blog/mml/output/sc_ariz_hom.png",res=100)
path.plot(dataprep.res = dataprep.ariz_hom,synth.res = ariz.synth_hom,Ylab="Homicide rate",Xlab='',Main='Synthetic control: Arizona, homicides',tr.intake=2010)
dev.off()

gaps.plot(dataprep.res = dataprep.ariz_hom,synth.res = ariz.synth_hom)

# rob
dataprep.ariz_rob <- 
  dataprep(
    foo = as.data.frame(ariz_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "rob_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 4,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:2009,
    time.optimize.ssr = 1994:2009,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

ariz.synth_rob <- synth(dataprep.ariz_rob)
# robbery: 61% on texas, 24% on 12 Florida, 15% on 56 Wyoming

png(file="/home/michael/Dropbox/blog/mml/output/sc_ariz_rob.png",res=100)
path.plot(dataprep.res = dataprep.ariz_rob,synth.res = ariz.synth_rob,Ylab="Robbery rate",Xlab='',Main='Synthetic control: Arizona, robberies',tr.intake=2010)
dev.off() 

gaps.plot(dataprep.res = dataprep.ariz_rob,synth.res = ariz.synth_rob)

# ass
dataprep.ariz_ass <- 
  dataprep(
    foo = as.data.frame(ariz_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "ass_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 4,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:2009,
    time.optimize.ssr = 1994:2009,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

ariz.synth_ass <- synth(dataprep.ariz_ass)
# 67% on texas, 14% on wyoming 56, 13% on florida 12

png(file="/home/michael/Dropbox/blog/mml/output/sc_ariz_ass.png",res=100)
path.plot(dataprep.res = dataprep.ariz_ass,synth.res = ariz.synth_ass,Ylab="Assault rate",Xlab='',Main='Synthetic control: Arizona, assaults',tr.intake=2010)
dev.off()

gaps.plot(dataprep.res = dataprep.ariz_ass,synth.res = ariz.synth_ass)

#----------------------------------------------
# new mexico: 35, treated in 2007
nmex_mml <- state_mml %>%
  filter(state_name=="New Mexico"|treat_year==9999)

# hom
dataprep.nmex_hom <- 
  dataprep(
    foo = as.data.frame(nmex_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "hom_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 35,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:2006,
    time.optimize.ssr = 1994:2006,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

nmex.synth_hom <- synth(dataprep.nmex_hom)
# homicide: 62% on West Virginia (54), 6% on Wyoming (56), 28% on Louisiana (22), 3% on 28 Mississippi

png(file="/home/michael/Dropbox/blog/mml/output/sc_nmex_hom.png",res=100)
path.plot(dataprep.res = dataprep.nmex_hom,synth.res = nmex.synth_hom,Ylab="Homicide rate",Xlab='',Main='Synthetic control: New Mexico, homicides',tr.intake=2007)
dev.off() 

gaps.plot(dataprep.res = dataprep.nmex_hom,synth.res = nmex.synth_hom)

# rob
dataprep.nmex_rob <- 
  dataprep(
    foo = as.data.frame(nmex_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "rob_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 35,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:2006,
    time.optimize.ssr = 1994:2006,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

nmex.synth_rob <- synth(dataprep.nmex_rob)
# robbery: 51% on 28 Mississippi, 21% on 22 Louisiana, 18% on texas, 7% on Wyoming 56

png(file="/home/michael/Dropbox/blog/mml/output/sc_nmex_rob.png",res=100)
path.plot(dataprep.res = dataprep.nmex_rob,synth.res = nmex.synth_rob,Ylab="Robbery rate",Xlab='',Main='Synthetic control: New Mexico, robberies',tr.intake=2007)
dev.off()

gaps.plot(dataprep.res = dataprep.nmex_rob,synth.res = nmex.synth_rob)


# ass
dataprep.nmex_ass <- 
  dataprep(
    foo = as.data.frame(nmex_mml),
    predictors = c("lnmedian_income","unempl_state","male_rate","black_rate","hisp_rate","age2024_rate","age1019_rate","lnpop","pov_state","dec_state"),
    predictors.op = "mean",
    dependent = "ass_rate",
    unit.variable = "state",
    time.variable = "year",
    treatment.identifier = 35,
    controls.identifier = c(1,5,12,13,16:22,24,25,27:29,31,33,36:40,42,45:49,51,54:56),
    time.predictors.prior = 1994:2006,
    time.optimize.ssr = 1994:2006,
    unit.names.variable = c("state_name"),
    time.plot = 1994:2012
  )

nmex.synth_ass <- synth(dataprep.nmex_ass)
# 55% on Texas, 34% on Florida (12), 11% on Louisiana (22)

png(file="/home/michael/Dropbox/blog/mml/output/sc_nmex_ass.png",res=100)
path.plot(dataprep.res = dataprep.nmex_ass,synth.res = nmex.synth_ass,Ylab="Assault rate",Xlab='',Main='Synthetic control: New Mexico, assaults',tr.intake=2007)
dev.off()

gaps.plot(dataprep.res = dataprep.nmex_ass,synth.res = nmex.synth_ass)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# spec curve (state-level)

# models
s_baseline <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml", "|state+year"))
  feols(formula,cluster=~state,data)
}
s_trend <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|state + factor(state)[year]"))
  feols(formula,cluster=~state,data)
}
s_weights <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state+year"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}
s_border <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state + factor(year)[border]"))
  feols(formula,cluster=~state,data)
}
s_trend_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|state + factor(state)[year]"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}
s_trend_border <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|state + factor(state)[year] + factor(year)[border]"))
  feols(formula,cluster=~state,data)
}
s_border_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state + factor(year)[border]"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}
s_trend_border_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state + factor(state)[year] + factor(year)[border]"))
  feols(formula,cluster=~state,weights=~population_weight,data)
}

s_pbaseline <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml", "|state+year"))
  fepois(formula,cluster=~state,data)
}
s_ptrend <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|state + factor(state)[year]"))
  fepois(formula,cluster=~state,data)
}
s_pweights <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state+year"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}
s_pborder <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state + factor(year)[border]"))
  fepois(formula,cluster=~state,data)
}
s_ptrend_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|state + factor(state)[year]"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}
s_ptrend_border <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+mml|state + factor(state)[year] + factor(year)[border]"))
  fepois(formula,cluster=~state,data)
}
s_pborder_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state + factor(year)[border]"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}
s_ptrend_border_weight <- function(formula,data) {
  formula <-as.formula(paste0(formula,"+ mml | state + factor(state)[year] + factor(year)[border]"))
  fepois(formula,cluster=~state,weights=~population_weight,data)
}

#-----------------------
### violent rate

# level
s_results_viol <- run_specs(df = state_mml, 
                            y = c("violent_rate"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars) 
)
s_results_viol <- s_results_viol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_baseline","baseline"),
         model = replace(model, model=="s_border","border"),
         model = replace(model, model=="s_trend","trend"),
         model = replace(model, model=="s_weights","weights"),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_viol_level.png",res=100)
plot_specs(s_results_viol,
           choices = c("model","controls"))
dev.off()

# log(y+1)
s_results_lviol <- run_specs(df = state_mml, 
                             y = c("lviol"),
                             x = c("mmlXborder"),
                             model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                             controls = c(scovars)
)
s_results_lviol <- s_results_lviol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_viol_log.png",res=100)
plot_specs(s_results_lviol,
           choices = c("model","controls"))
dev.off()

# poisson
s_results_pviol <- run_specs(df = state_mml, 
                             y = c("violent_level"),
                             x = c("mmlXborder"),
                             model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                             controls = c(scovars)
)
s_results_pviol <- s_results_pviol %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_ptrend_weight","trends + weights"),
         model = replace(model, model=="s_ptrend_border","trends + border"),
         model = replace(model, model=="s_pborder_weight","border + weights"),
         model = replace(model, model=="s_ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="s_pweights","weights"),
         model = replace(model, model=="s_ptrend","trends"),
         model = replace(model, model=="s_pborder","border"),
         model = replace(model, model=="s_pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_viol_pois.png",res=100)
plot_specs(s_results_pviol,
           choices = c("model","controls"))
dev.off()

#---------------
### homicides

# level
s_results_hom <- run_specs(df = state_mml, 
                           y = c("hom_rate"),
                           x = c("mmlXborder"),
                           model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                           controls = c(scovars) 
)
s_results_hom <- s_results_hom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_baseline","baseline"),
         model = replace(model, model=="s_border","border"),
         model = replace(model, model=="s_trend","trend"),
         model = replace(model, model=="s_weights","weights"),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_hom_level.png",res=100)
plot_specs(s_results_hom,
           choices = c("model","controls"))
dev.off()

# log(y+1)
s_results_lhom <- run_specs(df = state_mml, 
                            y = c("lhom"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars)
)
s_results_lhom <- s_results_lhom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_hom_log.png",res=100)
plot_specs(s_results_lhom,
           choices = c("model","controls"))
dev.off()

# poisson
s_results_phom <- run_specs(df = state_mml, 
                            y = c("hom_state"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars)
)
s_results_phom <- s_results_phom %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_ptrend_weight","trends + weights"),
         model = replace(model, model=="s_ptrend_border","trends + border"),
         model = replace(model, model=="s_pborder_weight","border + weights"),
         model = replace(model, model=="s_ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="s_pweights","weights"),
         model = replace(model, model=="s_ptrend","trends"),
         model = replace(model, model=="s_pborder","border"),
         model = replace(model, model=="s_pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_hom_pois.png",res=100)
plot_specs(s_results_phom,
           choices = c("model","controls"))
dev.off()

#-----------------------
### robberies

# level
s_results_rob <- run_specs(df = state_mml, 
                           y = c("rob_rate"),
                           x = c("mmlXborder"),
                           model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                           controls = c(scovars) 
)
s_results_rob <- s_results_rob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_baseline","baseline"),
         model = replace(model, model=="s_border","border"),
         model = replace(model, model=="s_trend","trend"),
         model = replace(model, model=="s_weights","weights"),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_rob_level.png",res=100)
plot_specs(s_results_rob,
           choices = c("model","controls"))
dev.off()

# log(y+1)
s_results_lrob <- run_specs(df = state_mml, 
                            y = c("lrob"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars)
)
s_results_lrob <- s_results_lrob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_rob_log.png",res=100)
plot_specs(s_results_lrob,
           choices = c("model","controls"))
dev.off()

# poisson
s_results_prob <- run_specs(df = state_mml, 
                            y = c("rob_state"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars)
)
s_results_prob <- s_results_prob %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_ptrend_weight","trends + weights"),
         model = replace(model, model=="s_ptrend_border","trends + border"),
         model = replace(model, model=="s_pborder_weight","border + weights"),
         model = replace(model, model=="s_ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="s_pweights","weights"),
         model = replace(model, model=="s_ptrend","trends"),
         model = replace(model, model=="s_pborder","border"),
         model = replace(model, model=="s_pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_rob_pois.png",res=100)
plot_specs(s_results_prob,
           choices = c("model","controls"))
dev.off()

#-----------------------
### assaults

# level
s_results_ass <- run_specs(df = state_mml, 
                           y = c("ass_rate"),
                           x = c("mmlXborder"),
                           model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                           controls = c(scovars) 
)
s_results_ass <- s_results_ass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_baseline","baseline"),
         model = replace(model, model=="s_border","border"),
         model = replace(model, model=="s_trend","trend"),
         model = replace(model, model=="s_weights","weights"),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_ass_level.png",res=100)
plot_specs(s_results_ass,
           choices = c("model","controls"))
dev.off()

# log(y+1)
s_results_lass <- run_specs(df = state_mml, 
                            y = c("lass"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars)
)
s_results_lass <- s_results_lass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_trend_weight","trends + weights"),
         model = replace(model, model=="s_trend_border","trends + border"),
         model = replace(model, model=="s_border_weight","border + weights"),
         model = replace(model, model=="s_trend_border_weight","trends + border + weights")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_ass_log.png",res=100)
plot_specs(s_results_lass,
           choices = c("model","controls"))
dev.off()

# poisson
s_results_pass <- run_specs(df = state_mml, 
                            y = c("ass_state"),
                            x = c("mmlXborder"),
                            model = c("s_baseline","s_trend","s_weights","s_border","s_trend_weight","s_trend_border","s_border_weight","s_trend_border_weight"),
                            controls = c(scovars)
)
s_results_pass <- s_results_pass %>%
  distinct(across(x:subsets)) %>%
  mutate(controls = replace(controls, controls == scovars, "all covariates" ),
         model = replace(model, model=="s_ptrend_weight","trends + weights"),
         model = replace(model, model=="s_ptrend_border","trends + border"),
         model = replace(model, model=="s_pborder_weight","border + weights"),
         model = replace(model, model=="s_ptrend_border_weight","trends + border + weights"),
         model = replace(model, model=="s_pweights","weights"),
         model = replace(model, model=="s_ptrend","trends"),
         model = replace(model, model=="s_pborder","border"),
         model = replace(model, model=="s_pbaseline","baseline")
  )
png(file="/home/michael/Dropbox/blog/mml/output/s_ass_pois.png",res=100)
plot_specs(s_results_prob,
           choices = c("model","controls"))
dev.off()
