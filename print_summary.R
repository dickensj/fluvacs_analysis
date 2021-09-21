library(dplyr)
library(rstan)
#library(stringr)




full_res <- readRDS('results/full_age_wane_res.rds')
ols_res <- readRDS('results/ols_age_wane_res.rds')

full_out <- summary(full_res, pars=c('beta_cov', 'beta_treat', 'beta_flu', 'beta_batch', 'k', 'beta_flu', 'tau', 'k0', 'sigma'), prob=.5)
ols_out <- summary(ols_res, pars=c('beta_cov', 'beta_treat', 'beta_flu', 'beta_batch', 'k', 'beta_flu', 'tau', 'k0', 'sigma'), prob=.5)
fs <- full_out$summary[,c(1, 3)]
os <- ols_out$summary[,c(1, 3)]
out <- data.frame(round(cbind(fs, os)[,c(1,3,2,4)], 2))
names(out) <- c('full_est', 'ols_est', 'full_se', 'ols_est')
print(out)

