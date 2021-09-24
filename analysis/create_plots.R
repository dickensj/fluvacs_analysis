library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(cmdstanr)
library(posterior)
rstan_options(auto_write=T)

setwd('C:/Users/15636/Documents/GitHub/fluvacs_analysis/')
# setwd('~/Dropbox (University of Michigan)/fluvacs_analysis/')

fr <- readRDS('results/full_age_wane_res.rds')
or <- readRDS('results/ols_age_wane_res.rds')

print(fr, 'alpha0')
print(fr, 'k')

# what are the boosting effects by age

bt <- rstan::extract(fr, 'beta_treat')[[1]]
treat_design <- data.frame(vax1 = c(1, 1, 0, 0),
                           vax2 = c(0, 0, 1, 1),
                           vax1age = c(0, 3, 0, 0),
                           vax2age = c(0, 0, 0, 3))
boosts <- data.matrix(treat_design) %*% t(bt)
mu <- apply(boosts, 1, mean)
treat_design$mu <- mu
treat_design$age <- c(18, 48, 18, 48)
treat_design$vaxcode <- c('vax1', 'vax1', 'vax2', 'vax2')

ggplot(treat_design, aes(age, mu, color=vaxcode)) +
    geom_line() +
    labs(x = 'Subject age',
         y = 'Post-vaccination antibody boost') +
    scale_color_discrete('Study arm arm', labels = c('Flumist', 'Fluvax')) +
    theme_classic()

# waning rate comparison
k <- rstan::extract(fr, 'k')[[1]]
k0 <- rstan::extract(fr, 'k0')[[1]]
vc <- data.frame('vaxcode' = 0:2)
age <- data.frame(age = seq(0, 3, by = .1))

design <- merge(vc, age) %>%
    mutate(vax1 = vaxcode == 1,
           vax2 = vaxcode == 2,
           vax1age = vax1 * age,
           vax2age = vax2 * age) %>%
    as_tibble()

wr <- select(design, age, vax1, vax2, vax1age, vax2age) %>%
    data.matrix() %*% t(k)
design$wane_rate <- apply(exp(wr), 1, function(x) mean(k0 * x))
design$wane_rate_lb <- apply(exp(wr), 1, function(x) quantile(k0 * x, .05))
design$wane_rate_ub <- apply(exp(wr), 1, function(x) quantile(k0 * x, .95))


ggplot(design, aes(age, wane_rate, color=factor(vaxcode))) +
    geom_line() +
    geom_ribbon(aes(x = age ,
                    ymin = wane_rate_lb,
                    ymax = wane_rate_ub,
                    fill = factor(vaxcode)), alpha = .25, linetype = 'blank') +
    labs(x = 'Subject age',
         y = 'Antibody waning rate (1 / yr)') +
    scale_fill_discrete('Study arm', labels = c('Control', 'Flumist', 'Fluvax')) +
    scale_color_discrete(guide = 'none') +
    theme_classic()


ggplot(design, aes(18 + 10 * age, log(2) / wane_rate, color=factor(vaxcode))) +
    geom_line() +
    # geom_ribbon(aes(ymin = log(2) /wane_rate_lb,
    #                 ymax = log(2) /wane_rate_ub,
    #                 fill = factor(vaxcode)), alpha = .25, linetype = 'blank') +
    labs(x = 'Subject age',
         y = 'Antibody half life (years)') +
    scale_color_discrete('Study arm', labels = c('Control', 'Flumist', 'Fluvax')) +
    # scale_color_discrete(guide = 'none') +
    scale_y_continuous(limits = c(0, NA), breaks = 0:6) +
    theme_classic()


# fitted values
dat <- read_rdump('data/fluvacs.data.R')
ymat <- dat$y

est_mu <- rstan::extract(fr, 'ev')[[1]]
mu_mat <- matrix(0, length(batch_eff), 5)
for (j in 1:nrow(mu_mat)) {
    for (i in 1:5) {
        mu_mat[j, i] <- mean(est_mu[,j,i])
    }
}

mu <- as.numeric(mu_mat[dat$w > 0])
resid <- as.numeric(ymat[dat$w > 0] - mu_mat[dat$w > 0])
df <- data.frame(mu, resid, y=as.factor(ymat[dat$w > 0]))

ggplot(df, aes(mu, resid)) +
    geom_point(aes(color=y)) +
    geom_smooth() +
    scale_y_continuous(breaks=-5:10)

ggplot(df, aes(resid, group=y)) +
    geom_density() +
    facet_wrap(. ~ y, ncol=4)


# tables

fr <- readRDS('results/full_age_wane_res.rds')
or <- readRDS('results/ols_age_wane_res.rds')

beta_cens <- rstan::extract(fr, c('beta_treat', 'k0', 'k')) %>%
    do.call(cbind, .)
mu_cens <- round(apply(beta_cens, 2, mean), 2)
se_cens <- round(apply(beta_cens, 2, sd), 2)
bds_cens <- round(t(apply(beta_cens, 2, quantile, p = c(.05, .95))), 2)

npar <- length(mu_cens)
cens_summary <- sapply(1:npar, function(i) {
    stringr::str_c(mu_cens[i], ' (', bds_cens[i,1], ', ', bds_cens[i, 2], ')')
})


beta_mem <- rstan::extract(or, c('beta_treat', 'k0', 'k')) %>%
    do.call(cbind, .)
mu_mem <- round(apply(beta_mem, 2, mean), 2)
se_mem <- round(apply(beta_mem, 2, sd), 2)
bds_mem <- round(t(apply(beta_mem, 2, quantile, p = c(.05, .95))), 2)


npar <- length(mu_mem)
mem_summary <- sapply(1:npar, function(i) {
    stringr::str_c(mu_mem[i], ' (', bds_mem[i,1], ', ', bds_mem[i, 2], ')')
})
xtable::xtable(cbind(cens_summary,mem_summary))
