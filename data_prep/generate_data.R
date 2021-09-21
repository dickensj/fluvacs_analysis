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


set_cmdstan_path("C:/Users/15636/Documents/.cmdstanr/cmdstan-2.26.1")

setwd('C:/Users/15636/Dropbox (University of Michigan)/hai-titer-models/src/validation-sims/fluvacs/')


time_df <- data.frame(sero_timing = factor(c("Prevaccination: Fall 2005",
                                             "Postvaccination / Preseason: Fall 2005",
                                             "Postseason: Spring 2006",
                                             "Preseason: Fall 2006",
                                             "Postseason: Spring 2007")),
                      time_int = 1:5)

time_since_illness <- function(df) {
    sero <- df$sero_date
    id <- df$illness_date
    if (all(is.na(id))) {
        time_since_illness <- rep(0, length(sero))
        illness <- rep(0, length(sero))
    } else {
        time_since_illness <- difftime(sero, id[1], units = 'days') / 365
        illness <- time_since_illness > 0
        time_since_illness <- ifelse(time_since_illness < 0, 0, time_since_illness)
    }
    df$flu <- illness
    df$time_since_illness <- time_since_illness
    df
}

illness_timing <- function(df) {
    sero <- df$sero_date
    id <- df$illness_date
    if (all(is.na(id))) {
        illness_timing <- rep(NA, length(sero))
    } else {
        illness_timing <- as.numeric(difftime(id[1], min(sero), units = 'days') / 365)
    }
    df$illness_timing <- illness_timing
    df
}

demo <- read.csv('data/fluvacs-demographics.csv')
names(demo) <- tolower(names(demo))
titers <- read.csv('data/fluvacsH3NYHAI2004-2007_unadjusted.csv')
names(titers) <- tolower(names(titers))
titers <- rename(titers, 'titer' = 'log2_hai_h3ny')
infections <- read.csv('data/fluvacsH3NYHAI2004-2007_infection.csv')
names(infections) <- tolower(names(infections))

infections <- filter(infections, flua == 1) %>%
    select(-flua, -flub) %>%
    as_tibble() %>%
    mutate(illness_date = ymd(illness_date))

fv <- right_join(demo, titers, by = 'studyid') %>%
    as_tibble() %>%
    full_join(infections, by= 'studyid') %>%
    full_join(time_df, by = 'sero_timing') %>%
    mutate(vaxcode = ifelse(vaxcode == 3, 0, vaxcode)) %>%
    arrange(vaxcode, studyid, time_int, batch) %>%
    mutate(sero_date = mdy(sero_date)) %>%
    group_by(studyid) %>%
    do(time_since_illness(.)) %>%
    do(illness_timing(.))

#####################
# create basic design
#####################

reg_df <- select(fv, -c('studysite')) %>%
    mutate(vax_date = mdy(vax_date),
           vax_status = vaxcode != 0 & time_int != 1,
           time_since_vax = as.numeric(difftime(sero_date, vax_date, units = 'days')) / 365) %>%
    arrange(vaxcode, studyid, time_int, batch)

# create design matrices

design <- reg_df %>%
    ungroup() %>%
    mutate(vax1 = vaxcode == 1,
           vax2 = vaxcode == 2,
           time = time_since_vax,
           age = (age - 18) / 10,
           obs = ifelse(titer== 11, 0, 1),
           batch2 = batch == 'Summer 2007',
           batch3 = batch == 'Summer 2009') %>%
    group_by(studyid) %>%
    mutate(nobs = length(unique(sero_date))) %>%
    filter(nobs == 5) %>%
    arrange(vaxcode, studyid, time_int, batch) %>% ungroup() %>%
    mutate(id = stringr::str_c(studyid, '_', batch))

id_df <- group_by(design, studyid) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(studyid) %>%
    mutate(id_num = 1:n())

design <- full_join(id_df, design, by='studyid')

nom_time_df <- group_by(design, sero_timing) %>%
    filter(row_number() == 1) %>% 
    select(time_int, sero_timing) %>% ungroup() %>%
    arrange(time_int) %>%
    mutate(yname = stringr::str_c('y', time_int))

ids_batch <- unique(design$id)

y <- lapply(ids_batch, function(ii) {
    df <- subset(design, id == ii) %>%
        select(studyid, id_num, vaxcode, batch, sero_timing, titer)
    df <- full_join(df, nom_time_df, by = 'sero_timing')
    df <- pivot_wider(df, id_cols=studyid:batch, names_from=yname, values_from = titer)
    keep <- !apply(df, 1, function(xx) all(is.na(xx)))
    df[keep,]
}) %>% do.call(rbind, .)

w <- y %>% pivot_longer(cols=y1:y5, names_to = 'titer_id', values_to='titer') %>%
    group_by(studyid, titer_id) %>% 
    mutate(w = as.numeric(!is.na(titer)) / sum(!is.na(titer))) %>%
    mutate(titer_id = stringr::str_replace(titer_id, 'y', 'w')) %>%
    select(studyid, id_num, vaxcode, batch, titer_id, w) %>%
    pivot_wider(id_cols=studyid:batch, names_from=titer_id, values_from=w)

obs <- lapply(ids_batch, function(ii) {
    df <- subset(design, id == ii) %>%
        select(studyid, id_num, vaxcode, batch, sero_timing, obs)
    df <- full_join(df, nom_time_df, by = 'sero_timing')
    df <- pivot_wider(df, id_cols=studyid:batch, names_from=yname, values_from = obs)
    keep <- !apply(df, 1, function(xx) all(is.na(xx)))
    df <- df[keep,]
}) %>% do.call(rbind, .)

design_unique <- group_by(design, id_num, sero_timing) %>%
    filter(row_number()==1) %>% ungroup()


ids <- sort(unique(design$id_num))

time_masking <- lapply(ids, function(ii) {
    df <- subset(design_unique, id_num == ii) %>%
        select(time, illness_timing)
    if (all(is.na(df$illness_timing))) {
        tv <- sort(c(df$time, 2/52, max(df$time)))
        mask <- c(1, 3:6)
    } else {
        tv <- sort(c(df$time, 2/52, df$illness_timing[1] - .0001)) # subtracting .0001 takes care of illnesses discovered on serodates
        mask <- which(is.element(tv, df$time))
    }
    list(tv=tv, mask=mask)
})

tvec <- lapply(time_masking, function(lis) lis$tv) %>%
    do.call(rbind, .)

mask <- lapply(time_masking, function(lis) lis$mask) %>%
    do.call(rbind, .)

tflu <- sapply(ids, function(ii) {
    df <- subset(design_unique, id_num == ii) %>%
        select(time, illness_timing)
    if (all(is.na(df$illness_timing))) {
        -1
    } else {
        df$illness_timing[1] - .0001 # to match earlier adjustment
    }
})

Xcov <- design_unique  %>% group_by(id_num) %>% filter(row_number() ==1) %>% ungroup() %>%
    select(id_num, male, age)
Kc <- ncol(Xcov)- 1

Xk <- design_unique  %>% group_by(id_num) %>% filter(row_number() ==1) %>% ungroup() %>%
    mutate(vax1age = vax1*age, vax2age = vax2 * age) %>%
    select(id_num, age, vax1, vax2, vax1age, vax2age)
Kw  = ncol(Xk) - 1

Xtreat <- design  %>% group_by(id_num) %>% filter(row_number() ==1) %>% ungroup() %>%
    mutate(vax1age = vax1*age, vax2age = vax2 * age) %>%
    select(id_num, vax1, vax2, vax1age, vax2age)

Kt <- ncol(Xtreat) - 1

Xflu <- design_unique  %>% group_by(id_num) %>% mutate(flu = max(flu)) %>%  filter(row_number() ==1) %>% ungroup() %>%
    group_by(studyid) %>% ungroup() %>%
    mutate(fluXage = flu*age, fluXvax1 = flu * vax1, fluXvax2 = flu * vax2) %>%
    select(id_num, flu, fluXage, fluXvax1, fluXvax2)
Kf <- ncol(Xflu) -1

all(Xflu$id_num == Xtreat$id_num)
all(Xflu$id_num == Xk$id_num)
all(Xflu$id_num == Xcov$id_num)
all(Xtreat$id_num == Xk$id_num)
all(Xtreat$id_num == Xcov$id_num)
all(Xcov$id_num == Xk$id_num)

mu_id <- y$id_num
N <- nrow(y)
ymat <- select(y, y1:y5) %>% data.matrix()
ymat <- apply(ymat, c(1,2), function(x) ifelse(is.na(x), -1000, x))
wmat <- ungroup(w) %>% select(w1:w5) %>% data.matrix()
obs_mat <- ungroup(obs) %>% select(y1:y5) %>% data.matrix()
obs_mat <- apply(obs_mat, c(1,2), function(x) ifelse(is.na(x), 0, x))

batch_design <- model.matrix(~ batch, data=y)[, 2:3]

Nt <- 5
Nsub <- nrow(Xcov)
N <- nrow(ymat)

dat <- list(N= N,
            Nt = Nt,
            Nsub = Nsub,
            y = ymat,
            w=wmat,
            mask = mask,
            tvec = tvec,
            tflu = tflu,
            Kc = Kc,
            Kt = Kt,
            Kf = Kf,
            Xcov = Xcov %>% select(-id_num) %>% data.matrix(),
            Xtreat = Xtreat %>% select(-id_num) %>% data.matrix(),
            Xflu = Xflu %>% select(-id_num) %>% data.matrix(),
            Xbatch = batch_design,
            Xk=Xk %>% select(-id_num) %>% data.matrix(),
            obs=obs_mat,
            mu_id=mu_id)

N= N
Nt = Nt
Nsub = Nsub
y = ymat
w=wmat
mask = mask
tvec = tvec
tflu = tflu
Kc = Kc
Kt = Kt
Kf = Kf
Xcov = Xcov %>% select(-id_num) %>% data.matrix()
Xtreat = Xtreat %>% select(-id_num) %>% data.matrix()
Xflu = Xflu %>% select(-id_num) %>% data.matrix()
Xbatch = batch_design
Xk=Xk %>% select(-id_num) %>% data.matrix()
obs=obs_mat
mu_id=mu_id
p_re = 3

stan_rdump(list("N", "Nt", "Nsub", "y" , "w", "mask", "tvec", "tflu", "Kc", "Kt", "Kf", "Kw", "Xcov", "Xtreat", "Xflu", "Xbatch", "Xk", "obs", "mu_id", 'p_re'), 'fluvacs.data.R')


sigma=1
alpha0=6
beta_treat = c(1, 2.5, 0, 0)
beta_flu=c(3, 0, -1, -2.5)
tau=c(1,1,1)
Psi = diag(3)
z = matrix(rnorm(dat$Nsub * 3), ncol=3)
k = rep(0, 5)

stan_rdump(list('sigma', 'alpha0', 'beta_treat', 'beta_flu', 'tau', 'Psi', 'z', 'k'), 'inits.data.R')





dat <- read_rdump('fluvacs.data.R')
