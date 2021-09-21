library(stringr)
library(rstan)
fi <- list.files('results', full=TRUE)
full_fi <- fi[str_detect(fi, 'full_age_wane')]
ols_fi <- fi[str_detect(fi, 'ols_age_wane')]



if (any(str_detect('full_age_wane_res.rds', fi))) {
	full_res <- readRDS('results/full_age_wane_res.rds')
	print(full_res, 'beta_treat')
} else {
	full_res <- rstan::read_stan_csv(full_fi)
	print(full_res)
	saveRDS(full_res, 'results/full_age_wane_res.rds')
}

if (any(str_detect('ols_age_wane_res.rds', fi))) {
	ols_res <- readRDS('results/ols_age_wane_res.rds')
} else {
	ols_res <- rstan::read_stan_csv(ols_fi)
	print(ols_res)
	saveRDS(ols_res, 'results/ols_age_wane_res.rds')
}
