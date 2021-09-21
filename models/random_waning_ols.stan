functions {
    real mu_evolve(real mu, real dt, real boost, real k) {
        return mu - dt * k + boost;
    }
    
    real mu_evolve2(real mu, real dt, real boost, real k) {
        return mu / (1 - dt * k * mu);
    }
    
    vector mu_calc(int Nt, real mu0, real vax_boost, real flu_boost, real t_flu, 
                   real k, vector t, int[] mask) {
        real t0 = 0.0;
        real dt = t[1] - t0;
        vector[Nt + 2] mu = rep_vector(0, Nt + 2);
        mu[1] = mu0;
        mu[2] =  mu_evolve(mu[1], dt, vax_boost, k);
        t0 = t[1];
        for (j in 3:(Nt + 2)) {
            dt = t[j] - t0;
            if (t[j] == t_flu) {
                mu[j] = mu_evolve(mu[j-1], dt, flu_boost, k);
            } else {
                mu[j] = mu_evolve(mu[j-1], dt, 0.0, k);
            }
            t0 = t[j];
        }
        return mu[mask]; 
    }
    
    real log_lik_fn(real y, real obs, real mu, real sigma) {
        if (y == 0) 
            return normal_lcdf(y | mu, sigma);  
        else if (obs == 1)
            return log_diff_exp(normal_lcdf(y  | mu, sigma), 
                                normal_lcdf((y - 1) | mu, sigma));
        else
            return normal_lccdf(y | mu, sigma);     
    }
}
data {
    int N;
    int<lower=0> Nt;
    int<lower=0> Nsub;
    vector[Nt] y[N];
    vector[Nt] w[N];
    int mask[Nsub, Nt];
    vector[Nt+2] tvec[Nsub];// 2nd time pt is boost and last tp is throw away for non-flu cases
    vector[Nsub] tflu;
    int<lower=0> Kc;
    int Kt;
    int Kf;
    int Kw;
    matrix[Nsub, Kc] Xcov;
    matrix[Nsub, Kt] Xtreat;
    matrix[Nsub, Kf] Xflu;
    matrix[Nsub, Kw] Xk;
    matrix[N, 2] Xbatch;
    vector[Nt] obs[N];
    int mu_id[N];
    int p_re;
}
transformed data {
    matrix [3, 3] I = rep_matrix(0.0, 3, 3);
    for (j in 1:3)
        I[j,j] = 1.0;
}
parameters {
    real alpha0;
    real<lower=0> sigma;
    real<lower=0> k0;
    vector[Kw] k;
    vector[Kc] beta_cov;
    vector[Kt] beta_treat;
    vector[Kf] beta_flu;
    vector[2] beta_batch;
    vector[p_re] z[Nsub];
    simplex[p_re] pii;
    real<lower=0> tau;
    cholesky_factor_corr[p_re] L;
}
transformed parameters {
    matrix[p_re, p_re] psi = diag_pre_multiply(tau * pii, L);
    vector[N] batch_eff = Xbatch * beta_batch;
    vector[Nt] mu[Nsub];
    {
	    vector[3] re_i[Nsub];
	    vector[Nsub] eta = alpha0 + Xcov * beta_cov;
	    vector[Nsub] vax_boost = Xtreat * beta_treat;
	    vector[Nsub] flu_boost = Xflu * beta_flu;
	    vector[Nsub] ki = k0 * exp(Xk * k);
	    for (j in 1:Nsub) {
	        re_i[j] = L * z[j];
	        mu[j] = mu_calc(Nt, eta[j] + re_i[j][1], re_i[j][2] + vax_boost[j],
	                                flu_boost[j], tflu[j], exp(re_i[j][3]) * ki[j], tvec[j], mask[j,]);
    	}
    }
}
model {
    for (j in 1:N) {
        for (i in 1:5) {
            if (w[j][i] > 0) {
                target += w[j][i] * normal_lpdf(y[j][i] | mu[mu_id[j]][i] + batch_eff[j], sigma);
            }
        }
    }
    alpha0 ~ normal(0, 3);
    sigma ~ std_normal();
    tau ~ std_normal();
    k0 ~ std_normal();
    k ~ std_normal();
    beta_cov ~ normal(0, 5);
    beta_treat ~ normal(0, 5);
    beta_flu ~ normal(0, 5); 
    beta_batch ~ std_normal();
    // random effects
    L ~ lkj_corr_cholesky(1);
    tau ~ exponential(1);
    pii ~ dirichlet(rep_vector(1.0, p_re));
    for (j in 1:Nsub)
		z[j] ~ std_normal();
}
