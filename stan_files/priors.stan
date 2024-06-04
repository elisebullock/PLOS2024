data{
    real stdcounts;
    real meancounts;
    real numObs;
    real scalescale;
    
    real alpha_A_data;
    real alpha_B_data;
    real delta_A_data;
    real delta_B_data;
    real gamma_data;
    real fs_data;
    real beta_data;
    real Source_data;
    real eff_data;
    real mu_data;
}

parameters{

    // priors for parameters
    real<lower = gamma_data/10,     upper = gamma_data*100>    gamma;
    real<lower = delta_A_data/2,    upper = delta_A_data*2>   delta_A; //death of fast
    real<lower = alpha_A_data/2,    upper = delta_A>          alpha_A; //division of fast
    real<lower = alpha_B_data/2,    upper = alpha_A>          alpha_B; //division of slow
    real<lower = alpha_B,           upper = delta_A>          delta_B; //death of slow
    real<lower = fs_data/2,         upper = 0.95>             fs;
    real<lower = beta_data/2,       upper = beta_data*2>      beta;
    real<lower = Source_lo,       upper = Source_hi>    Source;
    real<lower = eff_data/2,        upper = 1>                eff;
    real<lower = mu_data/2,         upper = mu_data*2>        mu;

    real<lower = 0>         scale;
    real<lower = 0>         counts_t0;
    real<lower = 0>         phi_inv;
}

model{
    
    //priors based on 4cm data
    alpha_A ~ lognormal(log(alpha_A_data),  0.5);
    alpha_B ~ lognormal(log(alpha_B_data),  1);
    delta_A ~ lognormal(log(delta_A_data),  0.5);
    delta_B ~ lognormal(log(delta_B_data),  1);
    beta    ~ lognormal(log(beta_data),     0.1);
    gamma   ~ lognormal(log(gamma_data),    1);
    fs      ~ lognormal(log(fs_data),       0.5);
    mu      ~ lognormal(log(mu_data),       0.1);
    eff     ~ lognormal(log(eff_data),      0.1);
    Source  ~ lognormal(log(Source_data),   0.5);
    
    scale ~ lognormal(log(stdcounts), scalescale);
    phi_inv ~ lognormal(log(0.5),scalescale);
    counts_t0  ~ lognormal(log(meancounts),    stdcounts*10/numObs);
}     

