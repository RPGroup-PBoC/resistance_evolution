functions{
    real logistic_growth(real t, real K, real lambda, real y_0){
        real num = K * y_0 * exp(lambda * t);
        real denom = K + y_0 * (exp(lambda * t) - 1);
        return num / denom;
    }

}

data{
    // Number of Data Points
    int N;
    // Time Stamps
    real t[N];
    // Measured OD
    real y[N];
    // Parameters for Priors
    real y_0_sigma;
    real K_sigma;
    real sigma_sigma;
    real log_lambda_mu;
    real log_lambda_sigma;
    // Number of time stamps for posterior predictive check
    int N_ppc;
    // Time Stamps for posterior predictive check
    real t_ppc[N_ppc];
}

parameters{
    // Carrying Capacity
    real <lower=0> K;
    // Growth Rate
    real log_lambda;
    // Initial OD
    real <lower=0> y_0;
    //
    real <lower=0>sigma;
}

transformed parameters {
   real lambda = 10^log_lambda; 
}

model{
    // Priors
    y_0 ~ normal(0, y_0_sigma);
    K ~ normal(0, K_sigma);
    sigma ~ normal(0, sigma_sigma);
    log_lambda ~ normal(log_lambda_mu, log_lambda_sigma);
    
    // Compute Theoritical Values
    real y_theo[N];
    for (i in 1:N){
        y_theo[i] = logistic_growth(t[i], K, lambda, y_0);
        y[i] ~ normal(y_theo[i], sigma);
    }

}

generated quantities {
    real y_predict[N_ppc];
    for (i in 1:N_ppc){
        y_predict[i] = normal_rng(logistic_growth(t_ppc[i], K, lambda, y_0), sigma);
    }
}