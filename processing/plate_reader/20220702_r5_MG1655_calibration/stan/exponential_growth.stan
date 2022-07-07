data{
    int N;
    vector[N] t;
    vector[N] y;

    int N_ppc;
    vector[N_ppc] t_ppc;
}

parameters{
    real log_10_tau;
    real<lower=0> y_0;
    real<lower=0> sigma_0;
}

transformed parameters {
    real tau = 10^log_10_tau;
    real lambda = log(2)/tau;
    vector[N] y_th;
    for (i in 1:N) 
        y_th[i] = log(y_0) + lambda * t[i];
}

model{
    log_10_tau ~ normal(1.75, 0.25);
    y_0 ~ normal(0, 0.025);
    sigma_0 ~ normal(0, 0.001);
    y ~ normal(y_th, sigma_0);
}

generated quantities {
   vector[N_ppc] y_ppc;
   for (i in 1:N_ppc)
      y_ppc[i] = log(y_0)+ lambda * t_ppc[i];
}