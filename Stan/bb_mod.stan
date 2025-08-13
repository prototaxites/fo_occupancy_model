functions {
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    real theta = (1 - phi) / phi;
    return beta_binomial_lpmf(y | T, mu * theta, (1 - mu) * theta);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    real theta = (1 - phi) / phi;
    return beta_binomial_rng(T, mu * theta, (1 - mu) * theta);
  }
}
data {
  int N;
  int N_cat;
  array[N] int Y;
  array[N] int Total;
  array[N] int Categories;
  
  matrix[2, 2] cor;
  vector[2] sigma;
}
parameters {
  matrix[2, N_cat] b_phi_cat;
}
model {
  matrix[2, 2] cov = diag_matrix(sigma) * cor * diag_matrix(sigma);
  for(i in 1:N_cat) {
    target += multi_normal_lpdf(b_phi_cat[ ,i] | [0, 0], cov);
  }

  for(i in 1:N) {
    target += beta_binomial2_lpmf(Y[i] | inv_logit(b_phi_cat[1, Categories[i]]), inv_logit(b_phi_cat[2, Categories[i]]), Total[i]);
  }
}
generated quantities {
  vector[N_cat] epred;
  array[N_cat] int ppred;
  array[N_cat] int Y_rep;
  
  for(i in 1:N_cat) {
    epred[i] = inv_logit(b_phi_cat[1, i]);
    ppred[i] = beta_binomial2_rng(inv_logit(b_phi_cat[1, i]), inv_logit(b_phi_cat[2, i]), 10000);
    Y_rep[i] = beta_binomial2_rng(inv_logit(b_phi_cat[1, i]), inv_logit(b_phi_cat[2, i]), Total[i]);
  }
  
  vector[N] log_lik;
  for(i in 1:N) {
    log_lik[i] = beta_binomial2_lpmf(Y[i] | inv_logit(b_phi_cat[1, Categories[i]]), inv_logit(b_phi_cat[2, Categories[i]]), Total[i]);
  }
}
