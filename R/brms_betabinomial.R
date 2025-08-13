## Beta binomial family for BRMS parameterised with an inverse phi dispersion parameter
## (lower phi = less dispersed) - easier to set priors on
## Identity link for both

beta_binomial2 <- brms::custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("identity", "identity"),
  lb = c(0, 0), ub = c(1, NA),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    real theta = (1 - phi) / phi;
    return beta_binomial_lpmf(y | T, mu * theta, (1 - mu) * theta);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    real theta = (1 - phi) / phi;
    return beta_binomial_rng(T, mu * theta, (1 - mu) * theta);
  }
"

stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

log_lik_beta_binomial2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  extraDistr::dbbinom(x = y, 
                      size = trials, 
                      alpha = mu * (1 / phi), 
                      beta = (1 - mu) * (1 / phi),
                      log = TRUE
  )
}

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  extraDistr::rbbinom(n = prep$ndraws, 
                      size = trials,
                      alpha = mu * (1 / phi), 
                      beta = (1 - mu) * (1 / phi)
  )
}

posterior_epred_beta_binomial2 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  mu
}

