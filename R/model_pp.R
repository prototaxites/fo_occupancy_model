# posterior_predict_phi <- function(ASV, draws, occ_data) {
#   draws <- posterior::as_draws_rvars(draws)
#   
#   with(draws, { occ_data %**% b_occ + occ_data %**% b_occ_species[, ASV] }) |> 
#     draws_of(with_chains = TRUE) |> 
#     plogis() |> 
#     rvar(with_chains = TRUE)
# }

posterior_predict_phi <- function(ASV, draws, env_data, sample_data, phi = NULL, theta = NULL) {
  draws <- posterior::as_draws_rvars(draws)
  
  b_env <- as.vector(draws$b_env[1,])
  b_env_species <- drop(draws$b_env_species[1,,])
  
  (env_data %**% b_env + env_data %**% b_env_species[, ASV]) |> 
    draws_of(with_chains = TRUE) |> 
    plogis() |> 
    rvar(with_chains = TRUE)
}

posterior_predict_theta <- function(ASV, draws, env_data, sample_data, phi = NULL,  theta = NULL) {
  draws <- posterior::as_draws_rvars(draws)
  
  if(is.null(phi)) {
    phi <- posterior_predict_phi(ASV, draws, env_data, sample_data) 
  }
  
  b <- c(as.vector(draws$b_env[2,]), as.vector(draws$b_sample[1,]))
  b_species <- rbind(drop(draws$b_env_species[2,,]), drop(draws$b_sample_species[1,,]))
  whole_data <- cbind(env_data, sample_data)
  
  ((whole_data %**% b + whole_data %**% b_species[, ASV]) |> 
    draws_of(with_chains = TRUE) |> 
    plogis() |> 
    rvar(with_chains = TRUE)) * phi
}

posterior_predict_lambda <- function(ASV, draws, env_data, sample_data, phi = NULL, theta = NULL) {
  draws <- posterior::as_draws_rvars(draws)
  
  if(is.null(theta)) {
    theta <- posterior_predict_theta(ASV, draws, env_data, sample_data, phi)
  }
  
  b <- c(as.vector(draws$b_env[3,]), as.vector(draws$b_sample[2,]))
  b_species <- rbind(drop(draws$b_env_species[3,,]), drop(draws$b_sample_species[2,,]))
  whole_data <- cbind(env_data, sample_data)
  
  exp(whole_data %**% b + whole_data %**% b_species[, ASV]) * theta
}

# posterior_predict_z <- function(ASV, draws, occ_data) {
#   draws <- posterior::as_draws_rvars(draws)
#   
#   with(draws, { occ_data %**% b_occ + occ_data %**% b_occ_species[, ASV] }) |> 
#     draws_of(with_chains = TRUE) |> 
#     plogis() |> 
#     apply(c(1, 2, 3, 4), FUN = \(x) rbinom(prob = x, size = 1, n = 1)) |>
#     rvar(with_chains = TRUE)
# }

# posterior_predict_u <- function(ASV, draws, Z, rep_data, rep_group_sizes) {
#   draws <- posterior::as_draws_rvars(draws)
#   Z <- draws_rvars(Z = Z)
#   draws <- bind_draws(draws, Z)
#   
#   N_sites <- length(rep_group_sizes)
#   
#   res <- NULL
#   start <- 1
#   for(i in 1:N_sites) {
#     end <- start + rep_group_sizes[i] - 1
#     site_rep_df <- rep_data[start:end,]
#     
#     site_res <- with(draws, {
#       Z[i] * (site_rep_df %**% b_rep + site_rep_df %**% b_rep_species[, ASV])
#     }) |> draws_of(with_chains = TRUE) |> 
#       plogis() |> 
#       apply(c(1, 2, 3, 4), FUN = \(x) rbinom(prob = x, size = 1, n = 1)) |>
#       rvar(with_chains = TRUE)
#     
#     start <- start + rep_group_sizes[i]
#     if(i == 1) {
#       res <- site_res
#     } else {
#       res <- rbind(res, site_res)
#     }
#   }
#   
#   return(res)
# }
# 
# posterior_predict_lambda <- function(ASV, draws, U, read_data, use_SF = FALSE) {
#   draws <- posterior::as_draws_rvars(draws)
#   U <- draws_rvars(U = U)
#   draws <- bind_draws(draws, U)
# 
#   if(use_SF == TRUE) {
#     SF <- 0
#   } else {
#     SF <- draws$SF
#   }
#   
#   with(draws, {
#     U * (SF + read_data %**% b_reads + read_data %**% b_reads_species[, ASV])
#   })
# }

big_rnbinom <- function(lambda, phi) {
  lambda <- lambda |> t() |> as.vector()
  rnbinom(n = length(lambda), mu = lambda, size = phi)
}

posterior_predict_counts <- function(ASV, draws, lambda) {
  draws <- posterior::as_draws_rvars(draws)
  
  phi <- draws$phi[ASV]
  lambda <- exp(lambda)
  n_iter <- niterations(lambda)
  n_chains <- nchains(lambda)
  n_samples <- nrow(lambda)
  n_species <- length(phi)
  
  arr_lambda <- lambda |> draws_of(with_chains = TRUE)
  arr_phi <- phi |> draws_of(with_chains = TRUE)
 
  map(1:n_iter, \(x) 
      map(1:n_chains, \(y)
        big_rnbinom(arr_lambda[x, y, , ], arr_phi[x, y, ])
      )
    ) |> 
  unlist() |>
  array(dim = c(n_species, n_samples, n_chains, n_iter)) |>
  aperm(c(4,3,2,1)) |>
  rvar(with_chains = TRUE)
}

