functions {
  real sample_detection_lpmf(int Y, real lambda, real theta, real phi){
    return bernoulli_logit_lpmf(1 | theta) +
      neg_binomial_2_log_lpmf(Y | lambda, phi);
  }
  
  real sample_nondetection_lpmf(int Y, real lambda, real theta, real phi) {
    return log_sum_exp(
      sample_detection_lpmf(Y | lambda, theta, phi), bernoulli_logit_lpmf(0 | theta));
  }
  
  real sample_occupancy_sum_lpmf(array[] int Y, vector lambda, vector theta, real phi) {
    int K = size(Y);
    real temp_target = 0;
    
    for(k in 1:K) {
      if(Y[k] > 0) {
        temp_target += sample_detection_lpmf(Y[k] | lambda[k], theta[k], phi);
      } else {
        temp_target += sample_nondetection_lpmf(0 | lambda[k], theta[k], phi);
      }
    }
    
    return(temp_target);
  }
  
  real partial_sum_lpmf(array[,] int Y, 
                        int start, int end,
                        array[,] int Q,
                        int N_sites,
                        int N_replicates,
                        matrix X_env,
                        matrix X_env_sample,
                        matrix X_sample,
                        array[] vector b_env,
                        array[] vector b_sample,
                        array[] matrix b_env_species,
                        array[] matrix b_sample_species,
                        vector phi,
                        vector SF,
                        array[] int rep_group_sizes) {
    
    int N = size(Y);
    real temp_target = 0;
    array[N, N_sites] int Q_local = Q[start:end,];
    vector[N] phi_local = phi[start:end];
    
    matrix[N_sites, N] psi = rep_matrix(X_env * b_env[1], N) + 
                             X_env * b_env_species[1,,start:end];
    matrix[N_replicates, N] theta = rep_matrix(X_env_sample * b_env[2], N) + 
                                    rep_matrix(X_sample * b_sample[1], N) +
                                    (X_env_sample * b_env_species[2,,start:end]) +
                                    (X_sample * b_sample_species[1,,start:end]);
    matrix[N_replicates, N] lambda = rep_matrix(SF, N) +
                                     rep_matrix(X_env_sample * b_env[3], N) +  
                                     rep_matrix(X_sample * b_sample[2], N) +  
                                     (X_env_sample * b_env_species[3,,start:end]) +
                                     (X_sample * b_sample_species[2,,start:end]);
                                    
    for(i in 1:N) {
      int loc_start=1; //replicates are ordered by site - starting at replicate 1
      for(j in 1:N_sites) {
        int loc_end = loc_start + rep_group_sizes[j] - 1; // add number of reps for site j and subtract 1
  
        if(Q_local[i, j] == 1) {
          temp_target += bernoulli_logit_lpmf(1 | psi[j, i]) + 
              sample_occupancy_sum_lpmf(Y[i,loc_start:loc_end] | lambda[loc_start:loc_end, i], theta[loc_start:loc_end, i], phi_local[i]);
        }
        if(Q_local[i, j] == 0) {
          temp_target += log_sum_exp(bernoulli_logit_lpmf(1 | psi[j, i]) + 
              sample_occupancy_sum_lpmf(Y[i, loc_start:loc_end] | lambda[loc_start:loc_end, i], theta[loc_start:loc_end, i], phi_local[i]),
              bernoulli_logit_lpmf(0 | psi[j, i]));
        }
        loc_start += rep_group_sizes[j]; // update start point for next site
      }
    }
    
    return temp_target;
  }
}
data {
  int N_species;                            // Number of species
  int N_sites;                              // Number of sites
  int N_replicates;                         // Number of replicates
  array[N_species, N_replicates] int Y;     // Raw data [OTU][Replicate]
  array[N_species, N_sites] int Q;          // Site occupation factor
  
  // Different number of replicates per site, so the input replicates are passed in site
  // order - then we can access them per site by looping over the group sizes
  array[N_sites] int rep_group_sizes;

  // Replicate occupancy predictors
  int N_pred_env;
  matrix[N_sites, N_pred_env] X_env;

  // Read abundance predictors
  int N_pred_sample;
  matrix[N_replicates, N_pred_sample] X_sample;
  
  // Parallelisation options
  int grainsize;
}
transformed data {
  matrix[N_replicates, N_pred_env] X_env_sample;
  {
    int start=1;
    for(i in 1:N_sites) {
      int end = start + rep_group_sizes[i] - 1;
      for(j in start:end) {
        X_env_sample[j,] = X_env[i,];
      }
      start += rep_group_sizes[i]; // update start point for next site
    }
  }
}
parameters {
  // Env parameters
  vector<lower=0>[N_pred_env] sigmas_env;
  array[N_pred_env] simplex[3] props_env;
  cholesky_factor_corr[N_pred_env] L_Rho_env;
  array[3] matrix[N_pred_env, N_species] z_env;
  array[3] vector[N_pred_env] b_env;

  // Sample parameters
  vector<lower=0>[N_pred_sample] sigmas_sample;
  array[N_pred_sample] simplex[2] props_sample;
  cholesky_factor_corr[N_pred_sample] L_Rho_sample;
  array[2] matrix[N_pred_sample, N_species] z_sample;
  array[2] vector[N_pred_sample] b_sample;
  vector<lower=0>[N_replicates] SF;

  // Read dispersion
  vector[N_species] z_phi;
  real<lower=0> sigma_phi;
  real a_phi;
}
transformed parameters {
  array[3] matrix[N_pred_env, N_species] b_env_species;
  array[2] matrix[N_pred_sample, N_species] b_sample_species;
  matrix[N_pred_env, 3] sigmas_env_split;
  matrix[N_pred_sample, 2] sigmas_sample_split;
  
  for(i in 1:N_pred_env) {
    sigmas_env_split[i,] = (sqrt(props_env[i]) * sqrt(3) * sigmas_env[i])';
  }
  
  for(i in 1:N_pred_sample) {
    sigmas_sample_split[i,] = (sqrt(props_sample[i]) * sqrt(2) * sigmas_sample[i])';
  }
  
  for(i in 1:3) {
    b_env_species[i] = (diag_pre_multiply(sigmas_env_split[,i], L_Rho_env) * z_env[i]);
    if(i != 3) {
      b_sample_species[i] = (diag_pre_multiply(sigmas_sample_split[,i], L_Rho_sample) * z_sample[i]);
    }
  }
  
  // Read dispersion
  vector<lower=0>[N_species] phi;
  phi = exp((a_phi + z_phi) * sigma_phi);
}
model {
  // Priors
  // Environment
  target += exponential_lpdf(sigmas_env | 2);
  target += lkj_corr_cholesky_lpdf(L_Rho_env | 3);
  for(i in 1:N_pred_env) {
    target += dirichlet_lpdf(props_env[i] | rep_vector(1, 3));
  }
  for(i in 1:3) {
    target += std_normal_lpdf(to_vector(z_env[i]));
    target += std_normal_lpdf(to_vector(b_env[i]));
  }
  
  //Sample
  target += exponential_lpdf(sigmas_sample | 2);
  target += lkj_corr_cholesky_lpdf(L_Rho_sample | 3);
  for(i in 1:N_pred_sample) {
    target += beta_lpdf(props_sample[i] | 1, 1);
  }
  for(i in 1:2) {
    target += std_normal_lpdf(to_vector(z_sample[i]));
    target += std_normal_lpdf(to_vector(b_sample[i]));
  }

  // Size factors
  target += lognormal_lpdf(SF | 0, 1);

  // Read dispersion
  target += std_normal_lpdf(a_phi);
  target += std_normal_lpdf(z_phi);
  target += exponential_lpdf(sigma_phi | 1);

  //array[,,] int Y, int grainsize, array[,] int Q, int N_sites, int N_replicates, matrix X_env,
  //matrix X_env_samples, matrix X_samples, vector b_env, vector b_samples, matrix b_env_species,
  //matrix b_samples_species, vector phi, vector SF, array[] int rep_group_sizes
  
  target += reduce_sum(partial_sum_lpmf, Y, grainsize, Q, N_sites, N_replicates, X_env, X_env_sample, X_sample, 
                       b_env, b_sample, b_env_species, b_sample_species, phi, SF, rep_group_sizes);

}
generated quantities {
  matrix[N_sites, N_species] log_lik;

  {
    matrix[N_sites, N_species] psi = rep_matrix(X_env * b_env[1], N_species) + 
                             X_env * b_env_species[1];
    matrix[N_replicates, N_species] theta = rep_matrix(X_env_sample * b_env[2], N_species) + 
                                    rep_matrix(X_sample * b_sample[1], N_species) +
                                    (X_env_sample * b_env_species[2]) +
                                    (X_sample * b_sample_species[1]);
    matrix[N_replicates, N_species] lambda = rep_matrix(SF, N_species) +
                                     rep_matrix(X_env_sample * b_env[3], N_species) +  
                                     rep_matrix(X_sample * b_sample[2], N_species) +  
                                     (X_env_sample * b_env_species[3]) +
                                     (X_sample * b_sample_species[2]);
    
    for(i in 1:N_species) {
      int start=1; //replicates are ordered by site - starting at replicate 1
      for(j in 1:N_sites) {
        int end = start + rep_group_sizes[j] - 1; // add number of reps for site j and subtract 1
        if(Q[i, j] == 1) {
          log_lik[j, i] = bernoulli_logit_lpmf(1 | psi[j, i]) + 
              sample_occupancy_sum_lpmf(Y[i,start:end] | lambda[start:end, i], theta[start:end, i], phi[i]);
        }
        if(Q[i, j] == 0) {
          log_lik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | psi[j, i]) + 
              sample_occupancy_sum_lpmf(Y[i, start:end] | lambda[start:end, i], theta[start:end, i], phi[i]),
              bernoulli_logit_lpmf(0 | psi[j, i]));
        }
        start += rep_group_sizes[j]; // update start point for next site
      }
    }
  }
}
