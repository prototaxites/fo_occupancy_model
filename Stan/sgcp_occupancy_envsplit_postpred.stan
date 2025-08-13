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
  
  int number_non_zero(array[] int ints) {
    int N = size(ints);
    int count = 0;
    for(i in 1:N) {
      if(i > 0) count += 1;
    }
    
    return count;
  }
  
  real shannon(array[] int theta) {
    int N = size(theta);
    vector[N] vec = to_vector(theta);
    
    return sum(vec .* log(vec));
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
}
generated quantities {
  array[N_sites, N_species] int Z;
  array[N_sites, N_species] int Z_cond;
  array[N_sites, N_species] int Z_cond_norm;
  array[N_replicates, N_species] int U;
  array[N_replicates, N_species] int U_cond;
  array[N_replicates, N_species] int U_cond_norm;
  array[N_replicates, N_species] int Y_rep;
  array[N_replicates, N_species] int Y_norm_cond;

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
    matrix[N_replicates, N_species] lambda_norm = rep_matrix(X_env_sample * b_env[3], N_species) +  
                                     rep_matrix(X_sample * b_sample[2], N_species) +  
                                     (X_env_sample * b_env_species[3]) +
                                     (X_sample * b_sample_species[2]);
  
    for(i in 1:N_species) {
      Z[,i] = bernoulli_logit_rng(psi[,i]);
  
      int start=1; //replicates are ordered by site - starting at replicate 1
      for(j in 1:N_sites) {
        int end = start + rep_group_sizes[j] - 1; // add number of reps for site j and subtract 1
        
        // Calculate conditional occupancy at site level
        if(Q[i, j] == 1) {
          Z_cond[j, i] = 1;
        } else {
          real log_p_occ = bernoulli_logit_lpmf(1 | psi[j, i]) + 
              sample_occupancy_sum_lpmf(Y[i,start:end] | lambda[start:end, i], theta[start:end, i], phi[i]);
          real log_p_nonocc = log_sum_exp(bernoulli_logit_lpmf(1 | psi[j, i]) + 
              sample_occupancy_sum_lpmf(Y[i, start:end] | lambda[start:end, i], theta[start:end, i], phi[i]),
              bernoulli_logit_lpmf(0 | psi[j, i]));
          real log_p_cond = log_p_occ - log_p_nonocc;
          Z_cond[j, i] = bernoulli_rng(exp(log_p_cond));
        }
        
        U[start:end, i] = bernoulli_logit_rng(theta[start:end,i]);
        for(k in start:end) {
          // Conditional occupancy of U
          if(Q[i, j] == 1) {
            if(Y[i, k] > 0) {
              U_cond[k, i] = 1;
            } else {
              real log_p_occ = sample_detection_lpmf(0 | lambda[k, i], theta[k, i], phi[i]);
              real log_p_nonocc = sample_nondetection_lpmf(0 | lambda[k, i], theta[k, i], phi[i]);
              real log_p_cond = log_p_occ - log_p_nonocc;
              U_cond[k, i] = bernoulli_rng(exp(log_p_cond));
            }
          } else {
            real log_p_occ = sample_detection_lpmf(0 | lambda[k, i], theta[k, i], phi[i]);
            real log_p_nonocc = sample_nondetection_lpmf(0 | lambda[k, i], theta[k, i], phi[i]);
            real log_p_cond = log_p_occ - log_p_nonocc;
            U_cond[k, i] = bernoulli_rng(exp(Z_cond[j, i] * log_p_cond));
          }
          
          Y_rep[k, i] = neg_binomial_2_log_rng(lambda[k, i] * U[k, i], phi[i]);
          Y_norm_cond[k, i] = neg_binomial_2_log_rng(lambda_norm[k, i] * U_cond[k, i], phi[i]);
        }
  
        start += rep_group_sizes[j]; // update start point for next site
      }
    }
  }
  
  array[N_sites] int Z_cond_richness;
  array[N_replicates] int U_cond_richness;

  for(i in 1:N_sites){
    Z_cond_richness[i] = sum(Z_cond[i,]);
  }
  
  for(i in 1:N_replicates){
    U_cond_richness[i] = sum(U_cond[i,]);
  }
  
  corr_matrix[N_pred_env] Rho_env = multiply_lower_tri_self_transpose(L_Rho_env);
  corr_matrix[N_pred_sample] Rho_sample = multiply_lower_tri_self_transpose(L_Rho_sample);
}
