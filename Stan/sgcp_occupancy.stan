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
                        matrix X_reads,
                        matrix X_rep,
                        matrix X_occ,
                        vector b_reads,
                        vector b_rep,
                        vector b_occ,
                        matrix b_reads_species,
                        matrix b_rep_species,
                        matrix b_occ_species,
                        vector phi,
                        vector SF,
                        array[] int rep_group_sizes) {
    
    int N = size(Y);
    real temp_target = 0;
    array[N, N_sites] int Q_local = Q[start:end,];
    vector[N] phi_local = phi[start:end];
    
    matrix[N_replicates, N] lambda = rep_matrix(SF, N) +
                                     rep_matrix(X_reads * b_reads, N) +  
                                     (X_reads * b_reads_species[,start:end]);
    matrix[N_replicates, N] theta = rep_matrix(X_rep * b_rep, N) + X_rep * b_rep_species[,start:end];
    matrix[N_sites, N] psi = rep_matrix(X_occ * b_occ, N) + X_occ * b_occ_species[,start:end];
    
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

  // Site occupancy predictors
  int N_pred_occ;                    // Number of predictors
  matrix[N_sites, N_pred_occ] X_occ; // Predictor matrix for occupancy
  
  // Replicate occupancy predictors
  int N_pred_rep;
  matrix[N_replicates, N_pred_rep] X_rep;

  // Read abundance predictors
  int N_pred_reads;
  matrix[N_replicates, N_pred_reads] X_reads;
  
  // Parallelisation options
  int grainsize;
  
  // Priors
  real prior_sigma;
  real prior_b;
  real prior_ljk;
}
parameters {
  // Site occupancy parameters
  vector[N_pred_occ] b_occ;
  cholesky_factor_corr[N_pred_occ] L_Rho_occ;
  vector<lower=0>[N_pred_occ] sigmas_occ;
  matrix[N_pred_occ, N_species] z_occ;

  // Sample occupancy parameters
  vector[N_pred_rep] b_rep;
  cholesky_factor_corr[N_pred_rep] L_Rho_rep;
  vector<lower=0>[N_pred_rep] sigmas_rep;
  matrix[N_pred_rep, N_species] z_rep;
  
  // Size factors/row effects
  vector<lower=0>[N_replicates] SF;

  // Read abundance parameters
  vector[N_pred_reads] b_reads;
  cholesky_factor_corr[N_pred_reads] L_Rho_reads;
  vector<lower=0>[N_pred_reads] sigmas_reads;
  matrix[N_pred_reads, N_species] z_reads;
  
  // Read dispersion
  vector[N_species] z_phi;
  real<lower=0> sigma_phi;
  real a_phi;
}
transformed parameters {
  // Site occupancy parameters
  matrix[N_pred_occ, N_species] b_occ_species;
  b_occ_species = (diag_pre_multiply(sigmas_occ, L_Rho_occ) * z_occ);
  
  // Sample occupancy parameters
  matrix[N_pred_rep, N_species] b_rep_species;
  b_rep_species = (diag_pre_multiply(sigmas_rep, L_Rho_rep) * z_rep);

  // Read abundance parametters
  matrix[N_pred_reads, N_species] b_reads_species;
  b_reads_species = (diag_pre_multiply(sigmas_reads, L_Rho_reads) * z_reads);
  
  // Read dispersion
  vector<lower=0>[N_species] phi;
  phi = exp((a_phi + z_phi) * sigma_phi);
}
model {
  // Priors
  // Site occupancy
  target += normal_lpdf(b_occ | 0, prior_b);
  target += lkj_corr_cholesky_lpdf(L_Rho_occ | prior_corr);
  target += exponential_lpdf(sigmas_occ | prior_sigma);
  target += std_normal_lpdf(to_vector(z_occ));

  // Sample occupancy
  target += normal_lpdf(b_occ | 0, prior_b);
  target += lkj_corr_cholesky_lpdf(L_Rho_rep | prior_corr);
  target += exponential_lpdf(sigmas_rep | prior_sigma);
  target += std_normal_lpdf(to_vector(z_rep));
  
  // Size factors
  target += lognormal_lpdf(SF | 0, 1);

  // Read abundance 
  target += std_normal_lpdf(b_reads);
  target += lkj_corr_cholesky_lpdf(L_Rho_reads | prior_corr);
  target += exponential_lpdf(sigmas_reads | prior_sigma);
  target += std_normal_lpdf(to_vector(z_reads));

  // Read dispersion
  target += std_normal_lpdf(a_phi);
  target += std_normal_lpdf(z_phi);
  target += exponential_lpdf(sigma_phi | 1);

  //array[,,] int Y, int grainsize, array[,] int Q, int N_sites, int N_replicates, matrix X_reads,
  //matrix X_rep,   matrix X_occ, vector b_reads, vector b_rep, vector b_occ, matrix b_reads_species,
  //matrix b_rep_species, matrix b_occ_species, vector phi, vector SF, array[] int rep_group_sizes
  
  target += reduce_sum(partial_sum_lpmf, Y, grainsize, Q, N_sites, N_replicates, X_reads, X_rep, X_occ, b_reads,
                       b_rep, b_occ, b_reads_species, b_rep_species, b_occ_species, phi, SF, rep_group_sizes);

}
generated quantities {
  matrix[N_sites, N_species] log_lik;

  {
    matrix[N_replicates, N_species] lambda = rep_matrix(SF, N_species) +
                                            rep_matrix(X_reads * b_reads, N_species) +
                                            (X_reads * b_reads_species);
    matrix[N_replicates, N_species] theta = rep_matrix(X_rep * b_rep, N_species) + X_rep * b_rep_species;
    matrix[N_sites, N_species] psi = rep_matrix(X_occ * b_occ, N_species) + X_occ * b_occ_species;
    
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
