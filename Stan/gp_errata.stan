parameters {
  vector<lower=0>[N_species] sigma_GP_occ;
  matrix<lower=0>[N_pred_occ, N_species] lscale_GP_occ;
  matrix[N_sites, N_species] z_gp_occ;
}
transformed parameters {
  matrix[N_sites, N_species] a_Site;
  for(j in 1:N_species){
    a_Site[,j] = gp(X_occ, sigma_GP_occ[j], lscale_GP_occ[,j], z_gp_occ[,j]);
  }
}
model {
  target += std_normal_lpdf(to_vector(z_gp_occ));
  target += exponential_lpdf(sigma_GP_occ | 1);
  target += inv_gamma_lpdf(to_vector(lscale_GP_occ) | 1.494197, 0.02481908);
}