for(i in 1:N_species) {
  array[N_replicates, 2] int Y_spp = Y[,,i];

  int start=1; //replicates are ordered by site - starting at replicate 1
  for(j in 1:N_sites) {
    int end = start + rep_group_sizes[j] - 1; // add number of reps for site j and subtract 1

    if(Q[i, j] == 1) {
      target += bernoulli_logit_lpmf(1 | psi[j, i]) +
          sample_occupancy_sum_lpmf(Y_spp[start:end,] | lambda[start:end, i], theta[start:end, i], phi[i]);
    }
    if(Q[i, j] == 0) {
      target += log_sum_exp(bernoulli_logit_lpmf(1 | [j, i]) +
          sample_occupancy_sum_lpmf(Y_spp[start:end,] | lambda[start:end, i], theta[start:end, i], phi[i]),
          bernoulli_logit_lpmf(0 | psi[j, i]));
    }
    start += rep_group_sizes[j]; // update start point for next site
  }
}