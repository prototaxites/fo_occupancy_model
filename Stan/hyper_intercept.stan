functions{
  real extended_hypergeometric_lpmf(int y, int m1, int m2, int n, real alpha) {
    int y_min = max(0, n - m2);
    int y_max = min(n, m1);
    int length = y_max - y_min + 1;
    vector[length] P;
    //real lomega = log(omega);
    int counter = 1;
    
    real ldist = lchoose(m1, y) + lchoose(m2, n - y) + y * alpha; 
    
    for (x in y_min:y_max) {
      P[counter] = lchoose(m1, x) + lchoose(m2, n - x) + x * alpha;
      counter += 1;
    }
  
    return ldist - log_sum_exp(P);
  }
}
data {
  // Counters
  int N;
  
  // Y
  array[N] int y;
  array[N] int mA;
  array[N] int mB;
  int n;
}
parameters {
  array[N] real alpha;
}
model {
  for(i in 1:N){
      target += normal_lpdf(alpha[i] | 0, 5);
      target += extended_hypergeometric_lpmf(y[i] | mA[i], mB[i], n, alpha[i]);
  }
}
