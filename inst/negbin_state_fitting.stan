data {
  int<lower=1> n_species; // number of species
  int<lower=1> n_steps;// number of time points
  int<lower=1> n_clust;//number of clusters
  int<lower=1> n_person;//number of persons
  array[(n_steps*n_person), n_species] int<lower=0> y; //microbiom abundance 
}
parameters {
 vector<lower=0>[n_species] sigma_2;
 //vector<lower=0>[n_species] phi; 
 matrix<lower=0>[n_clust, n_species]phi;
 real<lower=0.1,upper = 10> dirich;
 real<lower = 0, upper = 1> scale;
  matrix[n_clust, n_species] beta;
  array[(n_steps*n_person)] simplex[n_clust] theta;
}
model {
  for (l in 1:n_clust) {
    for (k in 1:n_species) {
      beta[l, k] ~ normal(0, sigma_2[k]);
    }
  }
  dirich ~ gamma(5,scale);//can be changed to 10 // can make it a user choice
  for (i in 1:(n_steps*n_person)) {
    theta[i] ~ dirichlet(rep_vector(dirich, n_clust));//set to 0.1 will be easier to threshould
  }
  vector[n_species] p;
  vector[n_species] phi_2;
  for (i in 1:(n_steps*n_person)) {
    for (k in 1:n_species) {
      p = beta' * to_vector(theta[i]);
      phi_2 = phi' * to_vector(theta[i]);
      y[i, k] ~ neg_binomial_2_log(p[k], phi_2[k]);
    }
  }
}
