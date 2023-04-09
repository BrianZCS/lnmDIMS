functions {
  vector phi_inv(vector mu, int K) {
    vector[K] e_mu;

    e_mu = exp(mu);
    return append_row(e_mu / (1 + sum(e_mu)), 1 / (1 + sum(e_mu)));
  }
}

data {
  int<lower=1> n_species; // number of species - 1 //Q1: why I need to use species - 1 sometimes
  int<lower=1> n_clust;//number of clusters
  int<lower=1> n_person;//number of persons
  array[n_person] int count;
  array[sum(count), n_species+1] int<lower=0> y; //microbiom abundance //change to array
}

parameters {
 vector<lower=0>[n_species]sigma_2;
 //matrix<lower=0>[n_clust, n_species] sigma; //coefficient for in/out clusters
 real<lower=0.1,upper = 10> dirich;
 real<lower = 0, upper = 1> scale;
  //real<lower=0> sigma;
  matrix[n_clust, n_species] beta; //coefficient for in/out clusters
  array[sum(count)] simplex[n_clust] theta; // in/out clusters
  matrix[n_person, n_species] person_effect;
}


model {
  vector[n_species + 1] p;
  for (l in 1:n_clust) {
    for (k in 1:n_species) {
      beta[l, k] ~ normal(0, sigma_2[k]);// mean and variance can be changed based on parameters
    }
  }

  dirich ~ gamma(5,scale);//can be changed to 10 // can make it a user choice
  for (i in 1:sum(count)) {
    theta[i] ~ dirichlet(rep_vector(dirich, n_clust));//set to 0.1 will be easier to threshould
  }

  vector[n_species] mu_;
  
  for (l in 1:n_person) {
    for (k in 1:n_species) {
      person_effect[l][k] ~ normal(0,1);// mean and variance can be changed based on parameters
    }
  }
  
  for (i in 1:n_person) {
    for(j in 1:count[n_person]){
      if(i == 1){
        mu_ = beta'* to_vector(theta[j]) + (person_effect[i])';
        p = phi_inv(mu_, n_species);
        y[j] ~ multinomial(to_vector(p)); 
      }
      else{
        mu_ = beta'* to_vector(theta[sum(count[1:i-1])+j]) + (person_effect[i])';
        p = phi_inv(mu_, n_species);
        y[sum(count[1:i-1])+j] ~ multinomial(to_vector(p)); 
      }
    }
  }
}




