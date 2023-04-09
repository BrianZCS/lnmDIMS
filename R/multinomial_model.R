phi_inverse <- function(mu) {
  sum_exp <- sum(exp(mu))
  c(exp(mu) / (1 + sum_exp), 1 / (1 + sum_exp))
}


#' Estimate a transition matrix from the given sequence of states
#'
#' @param clusters the alpha sequence of states
#' @export  
estimate_transitions <- function(clusters) {
  z <- clusters$cluster
  K <- n_distinct(z)
  P <- matrix(0, K, K)
  colnames(P) = unique(z)
  rownames(P) = unique(z)
  
  for (i in seq_len(nrow(clusters))) {
    if (i == nrow(clusters) || clusters$subject[i] != clusters$subject[i + 1]) next
    P[as.character(z[i]), as.character(z[i + 1])] <- P[as.character(z[i]), as.character(z[i + 1])] + 1
  }
  
  P / rowSums(P) 
}

#' Generate a sequence of states from the given transition matrix
#'
#' @param prob_matrix transition matrix among clusters
#' @param n_steps number of steps in teh markov chain
#' @param initial_state initial states for each person
#' @export  
markov_sample<-function(prob_matrix, n_steps, initial_state){
  
  result<-vector(length=n_steps)
  
  result[1]=initial_state
  for (i in 1:n_steps){
    
    result[i+1]=sample(x=1:nrow(prob_matrix),size=1, prob=prob_matrix[result[i],])
  }
  
  return(result)
}

#' 
#'
#' @export
cal_fit_clust_multinomial <- function(data_list) {
  npp_model <- cmdstan_model(system.file("lmn_state_fitting.stan", package = "lnmDIMS"))
  fit <-  npp_model$variational(data_list, iter = 10000, adapt_engaged = FALSE, eta = 0.1)
  a = fit$summary(variables = c("theta"), "mean")
  beta = fit$summary(variables = c("beta"),"mean")
  sigma = fit$summary(variables = c("sigma_2"),"mean")
  person_effect = fit$summary(variables = c("person_effect"),"mean")
  
  beta_matrix = matrix(nrow = data_list$n_clust, ncol = data_list$n_species)
  for (i in 1:(data_list$n_clust)) {
    for (j in 1:data_list$n_species) {
      beta_matrix[i, j] = beta$mean[beta$variable == paste("beta[", i, ",", j, "]", sep = "")]
    }
  }
  person_effect_matrix = matrix(nrow = data_list$n_person, ncol = data_list$n_species)
  for (i in 1:(data_list$n_person)) {
    for (j in 1:data_list$n_species) {
      person_effect_matrix[i, j] = person_effect$mean[person_effect$variable == paste("person_effect[", i, ",", j, "]", sep = "")]
    }
  }
  theta = data.frame()
  for (i in 1:(sum(data_list$count))) {
    for (j in 1:data_list$n_clust) {
      theta[i, j] = a$mean[a$variable == paste("theta[", i, ",", j, "]", sep = "")]
    }
  }
  sequence = c()
  for (k in 1:(sum(data_list$count))){
    sequence[k] =  which.max(theta[k, ])
  }
  
  temp_data <- tibble(
    subject = rep(1:data_list$n_person, data_list$count),
    cluster = sequence
  ) %>%
    mutate(subject = as.factor(subject))
  fit_matrix = estimate_transitions(temp_data)
  
  ##fix the problem of unvisited states
  if (nrow(fit_matrix) != data_list$n_clust || ncol(fit_matrix) != data_list$n_clust) {
    updated_matrix = matrix(nrow = data_list$n_clust, ncol = data_list$n_clust)
    for (j in 1:data_list$n_clust) {
      if (!(j %in% sequence)) {
        updated_matrix[1:ncol(updated_matrix), j] = 0
        updated_matrix[j, 1:ncol(updated_matrix)] = 1 / ncol(updated_matrix)
      }
    }
    for (j in rownames(fit_matrix)) {
      for (k in colnames(fit_matrix)) {
        updated_matrix[as.numeric(j), as.numeric(k)] = fit_matrix[j, k]
      }
    }
    return(list(theta = sequence, centers = beta_matrix, sigma = sigma$mean, person_effect = person_effect_matrix, ts_matrix = updated_matrix))
  }
  return(list(theta = sequence, centers = beta_matrix, sigma = sigma$mean, person_effect = person_effect_matrix, ts_matrix = fit_matrix))
}

#' Simulate clustering data given a set of user defined parameters and model parameters from fitting existing data
#'
#' @param n_depth sequence depth
#' @param n_species number of species
#' @param ts_matrix transition matrix of the markov-chain
#' @param n_person number of persons 
#' @param count a vector of count of samples per person
#' @param initial_state initial state of the markov-chain
#' @param centers the centroids of each cluster
#' @param sigma variance added inside the cluster
#' @param person_effect difference among persons in the simulation
#' 
#' @param obs pilot data
#' @return the simulated data
#' @export
sim_clust_multinomial<-function(n_depth, n_species, ts_matrix, n_person, count, initial_state = NULL, centers, sigma, person_effect){
  
  if(is.null(initial_state)){
    initial_state = sample(nrow(ts_matrix),n_person, replace = TRUE)
  }
  if(length(n_depth) == 1){
    n_depth = rep(n_depth, sum(count))
  }
  
  states = c()
  
  for(i in 1:length(count)){
    states = c(states,markov_sample(ts_matrix, count[i], initial_state[i])) 
  }
  
  data=matrix(nrow=sum(count), ncol=n_species)
  for (i in 1:n_person){
    effect = person_effect[sample(nrow(person_effect),1),]
    for(j in 1:count[i]){
      ## x = centers[states[i],] + mvrnorm(n=1, rep(0, n_species-1), diag(sigma[states[i],]))
      x = centers[states[i],] + effect
      if(i == 1){
        data[j,]=rmultinom(n=1, size = n_depth[sum(count[1:i-1])+j], prob = phi_inverse(x)) 
      }
      else{
        data[sum(count[1:i-1])+j,]=rmultinom(n=1, size = n_depth[sum(count[1:i-1])+j], prob = phi_inverse(x))  
      }
    }
  }
  data<-as.data.frame(data)
  
  return(data)
}
