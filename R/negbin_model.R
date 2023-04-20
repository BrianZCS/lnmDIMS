#' Data normalization
#' @references https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02104-1#MOESM1
normDESeq2 <- function(physeq, whichOTUs = NULL, method = c("poscounts","ratio"))
{
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  } else {}# END - if: whichOTUs
  
  #   otu_table(physeq) <- otu_table(otuTab, taxa_are_rows = TRUE)
  
  ## Calculate size factors
  if (method == "poscounts")
  {
    obj <- phyloseq_to_deseq2(physeq,design = ~ 1)
    normFacts <- sizeFactors(DESeq2::estimateSizeFactors(obj,type = "poscounts"))
  } else {
    otuTab <- as(otu_table(physeq), "matrix")
    if (any(otuTab == 0))
    {
      otuTab <- otuTab + 1L
    } else {}
    normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
  }
}
  
  
  
#' Fit data to negative binomial clustering model
#' 
#' @param data_list a list containing number of species, number of steps, number of clusters,
#'           number of persons, the abundance data
#' @return a list containing the sequence of states, centers and phi
#' @export  
cal_fit_negbin <- function(data_list) {
  npp_model <- cmdstan_model(system.file("negbin_state_fitting", package = "lnmDIMS"))
  fit <-  npp_model$variational(data_list, iter = 10000, adapt_engaged = FALSE, eta = 0.1)
  a = fit$summary(variables = c("theta"), "mean")
  beta = fit$summary(variables = c("beta"),"mean")
  phi = fit$summary(variables = c("phi"),"mean")
  beta_matrix = matrix(nrow = data_list$n_clust, ncol = data_list$n_species)
  phi_matrix = matrix(nrow = data_list$n_clust, ncol = data_list$n_species)
  for (i in 1:(data_list$n_clust)) {
    for (j in 1:data_list$n_species) {
      beta_matrix[i, j] = beta$mean[beta$variable == paste("beta[", i, ",", j, "]", sep = "")]
      phi_matrix[i, j] = phi$mean[phi$variable == paste("phi[", i, ",", j, "]", sep = "")]
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
  return(list(theta = sequence, centers = beta_matrix, phi = phi_matrix))
}  

#' Simulate clustering data given a set of user defined parameters and model parameters from fitting existing data through negative bionomial model
#'
#' @param n_species number of species
#' @param ts_matrix transition matrix of the markov-chain
#' @param n_person number of persons 
#' @param count a vector of count of samples per person
#' @param initial_state initial state of the markov-chain
#' @param centers the centroids of each cluster
#' @param sigma variance added inside the cluster
#' @param person_effect difference among persons in the simulation
#' @return the simulated data
#' @export
sim_clust_negbin<-function(n_species, ts_matrix, n_timepoints, initial_state, centers, phi){
  states = c()
  
  for(i in 1:length(n_timepoints)){
    states = c(states,markov_sample(ts_matrix, n_timepoints[i], initial_state[i])) 
  }
  
  data=matrix(nrow=sum(n_timepoints), ncol=n_species)
  for (i in 1:sum(n_timepoints)){
    for  (k in 1:n_species){
      mu = exp(centers[states[i], k])
      data[i, k]= rnbinom(mu = mu, size = phi[states[i], k], n=1)
    }
  }
  data<-as.data.frame(data)
  
  return(data)
}



#' Make experimental configurations for power tests
#' 
#' @param  n_timepoints number of timepoints for each person
#' @param n_person number of persons in the experiment
#' @export  
make_configurations_negbin <- function(n_timepoints, n_person) {
  
  configurations <- cross(list(n_timepoints = n_timepoints, n_person = n_person))
  
  configurations
}


#' Conduct power tests and estimate the power of each configuration
#'
#' @param config configurations with user defined number of timepoints and number of samples
#' @param n_reps number of replicates for each experimental configurations
#' @param n_species number of species
#' @param ts_matrix transition matrix
#' @param initial_state a vector containing the initial states of each person
#' @param centers the centroids of clusters
#' @param phi dispersion parameter of the negative bionomial model
#' @export  
estimate_stat_negbin <- function(config, n_reps, n_species, ts_matrix, initial_state, centers, phi){
  
  statistics = matrix(nrow = length(configurations), ncol = n_reps)
  
  for (i in 1:length(configurations)) {
    for (j in 1:n_reps) {
      samples = data.frame()
      count = rep(configurations[[i]]$n_timepoints, configurations[[i]]$n_person)
      samples = sim_negbin(n_species = n_species,ts_matrix = ts_matrix, n_person = configurations[[i]]$n_person, count = count, initial_state = initial_state, centers = centers, phi = phi)
      samples<-as.matrix(samples)
      sample_list <- list(n_species = n_species, n_clust = 5, y = samples, n_person = configurations[[i]]$n_person, count = count)
      result=cal_fit_negbin(sample_list)
      theta = result$theta
      sample_state <- tibble(
        subject = rep(1:configurations[[i]]$n_person, count),
        cluster = theta
      ) %>%
        mutate(subject = as.factor(subject))  
      
      sample_matrix = estimate_transitions(sample_state)
      diff=mean((sample_matrix-ts_matrix)^2)
      statistics[i,j] = diff
    }
  }
  statistics
}
