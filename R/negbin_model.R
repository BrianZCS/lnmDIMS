# library(lnmDIMS)
# library(tidyverse)
# library(cmdstanr)
# library(dplyr)
# library(MASS)
# library(MCMCpack)
# library(markovchain)
# library(tidybayes)
# library(phyloseq)
# library(DESeq2)

#' Data normalization
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
  
  
  
# Fit data to negative binomial clustering model
# data_list a list containing number of species, number of clusters,
#           number of persons, and a vector containing number of samples for each person
# return theta the sequence of states for the input data, 
#        centers the matrix representing the contribution of each species in each cluster
#        phi the matrix containing all phi parameters for each species in each cluster
cal_fit_negbin <- function(data_list) {
  npp_model <- cmdstan_model("./Negbin2.stan")
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


# Make experimental configurations for power tests
# n_timepoints number of timepoints for each person
# n_person number of persons in the experiment
make_configurations_negbin <- function(n_timepoints, n_person) {
  
  configurations <- cross(list(n_timepoints = n_timepoints, n_person = n_person))
  
  configurations
}


# Conduct power tests and estimate the power of each configuration
# config configurations with user defined number of timepoints and number of samples
# n_reps number of replicates for each experimental configurations
# n_species number of species
# ts_matrix transition matrix
# initial_state a vector containing the initial states of each person
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
