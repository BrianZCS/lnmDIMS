## multinomial model


#' configurations generator for multinomial model
#'
#' The function generate all configurations used in the power test
#' @param n_person the number of persons
#' @param n_depth sequence depth
#' @param n_timepoints number of time points
#' @return configurations the configurations used in the power test
#' @importFrom purrr cross
#' @export
make_configurations_mult <- function(n_person, n_depth, n_timepoints) {
  
  configurations <- cross(list(n_timepoints = n_timepoints, n_person = n_person, n_depth = n_depth))
  
  configurations
}

#' calculate the powers
#'
#' The function calculates the statistics of power test
#' @param config configurations
#' @param n_reps number of repetitions
#' @param ts_matrix transition matrix
#' @param initial_state initial state
#' @param centers centroid of cluster
#' @param person_effect difference among persons in the simulation
#' @return  power statistics
#' @export
estimate_stat_mult <- function(config, n_reps, ts_matrix, initial_state, centers, person_effect){
  
  statistics = matrix(nrow = length(configurations), ncol = n_reps)
  n_species = ncol(centers)+1
  for (i in 1:length(configurations)) {
    for (j in 1:n_reps) {
      samples = data.frame()
      count = rep(configurations[[i]]$n_timepoints, configurations[[i]]$n_person)
      samples = sim_clust(n_species = n_species,ts_matrix = ts_matrix, n_person = configurations[[i]]$n_person, count = count, initial_state = initial_state, centers = centers, person_effect = person_effect, n_depth = configurations[[i]]$n_depth)
      samples<-as.matrix(samples)
      sample_list <- list(n_species = n_species-1, n_clust = 5, y = samples, n_person = configurations[[i]]$n_person, count = count)
      result=cal_fit_clust(sample_list)
      sample_matrix = result$matrix
      diff=sum((sample_matrix-ts_matrix)^2)
      statistics[i,j] = diff
    }
  }
  statistics
}

#' transform the matrix to the dataframe
#'
#' The function transform the power statistics to the dataframe
#' @param matrix matrix of the power statistics
#' @return transformed dataframe
#' @export
power_matrix_transform = function(matrix){
  stat = as.data.frame(matrix)
  for(i in 1:ncol(stat)){
    colnames(stat)[i] = paste("rep" , i, sep = "")
  }
  for(i in 1:length(configurations)){
    stat$n_person[i] = configurations[[i]]$n_person
    stat$n_timepoints[i] = configurations[[i]]$n_timepoints
    stat$n_depth[i] = configurations[[i]]$n_depth
  }
  pivot_longer(stat, cols = starts_with("rep"),names_to = "rep",values_to = "power")
}





#' Make experimental configurations for power tests
#' 
#' Make experimental configurations for power tests for negative bionomial model
#' @param  n_timepoints number of timepoints for each person
#' @param n_person number of persons in the experiment
#' @export  
make_configurations_negbin <- function(n_timepoints, n_person) {
  
  configurations <- cross(list(n_timepoints = n_timepoints, n_person = n_person))
  
  configurations
}


#' Conduct power tests and estimate the power of each configuration
#'
#' Conduct power tests and estimate the power of each configuration for negative bionomial model
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
      samples = sim_clust_negbin(n_species = n_species,ts_matrix = ts_matrix, n_person = configurations[[i]]$n_person, count = count, initial_state = initial_state, centers = centers, phi = phi)
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




#' Alpha sequence generator
#'
#' The function generate alpha (perturbation intensity sequence) which can be used in the simulator
#' @param n_timepoints number of time points
#' @param interval the interval of perturbation (interval[1]-interval[2]): increase the permutation intensity.(interval[2]-interval[3]): achieve maximized perturbation intensity. (interval[3]-interval[4]): decrease the permutation intensity.
#' @param signal_size the maximum perturbation intensity (between 0 and 1)
#' @return  alpha the perturbation intensity seuqence, 0 means no perturbed, 1 means fully perturbed
#' @export
make_alphas <- function(n_timepoints, intervals, signal_size) {
  alpha <- rep(0, n_timepoints)
  alpha[intervals[1]:intervals[2]] <- seq(0, signal_size, length.out = intervals[2] - intervals[1] + 1)
  alpha[intervals[2]:intervals[3]] <- signal_size
  alpha[intervals[3]:intervals[4]] <- seq(signal_size, 0, length.out = intervals[4] - intervals[3] + 1)
  alpha
}

#' configurations generator for perturbation
#'
#' The function generate all configurations used in the power test
#' @param ns the number of samples
#' @param n_depth sequence depth
#' @param alpha the maximum perturbation intensity (between 0 and 1)
#' @param interval the interval of perturbation (interval[1]-interval[2]): increase the permutation intensity.(interval[2]-interval[3]): achieve maximized perturbation intensity. (interval[3]-interval[4]): decrease the permutation intensity.
#' @param n_timepoints number of time points
#' @return configurations the configurations used in the power test
#' @importFrom purrr cross
#' @export
make_configurations_perturb <- function(ns, n_depth, alpha, intervals = NULL, n_timepoints = NULL) {
  if (is.null(intervals)) {
    intervals <- c(5, 8, 14, 18)
  }
  if (is.null(n_timepoints)) {
    n_timepoints <- 20
  }
  configurations <- cross(list(n = ns, n_depth = n_depth, alpha = alpha))
  for (i in seq_along(configurations)) {
    configurations[[i]][["alpha"]] <- make_alphas(n_timepoints, intervals, configurations[[i]][["alpha"]]) %>%
      round(3)
  }
  configurations
}


#' Get the centroids of the fitted model
#'
#' The function get the centroids of the fitted model
#' @param data_list list used in the fitting list(alpha=alpha, n_species=ncol(obs)-1, n_timepoints = length(alpha), y = obs))
#' @return  fitted centroids (beta)
#' @export
get_centroids <- function(data_list){
  result = cal_fit_perturbation(data_list)
  result$beta
}

#' Modify centroids
#'
#' The function modifies the centroids to perturbed/unperturbed
#' @param centroids the centorids of the cluster
#' @return  modified centroids
#' @export
modify_centroids_fraction <- function(centroids, fraction_perturbed = NULL) {
  # estimated centroids from pilot # get_centroids(...)
  # sort effects from highest to lowest (reorder columns)
  # find threshold for fraction_perturbed
  # all columns after fraction_perturbed -> replace with midpoint of two centroids
  if(is.null(fraction_perturbed)){
    fraction_perturbed = 0.5
  }

  difference = sort(rank(abs(centroids[2,] - centroids[1,])))
  for(i in 1:(ncol(centroids)*(1-fraction_perturbed))){
    column = as.integer(substring(names(difference[i]),2))
    centroids[1,column] = (centroids[1,column] + centroids[2,column])/2
    centroids[2,column] = centroids[1,column]
  }
  centroids
}

#' transform the matrix to the dataframe
#'
#' The function transform the power statistics to the dataframe
#' @param matrix matrix of the power statistics
#' @return transformed dataframe
#' @export
power_matrix_transform_perturb = function(matrix){
  stat = as.data.frame(matrix)
  for(i in 1:ncol(stat)){
    colnames(stat)[i] = paste("rep" , i, sep = "")
  }
  for(i in 1:length(configurations)){
    stat$n[i] = configurations[[i]]$n
    stat$n_depth[i] = configurations[[i]]$n_depth
    stat$alpha[i] = max(configurations[[i]]$alpha) ##signal_size, need to be changed if the graph changed.
  }
  pivot_longer(stat, cols = starts_with("rep"),names_to = "rep",values_to = "power")
}

#' calculate the powers
#'
#' The function calculates the statistics of power test
#' @param config configurations
#' @param centroids centroids of the cluster
#' @param n_reps number of repetitions
#' @param simulator simulator we use
#' @param sig_level significance level
#' @param sim_opts true perturbed species
#' @return  power statistics
#' @export
estimate_power_perturb <-
  function(configurations,
           centroids,
           n_reps,
           simulator,
           sig_level,
           sim_opts) {
    ##Take average of samples, then do t-test
    statistics = matrix(nrow = length(configurations), ncol = n_reps)
    for (i in 1:length(configurations)) {
      for (j in 1:n_reps) {
        samples = data.frame()
        for (k in 1:configurations[[i]]$n) {
          samples = rbind(samples, simulator(
            alpha = configurations[[i]]$alpha,
            n_depth = configurations[[i]]$n_depth,
            beta = centroids
          ))
        }
        detected = 0
        truth = 0
        for (k in sim_opts) {
          perturbed = vector()
          unperturbed = vector()
          for (l in 1:(length(configurations[[i]]$alpha)*configurations[[i]]$n)) {
            if (configurations[[i]]$alpha[((l-1)%%length(configurations[[i]]$alpha))+1] == 0) {
              unperturbed = list.append(unperturbed, samples[l, k])
            }
            else{
              perturbed = list.append(perturbed, samples[l, k])
            }
          }
          truth = truth + 1
          if (t.test(perturbed, unperturbed)$p.value < sig_level) {
            detected = detected + 1
          }
        }
        statistics[i, j] = detected / truth
      }
    }
    statistics
  }





#' #' calculate the powers
#' #'
#' #' The function calculates the statistics of power test
#' #' @param configurations configurations
#' #' @param n_reps number of repetitions
#' #' @param simulator simulator we use
#' #' @param ts_matrix transition matrix
#' #' @param initial_state initial state
#' #' @param centers centroid of cluster
#' #' @return  power statistics
#' #' @export
#' estimate_power_clust <- function(configurations, n_reps, simulator, ts_matrix, inital_state, centers){
#'   n_species = ncol(centers)+1
#'   n_clust = nrow(centers)
#'   statistics = matrix(nrow = length(configurations), ncol = length(n_reps))
#'   for (i in 1:length(configurations)){
#'     
#'     statistics = matrix(nrow = length(configurations), ncol = n_reps)
#'     for (i in 1:length(configurations)) {
#'       for (j in 1:n_reps) {
#'         diff<-vector()
#'         for (k in 1:configurations[[i]]$n){
#'           samples = data.frame()
#'           samples = simulator(n_sample=1, n_species=n_species, n_depth=configurations[[i]]$n_depth,ts_matrix=ts_matrix, n_timepoints=configurations[[i]]$n_timepoints, initial_state=inital_state, centers = centers)
#'           samples<-as.matrix(samples)
#'           sample_list <- list(n_species = n_species-1, n_depth=configurations[[i]]$n_depth, n_steps = configurations[[i]]$n_timepoints, n_clust = n_clust, y = samples, n_person = 1)
#'           result=cal_fit(sample_list)
#'           diff[k]=mean(result$matrix-ts_matrix)^2
#'         }
#'         Error<-mean(diff)
#'         statistics[i,j] = Error
#'       }
#'     }
#'     statistics
#'   }
#'}
