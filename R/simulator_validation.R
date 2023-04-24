#' Calculate relative abundance
#' 
#' Calculate the relative abundance of the microbiome data
#' @param matrix microbiome abundance data
#' @return return the relative abundance
#' @export
relative_abundence = function(matrix){
  for(i in 1:nrow(matrix)){
    matrix[i,] = matrix[i,]/sum(matrix[i,])
  }
  return(matrix)
}

#' Compare simulated with real data
#' 
#' Violin plot to compare the mean and variance of simulated and real data.
#' @param sim_data simulated data
#' @param obs_datas observed data
#' @param relative_abundance compare the relative_abundance or not
#' @export
violin_plot<-function(sim_data, obs_data, relative_abundance){
  if (relative_abudance == TRUE){
    sim_data = relative_abundance(sim_data)
    obs_data = relative_abundance(obs_data)
  }
  mean_sim <- colMeans(sim_data)
  mean_sim_data<-data.frame("dataset" = rep("sim", ncol(sim_data)), "measure" = rep("mean", ncol(sim_data)), "value" = mean_sim)
  
  var_sim <- apply(sim_data,2,var)
  var_sim_data<-data.frame("dataset" = rep("sim", ncol(sim_data)), "measure" = rep("var", ncol(sim_data)), "value" = var_sim)
  
  mean_obs <- colMeans(obs_data)
  mean_obs_data<-data.frame("dataset" = rep("obs", ncol(obs_data)), "measure" = rep("mean", ncol(obs_data)), "value" = mean_obs)
  
  var_obs <- unname(apply(obs_data,2,var))
  var_obs_data <- data.frame("dataset" = rep("obs", ncol(obs_data)), "measure" = rep("var", ncol(obs_data)), "value" = mean_obs)
  
  stat<-bind_rows(mean_sim_data,var_sim_data, mean_obs_data,var_obs_data)
  
  ggplot(stat, aes(x = dataset, y = value)) +
    scale_y_log10()+
    geom_violin(scale = 'width', trim = TRUE) +
    facet_wrap(~measure, scales = "free", ncol = 3)
}

#' Compare simulated with real data through heatmap
#' 
#' heatmap to compare the structure and abundance of simulated and real data.
#' @param sim_data simulated data
#' @param obs_datas observed data
#' @param relative_abundance compare the relative_abundance or not
#' @param log_transform log transfrom or not
#' @export
heat_map<-function(sim_data, obs_data, relative_abundance, log_transform){
  if (relative_abudance == TRUE){
    sim_data = relative_abundance(sim_data)
    obs_data = relative_abundance(obs_data)
  }
  if(log_transform == TRUE){
    sim_data = log(1+sim_data)
    obs_data = log(1+obs_data)
  }
  par(mfrow=c(1,2))
  heatmap(t(sim_data))
  heatmap(t(obs_data))
  par(mfrow=c(1,1))
}

#' Compare simulated with real data
#' 
#' species-wise scatter plot to compare the abundance of simulated and real data.
#' @param sim_data simulated data
#' @param obs_datas observed data
#' @param relative_abundance compare the relative_abundance or not
#' @export
specieswise_plot<-function(sim_data,obs_data,relative_abundance){
  if (relative_abudance == TRUE){
    sim_data = relative_abundance(sim_data)
    obs_data = relative_abundance(obs_data)
  }
  sim.vs.real<-data.frame("species"=rep(1:ncol(sim_data), nrow(sim_data)), "sim"=c(log(1+t(sim_data))), "real"=c(log(1+t(obs_data))))
  ggplot(sim.vs.real)+
    geom_point(aes(sim_data, obs_data), alpha = 0.3)+
    facet_wrap(~species, scales = "free_y")+
    geom_abline()
}
