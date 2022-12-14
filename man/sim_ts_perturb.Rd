% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_without_data.R
\name{sim_ts_perturb}
\alias{sim_ts_perturb}
\title{Simulate microbiome perturbation data}
\usage{
sim_ts_perturb(
  endtime,
  n_species,
  n_depth,
  Sigma = NULL,
  beta,
  t1,
  t2,
  time,
  species,
  duration
)
}
\arguments{
\item{endtime}{the number of time points of the sample}

\item{n_species}{number of species}

\item{n_depth}{sequence depth}

\item{Sigma}{correlation matrix for muti-nomal distribution}

\item{beta}{effect of perturbation}

\item{t1}{time of reaching the maximum of perturbation}

\item{t2}{time of ending the maximum of perturbation}

\item{time}{the start time for perturbation}

\item{species}{perturbed species}

\item{duration}{the duration of perturbation}
}
\value{
the simulated data
}
\description{
To simulate microbiome longitudinal count data for perturbation experiments with user defined starting and ending time point, effect size and duration of perturbation
}
\examples{
sim_ts_perturb(endtime=30, n_species=30, 1000, Sigma=NULL, beta=1,species=c(1,2,3),time=4, t1=8,t2=12,duration = 16)
}
