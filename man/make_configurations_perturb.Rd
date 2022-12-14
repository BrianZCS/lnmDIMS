% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_test.R
\name{make_configurations_perturb}
\alias{make_configurations_perturb}
\title{configurations generator for perturbation}
\usage{
make_configurations_perturb(
  ns,
  n_depth,
  alpha,
  intervals = NULL,
  n_timepoints = NULL
)
}
\arguments{
\item{ns}{the number of samples}

\item{n_depth}{sequence depth}

\item{alpha}{the maximum perturbation intensity (between 0 and 1)}

\item{n_timepoints}{number of time points}

\item{interval}{the interval of perturbation (interval\link{1}-interval\link{2}): increase the permutation intensity.(interval\link{2}-interval\link{3}): achieve maximized perturbation intensity. (interval\link{3}-interval\link{4}): decrease the permutation intensity.}
}
\value{
configurations the configurations used in the power test
}
\description{
The function generate all configurations used in the power test
}
