% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Canon.R
\name{run_Canon}
\alias{run_Canon}
\title{Casual relationship identification using one sample instrumental variable model}
\usage{
run_Canon(
  x_expression = NULL,
  y_expression = NULL,
  gRNA = NULL,
  Gibbsnumber = 3000,
  burninproportion = 0.2,
  pi_beta_shape = 0.7,
  pi_beta_scale = 0.7,
  a_beta = 21,
  b_beta = 0.01,
  pi_eta_shape = 1,
  pi_eta_scale = 1,
  a_eta = 51,
  b_eta = 0.1,
  initial_beta = rep(0, ncol(gRNA)),
  initial_gamma = rep(0, ncol(gRNA)),
  sigma_df = 2,
  off_target = T,
  verbose = T
)
}
\arguments{
\item{x_expression}{the gene expression for the exposure gene}

\item{y_expression}{the gene expression for the outcome gene}

\item{gRNA}{the gRNA matrix}

\item{Gibbsnumber}{the number of Gibbs sampling iterations with the default to be 3000}

\item{burninproportion}{the proportion to burn in from Gibbs sampling iterations, with default to be 20%}

\item{pi_beta_shape}{the prior shape parameter for pi_beta with the default to be 0.7}

\item{pi_beta_scale}{the prior scale parameter for pi_beta with the default to be 0.7}

\item{a_beta}{the prior parameter for sigma_beta2 with the default to be 21}

\item{b_beta}{the prior parameter for sigma_beta2 with the default to be 0.01}

\item{pi_eta_shape}{the prior shape parameter for pi_eta with the default to be 1}

\item{pi_eta_scale}{the prior scale parameter for pi_eta with the default to be 1}

\item{a_eta}{the prior parameter for sigma_eta2 with the default to be 51}

\item{b_eta}{the prior parameter for sigma_eta2 with the default to be 0.1}

\item{initial_beta}{the initial value for the gRNA effect sizes}

\item{initial_gamma}{the initial indicator for whether the gRNA is a valid instrument or not}

\item{sigma_df}{the degree of freedom in the inverse wishart prior for covariance matrix}

\item{off_target}{whether the algorithm is using a model with off target effects}

\item{verbose}{whether to print out the fitting process}
}
\value{
A list of estimated parameters including the p values for the causal effect test
}
\description{
Canon is able to model an initial set of candidate gRNAs
and perform automated instrument selection to identify suitable gRNAs to serve as instrumental variables. 
Canon relies on a scalable sampling-based inference algorithm to identify genes that are potentially causally influenced by perturbed target genes across diverse sc-CRISPR platforms.
}
