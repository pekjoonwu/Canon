#' @title Casual relationship identification using one sample instrumental variable model
#' @description Canon is able to model an initial set of candidate gRNAs
#' and perform automated instrument selection to identify suitable gRNAs to serve as instrumental variables. 
#' Canon relies on a scalable sampling-based inference algorithm to identify genes that are potentially causally influenced by perturbed target genes across diverse sc-CRISPR platforms.
#' @param x_expression the gene expression for the exposure gene
#' @param y_expression the gene expression for the outcome gene
#' @param gRNA the gRNA matrix
#' @param Gibbsnumber the number of Gibbs sampling iterations with the default to be 3000
#' @param burninproportion  the proportion to burn in from Gibbs sampling iterations, with default to be 20%% 
#' @param pi_beta_shape the prior shape parameter for pi_beta with the default to be 0.7
#' @param pi_beta_scale the prior scale parameter for pi_beta with the default to be 0.7
#' @param a_beta the prior parameter for sigma_beta2 with the default to be 21
#' @param b_beta the prior parameter for sigma_beta2 with the default to be 0.01
#' @param pi_eta_shape the prior shape parameter for pi_eta with the default to be 1
#' @param pi_eta_scale the prior scale parameter for pi_eta with the default to be 1
#' @param a_eta the prior parameter for sigma_eta2 with the default to be 51
#' @param b_eta the prior parameter for sigma_eta2 with the default to be 0.1
#' @param initial_beta the initial value for the gRNA effect sizes
#' @param initial_gamma the initial indicator for whether the gRNA is a valid instrument or not
#' @param sigma_df the degree of freedom in the inverse wishart prior for covariance matrix
#' @param off_target whether the algorithm is using a model with off target effects
#' @param verbose whether to print out the fitting process
#' @return A list of estimated parameters where the main parameters are the causal effect estimate (causal_effect) and corresponding p-value (causal_pvalue)  

run_Canon <- function(x_expression = NULL, y_expression = NULL, 
                                    gRNA = NULL, Gibbsnumber = 3000, burninproportion = 0.2, 
                                    pi_beta_shape = 0.7, pi_beta_scale = 0.7,
                                    a_beta = 21, b_beta = 0.01,
                                    pi_eta_shape = 1, pi_eta_scale = 1,
                                    a_eta = 51, b_eta = 0.1,
                                    initial_beta = rep(0, ncol(gRNA)), initial_gamma = rep(0, ncol(gRNA)),
                                    sigma_df = 2, off_target = T,
                                    verbose = T) {
  ## Preprocessing so that the the mean is 0 and variance is 1 for expression and gRNA status
  x_expression <- scale(x_expression)
  y_expression <- scale(y_expression)
  gRNA <- scale(gRNA)

  ## Run the main sampling function and output the results
  ## initial_betain <- rep(0, ncol(gRNA))
  if (off_target == T) {
    results <- Canon_sampling_off_target(xin = x_expression, yin = y_expression, gRNAin = gRNA, 
                       Gibbsnumberin = Gibbsnumber, burninproportion = burninproportion,
                       initial_betain = initial_beta, initial_gammain = initial_gamma,
                       pi_beta_shape_in = pi_beta_shape, pi_beta_scale_in = pi_beta_scale,
                       a_betain = a_beta, b_betain = b_beta, 
                       pi_eta_shape_in = pi_eta_shape, pi_eta_scale_in = pi_eta_scale,
                       a_etain = a_eta, b_etain = b_eta, 
                       sigma_dfin = sigma_df,
                       verbosein = verbose)
  } else {
    results <- Canon_sampling(xin = x_expression, yin = y_expression, gRNAin = gRNA, 
                       Gibbsnumberin = Gibbsnumber, burninproportion = burninproportion,
                       initial_betain = initial_beta, initial_gammain = initial_gamma,
                       pi_beta_shape_in = pi_beta_shape, pi_beta_scale_in = pi_beta_scale,
                       a_betain = a_beta, b_betain = b_beta, 
                       sigma_dfin = sigma_df,
                       verbosein = verbose)
  }

  ## pvalue <- 2 * (1 - pnorm(abs(results$alpha/results$alpha_sd)))
  ## pvalue <- 2 * (pnorm(-abs(results$alpha/results$alpha_sd)))
  ## update the pvalue calculation by incorporating the prior distribution on alpha
  pvalue <- 2 * (pnorm(-abs(results$alpha/(results$alpha_sd * sqrt(1 - results$alpha_sd ^ 2)))))
  result <- list()
  result$causal_effect <- results$alpha
  result$causal_effect_sd <- results$alpha_sd
  result$causal_pvalue <- pvalue
  result$sigma_beta2 <- results$sigma_beta2_estimate
  result$sigma_beta2_samples <- results$sigma_beta2_sample
  result$pi_beta <- results$pi_beta_estimate
  result$pi_beta_samples <- results$pi_beta_sample
  result$Sigma <- results$cov_estimate
  result$gamma_ind <- results$Inclusion_indicator
  if (off_target == T) {
    result$sigma_eta2 <- results$sigma_eta2_estimate
    result$sigma_eta2_samples <- results$sigma_eta2_sample
    result$pi_eta <- results$pi_eta_estimate
    result$pi_eta_samples <- results$pi_eta_sample
  }
  return(result)
}






