# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Canon without off target option
#' 
#' @export
Canon_sampling <- function(xin, yin, gRNAin, Gibbsnumberin, burninproportion, initial_betain, initial_gammain, pi_beta_shape_in, pi_beta_scale_in, a_betain, b_betain, sigma_dfin, verbosein) {
    .Call(`_Canon_Canon_sampling`, xin, yin, gRNAin, Gibbsnumberin, burninproportion, initial_betain, initial_gammain, pi_beta_shape_in, pi_beta_scale_in, a_betain, b_betain, sigma_dfin, verbosein)
}

#' Canon with off target effects
#' 
#' @export
Canon_sampling_off_target <- function(xin, yin, gRNAin, Gibbsnumberin, burninproportion, initial_betain, initial_gammain, pi_beta_shape_in, pi_beta_scale_in, a_betain, b_betain, pi_eta_shape_in, pi_eta_scale_in, a_etain, b_etain, sigma_dfin, verbosein) {
    .Call(`_Canon_Canon_sampling_off_target`, xin, yin, gRNAin, Gibbsnumberin, burninproportion, initial_betain, initial_gammain, pi_beta_shape_in, pi_beta_scale_in, a_betain, b_betain, pi_eta_shape_in, pi_eta_scale_in, a_etain, b_etain, sigma_dfin, verbosein)
}

