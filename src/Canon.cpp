// Code for conducting one sample MR with gRNA auto selection
// Parameters including beta_j/ gamma_j/ sigma_alpha/ sigma_beta/ alpha/ Sigma
// alpha has a N(0, 1) prior distribution

//library('Rcpp')
#include <RcppDist.h> // to generate the random number from beta distribution
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

//*******************************************************************//
//                MAIN FUNC                        //
//*******************************************************************//
//' Canon without off target option
//' 
//' @export
// [[Rcpp::export]]

SEXP Canon_sampling(SEXP xin, SEXP yin, SEXP gRNAin, SEXP Gibbsnumberin, 
               SEXP burninproportion, SEXP initial_betain, SEXP initial_gammain, 
               SEXP pi_beta_shape_in, SEXP pi_beta_scale_in, SEXP a_betain,
               SEXP b_betain, SEXP sigma_dfin, SEXP verbosein) {
  try {
    // Input initialization
    const int Gibbs_number = Rcpp::as<int>(Gibbsnumberin); // Total Gibbs sampling number
    const double burnin_p = Rcpp::as<double>(burninproportion); // Burnin proportion
    const double burnin_number = ceil(burnin_p * Gibbs_number); // Burnin sample number
    const double pi_beta_shape1 = Rcpp::as<double>(pi_beta_shape_in); // Pi_beta prior shape
    const double pi_beta_shape2 = Rcpp::as<double>(pi_beta_scale_in); // Pi_beta prior scale
    const vec x_expression = as<arma::vec>(xin); // Exposure gene expression
    const vec y_expression = as<arma::vec>(yin); // Outcome gene expression
    const mat gRNA = as<arma::mat>(gRNAin); // The entire gRNA matrix before selection
    const vec initial_beta = as<arma::vec>(initial_betain); // Initial beta value
    const vec initial_gamma = as<arma::vec>(initial_gammain); // Initial gamma value
    const double a_beta = Rcpp::as<double>(a_betain); // Sigma_beta2 prior variance shape
    const double b_beta = Rcpp::as<double>(b_betain); // Sigma_beta2 prior variance scale
    const double Sigma_df = Rcpp::as<double>(sigma_dfin); // Covariance matrix prior df
    const bool verbose = Rcpp::as<bool>(verbosein); // whether to print out
    
    // Parameters initial value initialization
    
    // Beta related parameters
    int p = gRNA.n_cols;
    int n = gRNA.n_rows;
    double pi_beta = pi_beta_shape1/(pi_beta_shape1 + pi_beta_shape2);
    double sigma_beta2 = b_beta/(a_beta - 1.0);
    // alpha related parameters
    double alpha = 0;
    // Sigma related parameters
    mat Sigma = eye(2, 2);
    mat Sigma_scale = eye(2, 2);
    // int Sigma_df = 1;
    
    
    // Stored Sampling paramters and useful identities used during sampling
    vec sample_gamma = initial_gamma; // sampled gamma indicator vector during each sampling
    vec sample_beta = initial_beta; // sampled beta for selected gRNAs or 0 for not selected gRNAs during each sampling
    vec sample_alpha(Gibbs_number); // sampled alpha value for all sampling
    vec sample_pi_beta(Gibbs_number); /// sampled pi_beta value for all sampling
    vec sample_sigma_beta2(Gibbs_number); // sampled beta prior variance for all sampling
    cube sample_cov(2, 2, Gibbs_number); // sampled bivariate normal covariance matrix stored in this cube
    
    mat sample_inclusion_indicator(Gibbs_number, p); // sampled selection probability for all sampling
    
    // Precalculate some identities to speed up the calculation for each sampling
    rowvec gRNA_square_sum = sum((gRNA % gRNA), 0); // gRNA matrix square column sum
    umat indexpre(p - 1, p, fill::zeros); // Precreate index when extracting the columns for each snp
    vec index = linspace<vec>(0, p - 1, p); // linearly spaced vectors
    for (int i = 0; i < p; i++) {
      vec indexnew = index.elem(find(index != i));
      uvec index_minusi = conv_to<uvec>::from(indexnew);
      indexpre.col(i) = index_minusi;
    }
    
    // Sampling step
    for (int s_id = 0; s_id < Gibbs_number; s_id++) {
      // Precalculate some identities to speed up the calculation during each gRNA update
      mat Sigmainv = inv_sympd(Sigma); // Sigma inverse
      mat kernel_variance_beta(2, 1, fill::ones);
      kernel_variance_beta(1, 0) = alpha; // Kernel matrix for faster calculation in variance of beta
      vec z_beta_all = gRNA * sample_beta;
      vec sample_beta_copy = sample_beta; // Create copy so that we can use formulas to save time
      double kernel_variance_beta_total = as_scalar(kernel_variance_beta.t() * Sigmainv * 
        kernel_variance_beta); //Precalculate the kernel in the variance of beta
      
      // Looping across different SNPs
      for (int j = 0; j < p; j++) {
        // Update beta_j
        double pos_beta_var = 1.0/(as_scalar(kernel_variance_beta_total * gRNA_square_sum(j)) +
                              1.0/sigma_beta2); // Need to check if we need to use gamma here // can move this outside of the loop

          // Calculation for the posterior beta mean
        vec z_j_beta_j = gRNA.col(j) * sample_beta(j);
        vec z_beta = z_beta_all - z_j_beta_j;
        vec z_beta_alpha = alpha * z_beta;
        mat leftpart = join_rows(x_expression - z_beta, y_expression - z_beta_alpha);
        double pos_beta_mean = pos_beta_var * as_scalar(sum(((leftpart * Sigmainv * kernel_variance_beta) % 
                                         gRNA.col(j)), 0));
        
        // cout << pos_beta_var << endl << pos_beta_mean << endl;
        
        // Update the gamma_j - binary indicators to select the betas
        // double log_gamma_equal_1 = 0.5 * pos_beta_mean * pos_beta_mean/pos_beta_var +
        //                         0.5 * log(pos_beta_var) - 0.5 * log(sigma_beta2) +
        //                         log(pi_beta);
        // double log_gamma_equal_0 = log(1 - pi_beta);
        // double maxtwo = max(log_gamma_equal_1, log_gamma_equal_0);
        // double nume = exp(log_gamma_equal_1 - maxtwo);
        // double deno = exp(log_gamma_equal_1 + log_gamma_equal_0 - 2 * maxtwo);
        double gamma_equal_1 = exp(0.5 * pos_beta_mean * pos_beta_mean/pos_beta_var +
          0.5 * log(pos_beta_var) - 0.5 * log(sigma_beta2) +
          log(pi_beta));
        double gamma_equal_0 = exp(log(1.0 - pi_beta));
        // cout << gamma_equal_1 << endl << gamma_equal_0 << endl;
          // Probability of gamma being selected
        double prob_gamma = gamma_equal_1/(gamma_equal_1 + gamma_equal_0); 
        // cout << prob_gamma << endl;
        double rand_num = as_scalar(randu(1)); // random number to indicate selection
        if (rand_num <= prob_gamma) {
          sample_gamma(j) = 1;
          sample_beta(j) = as_scalar(randn(1) * sqrt(pos_beta_var) + pos_beta_mean);
          sample_inclusion_indicator(s_id, j) = 1;
        } else {
          sample_gamma(j) = 0;
          sample_beta(j) = 0;
          sample_inclusion_indicator(s_id, j) = 0;
        }
        z_beta_all += gRNA.col(j) * as_scalar(sample_beta(j) - sample_beta_copy(j));
      }
      // Updating other parameters shared across all the SNPs
      // Update alpha
        // Precalculate gRNA times sample_beta
        // wondering if we could just use z_beta_all here?
      vec z_beta_complete = gRNA * sample_beta;
      double post_alpha_variance = 1.0/(as_scalar((sum((z_beta_complete) % (z_beta_complete), 0) *
                                        Sigmainv(1, 1))) + 1.0);
      double post_alpha_mean = post_alpha_variance * (as_scalar((sum((x_expression - z_beta_complete) % 
                                (z_beta_complete), 0) * Sigmainv(0, 1))) + as_scalar((sum(y_expression % 
                                (z_beta_complete), 0) * Sigmainv(1, 1))));
      sample_alpha(s_id) = as_scalar(randn(1) * sqrt(post_alpha_variance) + post_alpha_mean);
      alpha = sample_alpha(s_id);
      
      // Update pi_beta
      double pi_beta_post_shape1 = sum(sample_gamma) + pi_beta_shape1;
      double pi_beta_post_shape2 = sum(1.0 - sample_gamma) + pi_beta_shape2;
      
      sample_pi_beta(s_id) = r_4beta(pi_beta_post_shape1, pi_beta_post_shape2, 0, 1);
      pi_beta = sample_pi_beta(s_id);
      
      // Update sigma_beta2
      double post_sigma_beta2_scale = 2/(as_scalar(sum(sample_gamma % sample_beta % 
                                         sample_beta) + 2 * b_beta));
      double post_sigma_beta2_shape = 0.5 * sum(sample_gamma) + a_beta;
      post_sigma_beta2_scale = abs(post_sigma_beta2_scale);
      
      sample_sigma_beta2(s_id) = 1.0/as_scalar(randg(1, 
                            distr_param(post_sigma_beta2_shape, post_sigma_beta2_scale)));
      sigma_beta2 = sample_sigma_beta2(s_id);
      
      // Update the covariance matrix
        // Calculate the kernel in the posterior scale matrix
      mat post_cov_kernel(2, 2, fill::zeros);
      vec x_zbeta = x_expression - z_beta_complete;
      vec y_zbetaalpha = y_expression - alpha * z_beta_complete;
      post_cov_kernel(0, 0) = as_scalar(sum(x_zbeta % x_zbeta, 0));
      post_cov_kernel(0, 1) =  as_scalar(sum(x_zbeta % y_zbetaalpha, 0));
      post_cov_kernel(1, 0) = post_cov_kernel(0, 1);
      post_cov_kernel(1, 1) = as_scalar(sum(y_zbetaalpha % y_zbetaalpha, 0));
      
      sample_cov.slice(s_id) = iwishrnd(post_cov_kernel + Sigma_scale, n + Sigma_df);
      Sigma = sample_cov.slice(s_id);
      // Sigma = eye(2, 2);
      // show the sampling procedure with updated number
      if (verbose) {
        if ((s_id + 1) % 100 == 0) {
          if (s_id < burnin_number) {
            cout << "\rBurnin iteration: " << s_id + 1 << 
              "/" << burnin_number << flush;
          } else {
            cout << "\rSampling iteration: " << 
              s_id - burnin_number + 1 << 
                "/" << Gibbs_number - burnin_number << flush;
          }
        }
      }
    }
    
    // Get the summary for the parameters after burn-in sample
    double alpha_estimate = mean(sample_alpha.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double alpha_sd = stddev(sample_alpha.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double sigma_beta2_estimate = mean(sample_sigma_beta2.subvec((burnin_number - 1), (Gibbs_number - 1)));
    // double sigma_beta2_sd = stddev(sample_sigma_beta2.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double pi_beta_estimate = mean(sample_pi_beta.subvec((burnin_number - 1), (Gibbs_number - 1)));
    // double pi_beta_sd = stddev(sample_pi_beta.subvec((burnin_number - 1), (Gibbs_number - 1)));
    mat cov_estimate = mean(sample_cov.slices((burnin_number - 1), (Gibbs_number - 1)), 2); 
    // dim = 2 is the element-wise mean across cube
    
    return List::create(Rcpp::Named("alpha") = alpha_estimate,
                        Rcpp::Named("alpha_sd") = alpha_sd, 
                        Rcpp::Named("sigma_beta2_estimate") = sigma_beta2_estimate,
                        Rcpp::Named("sigma_beta2_sample") = sample_sigma_beta2,
                        Rcpp::Named("pi_beta_estimate") = pi_beta_estimate,
                        Rcpp::Named("pi_beta_sample") = sample_pi_beta,
                        Rcpp::Named("Inclusion_indicator") = sample_inclusion_indicator,
                        Rcpp::Named("cov_estimate") = cov_estimate);
    //                           // Rcpp::Named("sample_alpha") = sample_alpha,
    // 						  // Rcpp::Named("sample_pi") = sample_pi,
    // 						  // Rcpp::Named("sample_pi_gamma_31") = sample_pi_gamma_31,
    // 						 //Rcpp::Named("sample_pi_gamma_21") = sample_pi_gamma_21,
    // 						  //Rcpp::Named("sample_pi_2") = sample_pi_2,
    // 						 //Rcpp::Named("sample_rho") = sample_rho,
    // 						 // Rcpp::Named("rho") = rho_estimate,
    // 						 Rcpp::Named("sigma2x") = sigma2y_estimate,
    // 						 Rcpp::Named("sigma2y") = sigma2z_estimate,
    // 						 Rcpp::Named("sigmabeta") = sigmaz1_estimate,
    // 						 Rcpp::Named("selectionProbability") = sample_inclusion_probability,
    // 						 Rcpp::Named("selectionIndicator") = sample_inclusion_indicator,
    // 						 // Rcpp::Named("sigmaeta") = sigma2gamma_1_estimate,
    //                              Rcpp::Named("sd") = alpha_sd);
    //                              //Rcpp::Named("sigmaX") = sigma2y,
    //                               //Rcpp::Named("sigmaY") = sigma2z,
    //                              //Rcpp::Named("sigmabeta") = sigma2beta,
    //                              //Rcpp::Named("loglik_seq") = loglik_out,
    //                              //Rcpp::Named("loglik") = loglik_max,
    //                              //Rcpp::Named("iteration") = Iteration-1);
	} catch (std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}


//*******************************************************************//
//                MAIN FUNC                        //
//*******************************************************************//
//' Canon with off target effects
//' 
//' @export
// [[Rcpp::export]]

SEXP Canon_sampling_off_target(SEXP xin, SEXP yin, SEXP gRNAin, SEXP Gibbsnumberin, 
               SEXP burninproportion, SEXP initial_betain, SEXP initial_gammain, 
               SEXP pi_beta_shape_in, SEXP pi_beta_scale_in, SEXP a_betain,
               SEXP b_betain, SEXP pi_eta_shape_in, SEXP pi_eta_scale_in, SEXP a_etain,
               SEXP b_etain, SEXP sigma_dfin, SEXP verbosein) {
  try {
    // Input initialization
    const int Gibbs_number = Rcpp::as<int>(Gibbsnumberin); // Total Gibbs sampling number
    const double burnin_p = Rcpp::as<double>(burninproportion); // Burnin proportion
    const double burnin_number = ceil(burnin_p * Gibbs_number); // Burnin sample number
    const double pi_beta_shape1 = Rcpp::as<double>(pi_beta_shape_in); // Pi_beta prior shape
    const double pi_beta_shape2 = Rcpp::as<double>(pi_beta_scale_in); // Pi_beta prior scale
    const double pi_eta_shape1 = Rcpp::as<double>(pi_eta_shape_in); // Pi_eta prior shape
    const double pi_eta_shape2 = Rcpp::as<double>(pi_eta_scale_in); // Pi_eta prior scale
    const vec x_expression = as<arma::vec>(xin); // Exposure gene expression
    const vec y_expression = as<arma::vec>(yin); // Outcome gene expression
    const mat gRNA = as<arma::mat>(gRNAin); // The entire gRNA matrix before selection
    const vec initial_beta = as<arma::vec>(initial_betain); // Initial beta value
    const vec initial_gamma = as<arma::vec>(initial_gammain); // Initial gamma value
    const double a_beta = Rcpp::as<double>(a_betain); // Sigma_beta2 prior variance shape
    const double b_beta = Rcpp::as<double>(b_betain); // Sigma_beta2 prior variance scale
    const double a_eta = Rcpp::as<double>(a_etain); // Sigma_eta2 prior variance shape
    const double b_eta = Rcpp::as<double>(b_etain); // Sigma_eta2 prior variance scale
    const double Sigma_df = Rcpp::as<double>(sigma_dfin); // Covariance matrix prior df
    const bool verbose = Rcpp::as<bool>(verbosein); // whether to print out
    
    // Parameters initial value initialization
    
    // Beta related parameters
    int p = gRNA.n_cols;
    int n = gRNA.n_rows;
    double pi_beta = pi_beta_shape1/(pi_beta_shape1 + pi_beta_shape2);
    double sigma_beta2 = b_beta/(a_beta - 1.0);
    // Eta related parameters
    double pi_eta = pi_eta_shape1/(pi_eta_shape1 + pi_eta_shape2);
    double sigma_eta2 = b_eta/(a_eta - 1.0);
    // alpha related parameters
    double alpha = 0;
    // Sigma related parameters
    mat Sigma = eye(2, 2);
    mat Sigma_scale = eye(2, 2);
    // int Sigma_df = 1;
    
    
    // Stored Sampling paramters and useful identities used during sampling
    vec sample_gamma = initial_gamma; // sampled gamma indicator vector during each sampling
    vec sample_tau = initial_gamma; // sampled tau indicator for eta vector during each sampling
    vec sample_beta = initial_beta; // sampled beta for selected gRNAs or 0 for not selected gRNAs during each sampling
    vec sample_eta = zeros(p); // sampled eta for selected gRNAs or 0 for not selected gRNAs during each sampling
    vec sample_alpha(Gibbs_number); // sampled alpha value for all sampling
    vec sample_pi_beta(Gibbs_number); /// sampled pi_beta value for all sampling
    vec sample_pi_eta(Gibbs_number); /// sampled pi_eta value for all sampling
    vec sample_sigma_beta2(Gibbs_number); // sampled beta prior variance for all sampling
    vec sample_sigma_eta2(Gibbs_number); // sampled eta prior variance for all sampling
    cube sample_cov(2, 2, Gibbs_number); // sampled bivariate normal covariance matrix stored in this cube
    
    mat sample_inclusion_indicator(Gibbs_number, p); // sampled selection probability for all sampling
    
    // Precalculate some identities to speed up the calculation for each sampling
    rowvec gRNA_square_sum = sum((gRNA % gRNA), 0); // gRNA matrix square column sum
    umat indexpre(p - 1, p, fill::zeros); // Precreate index when extracting the columns for each snp
    vec index = linspace<vec>(0, p - 1, p); // linearly spaced vectors
    for (int i = 0; i < p; i++) {
      vec indexnew = index.elem(find(index != i));
      uvec index_minusi = conv_to<uvec>::from(indexnew);
      indexpre.col(i) = index_minusi;
    }
    
    // Sampling step
    for (int s_id = 0; s_id < Gibbs_number; s_id++) {
      // Precalculate some identities to speed up the calculation during each gRNA update
      mat Sigmainv = inv_sympd(Sigma); // Sigma inverse
      mat kernel_variance_beta(2, 1, fill::ones);
      kernel_variance_beta(1, 0) = alpha; // Kernel matrix for faster calculation in variance of beta
      mat kernel_variance_eta(2, 1, fill::zeros);
      kernel_variance_eta(1, 0) = 1.0; // Kernel matrix for faster calculation in variance of eta
      vec z_beta_all = gRNA * sample_beta; // Z * beta for all gRNAs
      vec sample_beta_copy = sample_beta; // Create copy so that we can use formulas to save time
      vec z_eta_all = gRNA * sample_eta; // Z * eta for all gRNAs
      vec sample_eta_copy = sample_eta; // Create copy so that we can use formulas to save time
      double kernel_variance_beta_total = as_scalar(kernel_variance_beta.t() * Sigmainv * 
        kernel_variance_beta); //Precalculate the kernel in the variance of beta
      
      // Looping across different SNPs
      for (int j = 0; j < p; j++) {
        // Update beta_j
        double pos_beta_var = 1.0/(as_scalar(kernel_variance_beta_total * gRNA_square_sum(j)) +
                              1.0/sigma_beta2); // Need to check if we need to use gamma here

          // Calculation for the posterior beta mean
        vec z_j_beta_j = gRNA.col(j) * sample_beta(j); 
        vec z_beta = z_beta_all - z_j_beta_j; // Z_{-j} * beta_{-j}
        vec z_beta_alpha = alpha * z_beta;
        mat leftpart = join_rows(x_expression - z_beta, y_expression - z_beta_alpha - z_eta_all);
        double pos_beta_mean = pos_beta_var * as_scalar(sum(((leftpart * Sigmainv * kernel_variance_beta) % 
                                         gRNA.col(j)), 0));
        
        // cout << pos_beta_var << endl << pos_beta_mean << endl;
        
        // Update the gamma_j - binary indicators to select the betas
        // double log_gamma_equal_1 = 0.5 * pos_beta_mean * pos_beta_mean/pos_beta_var +
        //                         0.5 * log(pos_beta_var) - 0.5 * log(sigma_beta2) +
        //                         log(pi_beta);
        // double log_gamma_equal_0 = log(1 - pi_beta);
        // double maxtwo = max(log_gamma_equal_1, log_gamma_equal_0);
        // double nume = exp(log_gamma_equal_1 - maxtwo);
        // double deno = exp(log_gamma_equal_1 + log_gamma_equal_0 - 2 * maxtwo);
        double gamma_equal_1 = exp(0.5 * pos_beta_mean * pos_beta_mean/pos_beta_var +
          0.5 * log(pos_beta_var) - 0.5 * log(sigma_beta2) +
          log(pi_beta) + sample_tau(j) * log(pi_eta) + (1 - sample_tau(j)) * log(1 - pi_eta));
        double gamma_equal_0 = exp(log(1.0 - pi_beta));
        // cout << gamma_equal_1 << endl << gamma_equal_0 << endl;
          // Probability of gamma being selected
        double prob_gamma = gamma_equal_1/(gamma_equal_1 + gamma_equal_0); 
        // cout << prob_gamma << endl;
        double rand_num = as_scalar(randu(1)); // random number to indicate selection
        if (rand_num <= prob_gamma) {
          sample_gamma(j) = 1;
          sample_beta(j) = as_scalar(randn(1) * sqrt(pos_beta_var) + pos_beta_mean);
          sample_inclusion_indicator(s_id, j) = 1;
        } else {
          sample_gamma(j) = 0;
          sample_beta(j) = 0;
          sample_inclusion_indicator(s_id, j) = 0;
        }
        z_beta_all += gRNA.col(j) * as_scalar(sample_beta(j) - sample_beta_copy(j)); // update Z * beta for all gRNAs

        // Update eta_j only if beta_j != 0
        if (sample_gamma(j) == 1) {
          double pos_eta_var = 1.0/(as_scalar(Sigmainv(1, 1) * gRNA_square_sum(j)) +
                              1.0/sigma_eta2); // Need to check if we need to use gamma here

          // Calculation for the posterior eta mean
          vec z_j_eta_j = gRNA.col(j) * sample_eta(j);
          vec z_eta = z_eta_all - z_j_eta_j;
          mat leftpart_eta = join_rows(x_expression - z_beta_all, y_expression - z_beta_all * alpha - z_eta);
          double pos_eta_mean = pos_eta_var * as_scalar(sum(((leftpart_eta * Sigmainv * kernel_variance_eta) % 
                                           gRNA.col(j)), 0));

          double eta_equal_1 = exp(0.5 * pos_eta_mean * pos_eta_mean/pos_eta_var +
            0.5 * log(pos_eta_var) - 0.5 * log(sigma_eta2) +
            log(pi_eta));
          double eta_equal_0 = exp(log(1.0 - pi_eta));
          // cout << gamma_equal_1 << endl << gamma_equal_0 << endl;
            // Probability of gamma being selected
          double prob_eta = eta_equal_1/(eta_equal_1 + eta_equal_0); 
          // cout << prob_gamma << endl;
          double rand_num_eta = as_scalar(randu(1)); // random number to indicate selection
          if (rand_num_eta <= prob_eta) {
            sample_tau(j) = 1;
            sample_eta(j) = as_scalar(randn(1) * sqrt(pos_eta_var) + pos_eta_mean);
          } else {
            sample_tau(j) = 0;
            sample_eta(j) = 0;
          }
        } else {
            sample_tau(j) = 0;
            sample_eta(j) = 0;
        }
        
        z_eta_all += gRNA.col(j) * as_scalar(sample_eta(j) - sample_eta_copy(j)); // update Z * eta for all gRNAs
      }
      // Updating other parameters shared across all the SNPs
      // Update alpha
        // Precalculate gRNA times sample_beta
        // wondering if we could just use z_beta_all here?
      vec z_beta_complete = gRNA * sample_beta;
      vec z_eta_complete = gRNA * sample_eta;
      double post_alpha_variance = 1.0/(as_scalar((sum((z_beta_complete) % (z_beta_complete), 0) *
                                        Sigmainv(1, 1))) + 1.0);
      double post_alpha_mean = post_alpha_variance * (as_scalar((sum((x_expression - z_beta_complete) % 
                                (z_beta_complete), 0) * Sigmainv(0, 1))) + as_scalar((sum((y_expression - z_eta_complete) % 
                                (z_beta_complete), 0) * Sigmainv(1, 1))));
      sample_alpha(s_id) = as_scalar(randn(1) * sqrt(post_alpha_variance) + post_alpha_mean);
      alpha = sample_alpha(s_id);
      
      // Update pi_beta
      double pi_beta_post_shape1 = sum(sample_gamma) + pi_beta_shape1;
      double pi_beta_post_shape2 = sum(1.0 - sample_gamma) + pi_beta_shape2;
      
      sample_pi_beta(s_id) = r_4beta(pi_beta_post_shape1, pi_beta_post_shape2, 0, 1);
      pi_beta = sample_pi_beta(s_id);

      // Update pi_eta
      double pi_eta_post_shape1 = sum(sample_gamma % sample_tau) + pi_eta_shape1;
      double pi_eta_post_shape2 = sum((1.0 - sample_tau) % sample_gamma) + pi_eta_shape2;
      
      sample_pi_eta(s_id) = r_4beta(pi_eta_post_shape1, pi_eta_post_shape2, 0, 1);
      pi_eta = sample_pi_eta(s_id);
      
      // Update sigma_beta2
      double post_sigma_beta2_scale = 2/(as_scalar(sum(sample_gamma % sample_beta % 
                                         sample_beta) + 2 * b_beta));
      double post_sigma_beta2_shape = 0.5 * sum(sample_gamma) + a_beta;
      post_sigma_beta2_scale = abs(post_sigma_beta2_scale);
      
      sample_sigma_beta2(s_id) = 1.0/as_scalar(randg(1, 
                            distr_param(post_sigma_beta2_shape, post_sigma_beta2_scale)));
      sigma_beta2 = sample_sigma_beta2(s_id);

       // Update sigma_eta2
      double post_sigma_eta2_scale = 2/(as_scalar(sum(sample_tau % sample_eta % 
                                         sample_eta) + 2 * b_eta));
      double post_sigma_eta2_shape = 0.5 * sum(sample_tau) + a_eta;
      post_sigma_eta2_scale = abs(post_sigma_eta2_scale);
      
      sample_sigma_eta2(s_id) = 1.0/as_scalar(randg(1, 
                            distr_param(post_sigma_eta2_shape, post_sigma_eta2_scale)));
      sigma_eta2 = sample_sigma_eta2(s_id);
      
      // Update the covariance matrix
        // Calculate the kernel in the posterior scale matrix
      mat post_cov_kernel(2, 2, fill::zeros);
      vec x_zbeta = x_expression - z_beta_complete;
      vec y_zbetaalpha_zeta = y_expression - alpha * z_beta_complete - z_eta_complete;
      post_cov_kernel(0, 0) = as_scalar(sum(x_zbeta % x_zbeta, 0));
      post_cov_kernel(0, 1) =  as_scalar(sum(x_zbeta % y_zbetaalpha_zeta, 0));
      post_cov_kernel(1, 0) = post_cov_kernel(0, 1);
      post_cov_kernel(1, 1) = as_scalar(sum(y_zbetaalpha_zeta % y_zbetaalpha_zeta, 0));
      
      sample_cov.slice(s_id) = iwishrnd(post_cov_kernel + Sigma_scale, n + Sigma_df);
      Sigma = sample_cov.slice(s_id);
      // Sigma = eye(2, 2);
      // show the sampling procedure with updated number
      if (verbose) {
        if ((s_id + 1) % 100 == 0) {
          if (s_id < burnin_number) {
            cout << "\rBurnin iteration: " << s_id + 1 << 
              "/" << burnin_number << flush;
          } else {
            cout << "\rSampling iteration: " << 
              s_id - burnin_number + 1 << 
                "/" << Gibbs_number - burnin_number << flush;
          }
        }
      }
    }
    
    // Get the summary for the parameters after burn-in sample
    double alpha_estimate = mean(sample_alpha.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double alpha_sd = stddev(sample_alpha.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double sigma_beta2_estimate = mean(sample_sigma_beta2.subvec((burnin_number - 1), (Gibbs_number - 1)));
    // double sigma_beta2_sd = stddev(sample_sigma_beta2.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double pi_beta_estimate = mean(sample_pi_beta.subvec((burnin_number - 1), (Gibbs_number - 1)));
    // double pi_beta_sd = stddev(sample_pi_beta.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double sigma_eta2_estimate = mean(sample_sigma_eta2.subvec((burnin_number - 1), (Gibbs_number - 1)));
    double pi_eta_estimate = mean(sample_pi_eta.subvec((burnin_number - 1), (Gibbs_number - 1)));
    mat cov_estimate = mean(sample_cov.slices((burnin_number - 1), (Gibbs_number - 1)), 2); 
    // dim = 2 is the element-wise mean across cube
    
    return List::create(Rcpp::Named("alpha") = alpha_estimate,
                        Rcpp::Named("alpha_sd") = alpha_sd, 
                        Rcpp::Named("sigma_beta2_estimate") = sigma_beta2_estimate,
                        Rcpp::Named("sigma_beta2_sample") = sample_sigma_beta2,
                        Rcpp::Named("pi_beta_estimate") = pi_beta_estimate,
                        Rcpp::Named("pi_beta_sample") = sample_pi_beta,
                        Rcpp::Named("sigma_eta2_estimate") = sigma_eta2_estimate,
                        Rcpp::Named("sigma_eta2_sample") = sample_sigma_eta2,
                        Rcpp::Named("pi_eta_estimate") = pi_eta_estimate,
                        Rcpp::Named("pi_eta_sample") = sample_pi_eta,
                        Rcpp::Named("Inclusion_indicator") = sample_inclusion_indicator,
                        Rcpp::Named("cov_estimate") = cov_estimate);
    //                           // Rcpp::Named("sample_alpha") = sample_alpha,
    //              // Rcpp::Named("sample_pi") = sample_pi,
    //              // Rcpp::Named("sample_pi_gamma_31") = sample_pi_gamma_31,
    //             //Rcpp::Named("sample_pi_gamma_21") = sample_pi_gamma_21,
    //              //Rcpp::Named("sample_pi_2") = sample_pi_2,
    //             //Rcpp::Named("sample_rho") = sample_rho,
    //             // Rcpp::Named("rho") = rho_estimate,
    //             Rcpp::Named("sigma2x") = sigma2y_estimate,
    //             Rcpp::Named("sigma2y") = sigma2z_estimate,
    //             Rcpp::Named("sigmabeta") = sigmaz1_estimate,
    //             Rcpp::Named("selectionProbability") = sample_inclusion_probability,
    //             Rcpp::Named("selectionIndicator") = sample_inclusion_indicator,
    //             // Rcpp::Named("sigmaeta") = sigma2gamma_1_estimate,
    //                              Rcpp::Named("sd") = alpha_sd);
    //                              //Rcpp::Named("sigmaX") = sigma2y,
    //                               //Rcpp::Named("sigmaY") = sigma2z,
    //                              //Rcpp::Named("sigmabeta") = sigma2beta,
    //                              //Rcpp::Named("loglik_seq") = loglik_out,
    //                              //Rcpp::Named("loglik") = loglik_max,
    //                              //Rcpp::Named("iteration") = Iteration-1);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
  }
  return R_NilValue;
}
