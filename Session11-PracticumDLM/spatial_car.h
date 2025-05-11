#ifndef SPATIAL_CAR_H
#define SPATIAL_CAR_H

#include <RcppArmadillo.h> // Include necessary headers within the header guard
#include <cmath>           // For std::log, std::exp, std::isfinite

// Use inline keyword to prevent multiple definition errors if this header
// is included in multiple source files (though less likely in typical Rcpp packages)

//--------------------------------------------------------------------------
// Helper function for the log-sum-exp trick
//--------------------------------------------------------------------------
inline double log_sum_exp(const arma::vec& x) {
    if (x.is_empty()) {
        // Return a value indicating error or empty input, e.g., NaN or -Inf
        // Rcpp::stop might be too harsh for a helper in a header.
        return -arma::datum::inf;
    }
    double max_val = x.max();
    if (!std::isfinite(max_val)) {
         // Handle cases where max is Inf or -Inf or NaN
         if (max_val == arma::datum::inf) return arma::datum::inf;
         // If max is -Inf, all elements are -Inf, result is -Inf
         // If max is NaN, result is NaN
         return max_val;
    }
    double sum_exp = arma::accu(arma::exp(x - max_val));
     if (sum_exp <= 0) { // Handle underflow or non-positive values
         return -arma::datum::inf;
     }
    return max_val + std::log(sum_exp);
}


//--------------------------------------------------------------------------
// Function to sample spatial correlation parameter for CAR model
//--------------------------------------------------------------------------
inline double MH_spatial_correlation_CAR(
    const arma::sp_mat& W,  // Adjacency matrix (p x p)
    const arma::mat& Z,     // Data matrix (p x T) - Assumed centered
    double rho1,            // Conditional variance factor
    double min_tau,         // Min value for spatial parameter (rho2)
    double max_tau          // Max value for spatial parameter (rho2)
) {

    // Input checks (optional in header, can be done in calling function or kept here)
    if (W.n_rows != W.n_cols) {
        Rcpp::stop("MH_spatial_correlation_CAR_cpp: W must be a square matrix.");
    }
    if (W.n_rows != Z.n_rows) {
        Rcpp::stop("MH_spatial_correlation_CAR_cpp: Number of rows in W must match number of rows in Z.");
    }

    arma::uword N = Z.n_rows; // Number of spatial locations (p)
    arma::uword T = Z.n_cols; // Number of time points or samples

    // Define the grid for the spatial parameter rho2
    arma::vec allowed_range = arma::linspace<arma::vec>(min_tau, max_tau, 15);
    arma::uword nr = allowed_range.n_elem;
    arma::vec lpos(nr, arma::fill::zeros);

    arma::sp_mat D(N, N);
	D.diag() += arma::sum(W, 1);

    for (arma::uword i = 0; i < nr; ++i) {
        double range = allowed_range(i);
        arma::sp_mat invS_unscaled = D - range * W;
        invS_unscaled = 0.5 * (invS_unscaled + invS_unscaled.t());
        arma::sp_mat invS = invS_unscaled * (1.0 / rho1);
        arma::mat invS_dense = arma::mat(invS);

         if (!invS_dense.is_sympd(1e-8)) { // Use tolerance for sympd check
             lpos(i) = -arma::datum::inf;
             continue;
         }

        double appo = 0;
        for (arma::uword t = 0; t < T; ++t) {
            appo += arma::as_scalar(Z.col(t).t() * invS_dense * Z.col(t));
        }

        arma::vec eigval;
        bool eig_success = arma::eig_sym(eigval, invS_dense);

        if (!eig_success) {
            Rcpp::warning("MH_spatial_correlation_CAR_cpp: Eigenvalue decomposition failed for range = %.4f.", range);
            lpos(i) = -arma::datum::inf;
            continue;
        }

        arma::uvec positive_idx = arma::find(eigval > 1e-12);
        if (positive_idx.n_elem < eigval.n_elem) {
             lpos(i) = -arma::datum::inf;
             continue;
        }
        if (positive_idx.n_elem == 0) {
             lpos(i) = -arma::datum::inf;
             continue;
        }

        double log_det_invS = arma::sum(arma::log(eigval(positive_idx)));
        if (!std::isfinite(log_det_invS)) {
             lpos(i) = -arma::datum::inf;
             continue;
        }

        lpos(i) = 0.5 * static_cast<double>(T) * log_det_invS - 0.5 * appo;
        if (!std::isfinite(lpos(i))) {
              lpos(i) = -arma::datum::inf;
        }
    }

    if (arma::all(lpos == -arma::datum::inf)) {
        Rcpp::warning("MH_spatial_correlation_CAR_cpp: All log-posterior values are -Inf. Returning middle of range.");
        return (min_tau + max_tau) / 2.0;
    }

    double log_normalizer = log_sum_exp(lpos);
    arma::vec log_probx = lpos - log_normalizer;
    arma::vec probx = arma::exp(log_probx);
    probx.elem(arma::find_nonfinite(probx)).zeros();

    double sum_probx = arma::accu(probx);
    if (sum_probx <= 1e-12 || !std::isfinite(sum_probx)) {
        Rcpp::warning("MH_spatial_correlation_CAR_cpp: Probabilities sum close to zero or non-finite (sum=%.2e). Using uniform weights.", sum_probx);
        probx.fill(1.0 / static_cast<double>(nr));
    } else {
        probx /= sum_probx;
    }

    Rcpp::NumericVector range_vec = Rcpp::wrap(allowed_range);
    Rcpp::NumericVector prob_vec = Rcpp::wrap(probx);

    if (Rcpp::is_true(Rcpp::any(Rcpp::is_na(prob_vec)))) {
         Rcpp::warning("MH_spatial_correlation_CAR_cpp: NA values found in sampling probabilities. Using uniform weights.");
         prob_vec = Rcpp::NumericVector(nr, 1.0/nr);
    }

    Rcpp::NumericVector sampled_value = Rcpp::sample(range_vec, 1, false, prob_vec);
    return sampled_value[0];
}


//--------------------------------------------------------------------------
// Function to sample conditional variance parameter rho1 for CAR model
// Assumes rho1 ~ IG(a1, b1) prior => 1/rho1 ~ Gamma(a1, b1) prior
//--------------------------------------------------------------------------
inline double posterior_conditional_variance(
    const arma::mat& Z,         // Data matrix (p x T, e.g., B or B_diff)
    const arma::sp_mat& invCorrZ, // Spatial precision structure (p x p, e.g., D - rho2*W)
    double a1,                  // Prior shape for precision (1/rho1)
    double b1,                  // Prior rate for precision (1/rho1)
    int nstar,                  // Effective spatial dimension (p or p-1)
    int tstar                   // Effective temporal dimension (T or T-1)
) {
    // Input checks
    if (invCorrZ.n_rows != invCorrZ.n_cols) {
        Rcpp::stop("posterior_conditional_variance_cpp: invCorrZ must be a square matrix.");
    }
    if (invCorrZ.n_rows != Z.n_rows) {
        Rcpp::stop("posterior_conditional_variance_cpp: Spatial dimensions of Z and invCorrZ must match.");
    }
    if (a1 <= 0 || b1 <= 0) {
         Rcpp::warning("posterior_conditional_variance_cpp: Prior parameters a1, b1 must be positive.");
         // Provide default safe values or stop? Let's use defaults for now.
         if (a1 <= 0) a1 = 0.01;
         if (b1 <= 0) b1 = 0.01;
    }
    if (nstar <= 0 || tstar < 0) { // tstar can be 0 if T=1 and no diff? Check usage. T is columns of Z. If T=1, tstar=1? Original code uses T, not T-1. Let's assume tstar >= 1.
         Rcpp::warning("posterior_conditional_variance_cpp: nstar must be positive, tstar must be non-negative.");
         if (nstar <= 0) nstar = Z.n_rows; // Default to p
         if (tstar < 0) tstar = Z.n_cols; // Default to T
    }


    arma::uword N = Z.n_rows; // p
    arma::uword T = Z.n_cols; // T

    // Posterior shape for precision (1/rho1)
    double g1_pos = static_cast<double>(nstar * tstar) / 2.0 + a1;

    // Calculate quadratic form sum: sum_{t=1}^T Z(:,t)' * invCorrZ * Z(:,t)
    double appo = 0;
    arma::mat invCorrZ_dense = arma::mat(invCorrZ); // Dense for multiplication

    for (arma::uword t = 0; t < T; ++t) {
        appo += arma::as_scalar(Z.col(t).t() * invCorrZ_dense * Z.col(t));
    }

    // Posterior rate for precision (1/rho1)
    // MATLAB: g2_pos = b1 + 0.5 * sum(diag(appo)); -> Since appo is scalar, this is just b1 + 0.5*appo
    double g2_pos = b1 + 0.5 * appo;

     if (g1_pos <= 0 || g2_pos <= 0 || !std::isfinite(g1_pos) || !std::isfinite(g2_pos)) {
        Rcpp::warning("posterior_conditional_variance_cpp: Invalid posterior parameters for Gamma distribution (shape=%.2e, rate=%.2e). Returning prior mean variance.", g1_pos, g2_pos);
        // Return prior mean variance E[rho1] = b1 / (a1 - 1) if a1 > 1, else maybe b1/a1 or stop.
        // Let's return a small fixed variance as fallback.
        return 1e-6;
    }


    // Sample precision h_eta_star ~ Gamma(shape=g1_pos, scale=1/g2_pos)
    // Assuming gamm_rnd(1,1,g1,g2) samples Gamma(shape=g1, scale=1/g2)
    double h_eta_star = R::rgamma(g1_pos, 1.0 / g2_pos);

    if (!std::isfinite(h_eta_star) || h_eta_star <= 1e-12) { // Check for non-positive or non-finite precision
        Rcpp::warning("posterior_conditional_variance_cpp: Sampled precision is non-positive or non-finite (%.2e). Returning small variance.", h_eta_star);
        return 1e-6; // Return a small variance
    }


    // Return variance rho1_star = 1 / h_eta_star
    double rho1_star = 1.0 / h_eta_star;

     if (!std::isfinite(rho1_star) || rho1_star <= 0) {
        Rcpp::warning("posterior_conditional_variance_cpp: Calculated variance is non-positive or non-finite (%.2e). Returning small variance.", rho1_star);
        return 1e-6; // Return a small variance
     }

    return rho1_star;
}

#endif // SPATIAL_CAR_H