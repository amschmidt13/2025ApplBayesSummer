// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <simulation_smoother.h>
#include <cmath> // For std::sqrt, std::log, std::abs


double gamm_rnd(double shape, double scale_inv) {
    // Using R's RNGs via Rcpp
    return R::rgamma(shape, 1.0 / scale_inv);
}



// Placeholder for prctile(Y, 95, 2) -> calculates 95th percentile row-wise
arma::vec cpp_prctile(const arma::mat& Y, const arma::vec& P) {
	arma::vec quantiles_by_row(Y.n_rows);
	for(arma::uword j = 0; j < Y.n_rows; ++j) {
        arma::vec y = Y.row(j).t();
        quantiles_by_row(j) = arma::as_scalar(arma::quantile(y.elem(arma::find_finite(y)), P));
    }
    
    return quantiles_by_row;
}


// Returns OLS coefficient vector
arma::vec lmfit(const arma::mat& X_data, const arma::vec& y_data) {
    // Simple OLS solution: (X'X)^(-1) * X'y
    return arma::solve(X_data, y_data);
}


arma::cube repcube(const arma::cube& A, int s1, int s2 ) {
  int Nr = A.n_rows;
  int Nc = A.n_cols;
  int Ns = A.n_slices;
  arma::cube out(Nr*s1, Nc*s2, Ns);
  for (int i = 0; i < Ns; ++i) {
    out.slice(i) = arma::repmat(A.slice(i), s1, s2);
  }
  return out;
}



// [[Rcpp::export]]
Rcpp::List dlm_cpp(
    const arma::mat& Y,           // Input matrix p x t
    const arma::cube& X,          // Input cube p x t x ncovx
    const Rcpp::Nullable<arma::cube> Z_nullable, // Input cube p x t x ncovz (optional)
    const int& nrep,
    const int& nburn,
    const int& thin,
	const double& V_beta_0,		  // Prior variance of initial state
	const double& V_gamma,		  // Prior variance of constant coefficients
    const double& a1,             // Prior shape for temporal variance
    const double& b1,             // Prior rate for temporal variance
    const double& s2_a,           // Prior shape for measurement error variance
    const double& s2_b            // Prior rate for measurement error variance
) {
    // --- Input Processing ---
    arma::uword p = X.n_rows;
    arma::uword t = X.n_cols;
    arma::uword ncovx = X.n_slices;
	
    arma::mat Z_mat; // Will be p*t x ncovz
    arma::uword ncovz = 0;
    bool Z_is_null = Z_nullable.isNull();

    if (!Z_is_null) {
        arma::cube Z_cube(Rcpp::as<arma::cube>(Z_nullable));
        ncovz = Z_cube.n_slices;
        if (Z_cube.n_rows != p || Z_cube.n_cols != t) {
            Rcpp::stop("Dimensions of Z must be p x t x ncovz");
        }
        // Reshape Z: p x t x ncovz -> (pt) x ncovz
        Z_mat.set_size(p * t, ncovz);
        for(arma::uword k = 0; k < ncovz; ++k) {
            arma::mat Z_slice = Z_cube.slice(k);
            Z_mat.col(k) = arma::vectorise(Z_slice);
        }
    }

    
    // Handle NaNs
    arma::umat nan_mask = arma::find_nonfinite(Y); // Indices of NaNs
    arma::mat Y_filled = Y; // Create a copy to fill NaNs
    double mean_Y_valid = 0;
    double std_Y_valid = 1.0;
    arma::vec Y_valid_vec = Y.elem(arma::find_finite(Y));
    if (Y_valid_vec.n_elem > 0) {
         mean_Y_valid = arma::mean(Y_valid_vec);
         std_Y_valid = arma::stddev(Y_valid_vec);
    }
    Y_filled.elem(nan_mask) = arma::randn(nan_mask.n_elem) * std_Y_valid + mean_Y_valid;
    arma::mat eta_tilde = Y_filled; // Initialize eta_tilde

    // MCMC parameters
    arma::vec Q1invdraw_time = 100000.0 * arma::ones<arma::vec>(ncovx); // Precision for T-beta
    double s2_err_mis_draw = 0.000001; // Measurement error variance

    arma::cube Btimedraw = arma::zeros<arma::cube>(1, t, ncovx); // Time effect T-beta
    arma::vec gamma_draw; // Coefficients for Z
    arma::mat meanZ = arma::zeros<arma::mat>(p, t); // Contribution of Z covariates

    


    // --- Precompute Matrices ---
    arma::uword q = 1; // Dimension of time states
    arma::uword Tq = t * q;

    // Difference operator H for temporal effects
    arma::sp_mat H(Tq, Tq);
    H.eye();
    for (arma::uword i = q; i < Tq; ++i) {
        H(i, i - q) = -1.0;
    }


    // Constructing big observation matrices (potentially large and sparse)
    // These map state vectors to observations
    std::vector<arma::sp_mat> bigG(ncovx);  // For T-beta (time)

    for (arma::uword k = 0; k < ncovx; ++k) {
        //arma::field<arma::sp_mat> temp_list_G(t);
        arma::field<arma::vec> temp_list_G(t);

        for (arma::uword i = 0; i < t; ++i) {
            arma::vec x_col = X.slice(k).col(i); // p x 1
            // For G: X maps a single time state (beta_t) to p observations
            //temp_list_G(i) = arma::sp_mat(x_col); // p x 1 sparse matrix
            temp_list_G(i) = x_col; // p x 1 sparse matrix
        }
		bigG[k] = arma::sp_mat(p * t, t); // Initialize sparse
		for(arma::uword i=0; i<t; ++i) {
			bigG[k].submat(i*p, i, (i+1)*p - 1, i) = temp_list_G(i);
		}
    }


    // --- Storage Initialization ---
	if (nrep % thin != 0) {
      Rcpp::warning("nrep is not a multiple of thin, some iterations will be discarded.");
    }
    int collections = floor(nrep / thin); // Assicura intero
    int MCMC_samples = nrep + nburn;
    int collect_count = 0;
	

    // Output list structure
    Rcpp::List out_results;
    arma::mat S2_err_mis_(1, collections);
    arma::mat Q1inv_time_out(ncovx, collections); // Store precision for T-beta
    arma::mat G_out; // For Z coefficients
    if (ncovz > 0) {
        G_out.set_size(ncovz, collections);
    }
    arma::cube store_llike(p, t, collections);
    arma::cube YPRED_out(p, t, collections);


    // Averaging structure
    arma::mat Ypred_mean(p, t, arma::fill::zeros);
    arma::mat Ypred2_mean(p, t, arma::fill::zeros); // Sum of squares for variance
    arma::mat Eta_tilde_mean(p, t, arma::fill::zeros);
    arma::mat thetay_mean(p, t, arma::fill::zeros);
    arma::cube Btime_postmean(1, t, ncovx, arma::fill::zeros);
    arma::cube Btime2_postmean(1, t, ncovx, arma::fill::zeros);



    // --- MCMC Loop ---
    Rcpp::Rcout << "Starting MCMC (" << MCMC_samples << " iterations)..." << std::endl;
    int print_interval = 100;

    for (int irep = 0; irep < MCMC_samples; ++irep) {

        // Print progress
		if ((irep + 1) % print_interval == 0) {
			Rcpp::Rcout << "Iteration: " << irep + 1 << " / " << MCMC_samples << std::endl;
			Rcpp::checkUserInterrupt(); // Allow user to interrupt
		}

        arma::mat current_meanY1 = arma::zeros<arma::mat>(p,t);
        // --- STEP I.a: Sample T-beta (Time effects) ---
        arma::vec y2_vec = arma::vectorise(eta_tilde - meanZ); // pt x 1

        for (arma::uword k = 0; k < ncovx; ++k) {
            // Adjust residual: remove effects of other time components
            arma::vec y_k = y2_vec;
            for (arma::uword k2 = 0; k2 < ncovx; ++k2) {
                if (k == k2) continue;
                y_k -= bigG[k2] * arma::vectorise(Btimedraw.slice(k2).row(0)); // Btimedraw is 1 x t x ncovx
            }

            // Precision components
            double invOmega22_val = Q1invdraw_time(k); // Precision of innovations
            // arma::sp_mat invOmega22 = arma::sp_mat(q, q); // q=1
            // invOmega22(0,0) = invOmega22_val;
            arma::sp_mat invS_diag(Tq, Tq); // Precision of state vector [beta_0, ..., beta_{t-1}]
            invS_diag.diag() += invOmega22_val;
            invS_diag(0,0) = 1.0 / V_beta_0; // Prior for beta_0 (diffuse N(0, 10)) -> precision 1/10

            // Total precision K = H' * invS * H
            arma::sp_mat K = H.t() * invS_diag * H;
            arma::sp_mat GinvOmega11 = bigG[k].t() / s2_err_mis_draw;
            arma::sp_mat GinvOmega11G = GinvOmega11 * bigG[k];

            // Posterior precision invP = K + G'SigmaInvG
            arma::sp_mat invP = K + GinvOmega11G;
            invP = 0.5 * (invP + invP.t()); // Ensure symmetry

			arma::vec tmp = GinvOmega11 * y_k;
			arma::vec bb = arma::vectorise(rmvnorm_arma_solve(1, arma::mat(invP), tmp));

            // Reshape
            arma::mat Bdrawc = arma::reshape(bb, q, t); // Should be 1 x t
            Btimedraw.slice(k) = Bdrawc;
        }


        // --- STEP I.d: Sample constant coefficients (Z-gamma) ---


        // Construct combined regressor matrix [Z, X_1, X_2, ...]
        if (ncovz > 0) {
			arma::vec y2_reg = arma::vectorise(eta_tilde); // pt x 1
			for (arma::uword k = 0; k < ncovx; ++k) {
				y2_reg -= bigG[k] * Btimedraw.slice(k).row(0).t();
			}
			
			// Prior precision (diagonal, small values for diffuse prior)
			arma::mat K_reg(ncovz, ncovz);
			K_reg.eye();
			K_reg *= (1.0/V_gamma);
			arma::mat GinvOmega11_reg = Z_mat.t() / s2_err_mis_draw;
			arma::mat GinvOmega11G_reg = GinvOmega11_reg * Z_mat;

			// Posterior precision invP = K + G'SigmaInvG
			arma::mat invP_reg = K_reg + GinvOmega11G_reg;
			invP_reg = 0.5 * (invP_reg + invP_reg.t()); // Ensure symmetry
			
			arma::vec tmp_reg = GinvOmega11_reg * y2_reg;
			arma::mat P_reg = arma::inv(invP_reg);
			gamma_draw = arma::vectorise(arma::mvnrnd(P_reg*tmp_reg, P_reg));
			arma::vec Zgamma = Z_mat * gamma_draw;
            meanZ = arma::reshape(Zgamma, p, t); // Update mean contribution from Z
        }

        

        // --- Calculate current mean component from X ---
        current_meanY1.zeros();
		for (arma::uword k = 0; k < ncovx; ++k) {
			// Create full p x t matrices for each component for covariate k
			current_meanY1 += X.slice(k) % (Btimedraw.slice(k).row(0)); // Element-wise product
		}


        // --- STEP II: Sample variances of state vectors ---
        for (arma::uword k = 0; k < ncovx; ++k) {
            // TIME (T-beta): Sample precision Q1invdraw_time(k) from Gamma
            // Prior: IG(a1, b1) for variance -> Gamma(a1, b1) for precision
            double nu02 = a1;
            double S02 = b1;
            arma::mat Btime_k = Btimedraw.slice(k); // 1 x t
            arma::mat e2 = Btime_k.cols(1, t - 1) - Btime_k.cols(0, t - 2); // Differences (1 x t-1)
            double newnu2 = nu02 + static_cast<double>(t - 1) / 2.0; // Shape (assuming q=1)
            double newS2 = S02 + 0.5 * arma::accu(arma::pow(e2, 2));   // Rate parameter
			Q1invdraw_time(k) = R::rgamma(newnu2, 1.0 / newS2); // Sample precision
        }

        // --- Step III: Sample variance of the measurement error (s2_err_mis_draw) ---
        // Prior: IG(s2_a, s2_b) -> Precision ~ Gamma(s2_a, s2_b)
        arma::mat thetay_draw = current_meanY1 + meanZ; // Current estimate of the mean E[Y]
        arma::mat yhat = eta_tilde - thetay_draw; // Residuals
        double sse_2 = arma::accu(arma::pow(yhat, 2)); // Sum of squared errors

        double g1_pos = s2_a + static_cast<double>(p * t) / 2.0; // Posterior shape for precision
        double g2_pos = s2_b + 0.5 * sse_2;                   // Posterior rate for precision

        double precision_draw = R::rgamma(g1_pos, 1.0 / g2_pos);
        s2_err_mis_draw = 1.0 / precision_draw;


        // --- Step IV: Handle missing data (impute Y based on current model) ---
        // Ypred ~ N(thetay_draw, s2_err_mis_draw * I)
        double current_sd = std::sqrt(s2_err_mis_draw);
        arma::mat Ypred = arma::randn<arma::mat>(p, t) * current_sd + thetay_draw;

        // Fill only the originally missing values in eta_tilde
        if (nan_mask.n_elem > 0) {
			eta_tilde.elem(nan_mask) = Ypred.elem(nan_mask);
        }
		
		// Calculate point-wise log-likelihood
		arma::mat term1 = -0.5 * arma::pow(yhat, 2) / s2_err_mis_draw; // -0.5 * (eta-theta)' * inv(s2*I) * (eta-theta)
        arma::mat term2 = arma::ones<arma::mat>(p, t) * (-0.5 * std::log(2.0 * arma::datum::pi)); // Constant part
        arma::mat term3 = arma::ones<arma::mat>(p, t) * (-0.5 * std::log(s2_err_mis_draw));
		store_llike.slice(collect_count) = term1 + term2 + term3;

        // --- Store results post-burn-in and thinning ---
        if (irep >= nburn && (irep - nburn + 1) % thin == 0) {
            // Store posterior samples in 'out' structure
            S2_err_mis_(0, collect_count) = s2_err_mis_draw;
            Q1inv_time_out.col(collect_count) = Q1invdraw_time; // Store precision

            if (ncovz > 0) {
                G_out.col(collect_count) = gamma_draw;
            }

            YPRED_out.slice(collect_count) = Ypred;

            // Update averages ('ave' structure)
            Ypred_mean += Ypred;
            Ypred2_mean += arma::pow(Ypred, 2);
            Eta_tilde_mean += eta_tilde;
            thetay_mean += thetay_draw;
            Btime_postmean += Btimedraw;
            Btime2_postmean += arma::pow(Btimedraw, 2);

            collect_count++; // Increment collection counter
        }

    } // End MCMC loop

    Rcpp::Rcout << "MCMC finished." << std::endl;

    // --- Post-processing and calculating averages ---
    double n_samples_collected = static_cast<double>(collect_count);

    // Finalize 'ave' list
    Rcpp::List ave_results;
    ave_results["Ypred_mean"] = Ypred_mean / n_samples_collected;
    ave_results["Ypred2_mean"] = Ypred2_mean / n_samples_collected;
    ave_results["Eta_tilde_mean"] = Eta_tilde_mean / n_samples_collected;
    ave_results["thetay_mean"] = thetay_mean / n_samples_collected;
    ave_results["Btime_postmean"] = Btime_postmean / n_samples_collected;
    ave_results["Btime2_postmean"] = Btime2_postmean / n_samples_collected;


    // Finalize 'out' list
    out_results["S2_err_mis_"] = S2_err_mis_;
    out_results["Q1inv_time"] = Q1inv_time_out; // Store precision samples
    if (ncovz > 0) {
        out_results["G_"] = G_out;
    }
	out_results["YPRED_out"] = YPRED_out;
	out_results["store_llike"] = store_llike;

    // Store final state for potential restart
    out_results["eta_tilde"] = eta_tilde; // Last state of imputed data
    // Store last sampled values of parameters
    out_results["Btimedraw"] = Btimedraw;


    // Return both lists
    return Rcpp::List::create(Rcpp::Named("ave") = ave_results,
                              Rcpp::Named("out") = out_results);
}

