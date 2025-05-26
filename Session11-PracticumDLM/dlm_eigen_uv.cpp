#include <RcppEigen.h>
#include <algorithm>
#include <variant>
#include <vector>
#include <simulation_smoother_eigen.h>
#include <cmath> // For std::sqrt, std::log, std::abs
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]



double gamm_rnd(double shape, double scale_inv) {
  // Using R's RNGs via Rcpp
  return R::rgamma(shape, 1.0 / scale_inv);
}


Eigen::VectorXd cpp_prctile(const Eigen::MatrixXd& Y, double P) {
  const int n_rows = Y.rows();
  const int n_cols = Y.cols();
  Eigen::VectorXd quantiles(n_rows);
  
  for (int i = 0; i < n_rows; ++i) {
    std::vector<double> row;
    for (int j = 0; j < n_cols; ++j) {
      double val = Y(i, j);
      if (std::isfinite(val)) {
        row.push_back(val);
      }
    }
    
    if (row.empty()) {
      quantiles(i) = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    
    std::sort(row.begin(), row.end());
    double index = (row.size() - 1) * P;
    size_t lo = static_cast<size_t>(std::floor(index));
    size_t hi = static_cast<size_t>(std::ceil(index));
    double q_lo = row[lo];
    double q_hi = row[hi];
    
    if (index > lo && q_hi != q_lo) {
      double h = index - lo;
      quantiles(i) = (1.0 - h) * q_lo + h * q_hi;
    } else {
      quantiles(i) = q_lo;
    }
  }
  
  return quantiles;
}



Eigen::VectorXd lmfit(const Eigen::MatrixXd& X_data, const Eigen::VectorXd& y_data) {
  return X_data.colPivHouseholderQr().solve(y_data);
}



Eigen::VectorXd extract_last_col(Rcpp::NumericMatrix mat)  {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Eigen::VectorXd out(nrow);
  for (int i = 0; i < nrow; ++i) {
    out(i) = mat(i, ncol - 1);
  }
  return out;
}


Rcpp::NumericVector cube_to_array(const std::vector<Eigen::MatrixXd>& cube_like, const double& n_samples_collected)  {
  int rows = cube_like[0].rows();
  int cols = cube_like[0].cols();
  int slices = static_cast<int>(cube_like.size());
  Rcpp::NumericVector arr(Rcpp::Dimension(rows, cols, slices));
  for (int k = 0; k < slices; ++k) {
    const Eigen::MatrixXd& mat = cube_like[k];
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        int idx = i + j * rows + k * rows * cols;
        arr[idx] = mat(i, j) / n_samples_collected;
      }
    }
  }
  return arr;
}


Rcpp::List minMaxUpperTriangular(const Eigen::MatrixXd& mat) {
  int rows = mat.rows();
  int cols = mat.cols();
  double maxVal = std::numeric_limits<double>::lowest();
  double minNonZero = std::numeric_limits<double>::max();
  bool foundNonZero = false;
  
  for (int i = 0; i < rows; ++i) {
    for (int j = i; j < cols; ++j) {
      double val = mat(i, j);
      if (val > maxVal) {
        maxVal = val;
      }
      if (val != 0.0 && val < minNonZero) {
        minNonZero = val;
        foundNonZero = true;
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("max") = maxVal,
    Rcpp::Named("min_nz") = foundNonZero ? minNonZero : (maxVal / 100.0)
  );
}




// [[Rcpp::export]]
Rcpp::List dlm_cpp(
    const Eigen::MatrixXd& Y,
    const std::vector<Eigen::MatrixXd>& X, // 3D cube: p x t x ncovx
    const Rcpp::Nullable<Rcpp::List> Z_nullable,
    const int& nrep,
    const int& nburn,
    const int& thin,
	const int& print_interval,
    const double& V_beta_0,		  // Prior variance of initial state
    const double& V_gamma,		  // Prior variance of constant coefficients
    const double& a1,             // Prior shape for temporal variance
    const double& b1,             // Prior rate for temporal variance
    const double& s2_a,           // Prior shape for measurement error variance
    const double& s2_b           // Prior rate for measurement error variance
) {
  const int p = X[0].rows();
  const int t = X[0].cols();
  const int ncovx = X.size();
  
  
  Eigen::MatrixXd Z_mat;
  int ncovz = 0;
  
  if (Z_nullable.isNotNull()) {
    Rcpp::List Z_list(Z_nullable);
    ncovz = Z_list.size();
    Z_mat.resize(p * t, ncovz);
    for (int k = 0; k < ncovz; ++k) {
      Eigen::Map<Eigen::MatrixXd> Z_k(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Z_list[k]));
      if (Z_k.rows() != p || Z_k.cols() != t) {
        Rcpp::stop("Dimensions of Z must be p x t x ncovz");
      }
      Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>(Z_k.data(), Z_k.size());
      Z_mat.col(k) = vec;
    }
  }
  
  // Handle NaNs
  Eigen::MatrixXd Y_filled = Y;
  std::vector<std::pair<int, int>> nan_indices;
  std::vector<double> valid_values;
  
  for (int i = 0; i < Y.rows(); ++i) {
    for (int j = 0; j < Y.cols(); ++j) {
      double val = Y(i, j);
      if (std::isfinite(val)) {
        valid_values.push_back(val);
      } else {
        nan_indices.emplace_back(i, j);
      }
    }
  }
  int n_obs = p * t - nan_indices.size();
  
  double mean_Y_valid = 0.0, std_Y_valid = 1.0;
  if (!valid_values.empty()) {
    double sum = std::accumulate(valid_values.begin(), valid_values.end(), 0.0);
    mean_Y_valid = sum / valid_values.size();
    double sq_sum = std::inner_product(valid_values.begin(), valid_values.end(), valid_values.begin(), 0.0);
    std_Y_valid = std::sqrt(sq_sum / valid_values.size() - mean_Y_valid * mean_Y_valid);
  }
  
  if (nan_indices.size() > 0){
	  for (const auto& [i, j] : nan_indices) {
		Y_filled(i, j) = R::rnorm(mean_Y_valid, std_Y_valid);
	  }
  }
  
  
  Eigen::MatrixXd eta_tilde = Y_filled;
  
  
  
  
  // MCMC parameters
  Eigen::VectorXd Q1invdraw_time = Eigen::VectorXd::Constant(ncovx, 10.0); // Precision for T-beta
  double s2_err_mis_draw = 0.01; // Measurement error variance
  
  // Effects: each slice is a matrix
  std::vector<Eigen::MatrixXd> Btimedraw(ncovx, Eigen::MatrixXd::Zero(1, t));   // Time effect T-beta
  
  // Z coefficients and covariate contributions
  Eigen::VectorXd gamma_draw; // Will be resized later if ncovz > 0
  Eigen::MatrixXd meanZ = Eigen::MatrixXd::Zero(p, t); // Contribution of Z covariates
  Eigen::MatrixXd offset = Eigen::MatrixXd::Zero(p, t); // Offset term (unused in Gibbs steps)
  
  Eigen::VectorXd eta_tilde_vec = Eigen::Map<const Eigen::VectorXd>(eta_tilde.data(), eta_tilde.size());
  Eigen::VectorXd meanZ_vec = Eigen::Map<const Eigen::VectorXd>(meanZ.data(), meanZ.size());
  
  
  
  
  
  // Dimensions
  const int q = 1;
  const int Tq = t * q;
  
  
  
  // Construct sparse difference matrix H (Tq x Tq)
  Eigen::SparseMatrix<double> H(Tq, Tq);
  std::vector<Eigen::Triplet<double>> triplets_H;
  for (int i = 0; i < Tq; ++i) {
    triplets_H.emplace_back(i, i, 1.0);
    if (i >= q) {
      triplets_H.emplace_back(i, i - q, -1.0);
    }
  }
  H.setFromTriplets(triplets_H.begin(), triplets_H.end());
  
  
  
  std::vector<Eigen::SparseMatrix<double>> bigG(ncovx);   // pt x t
  
  for (int k = 0; k < ncovx; ++k) {
    std::vector<Eigen::Triplet<double>> triplets_G;
    
    for (int i = 0; i < t; ++i) {
      const Eigen::VectorXd& x_col = X[k].col(i);  // p x 1
      
      // bigG: pt x t (each column is stacked x_col)
      for (int j = 0; j < p; ++j) {
        triplets_G.emplace_back(i * p + j, i, x_col(j));
      }
      
    }
    
    Eigen::SparseMatrix<double> G(p * t, t);
    G.setFromTriplets(triplets_G.begin(), triplets_G.end());
    bigG[k] = G;
    
  }
  
  
  // --- Storage Initialization ---
  if (nrep % thin != 0) {
    Rcpp::warning("nrep is not a multiple of thin, some iterations will be discarded.");
  }
  int collections = nrep / thin;
  int MCMC_samples = nrep + nburn;
  int collect_count = 0;
  
  // Storage matrices
  Rcpp::List out_results;
  Eigen::MatrixXd S2_err_mis_(1, collections);
  Eigen::MatrixXd Q1inv_time_out(ncovx, collections);
  Eigen::MatrixXd G_out;  // Will be resized only if ncovz > 0
  // Eigen::VectorXd store_llike = Eigen::VectorXd::Zero(collections);
	std::vector<Eigen::MatrixXd> store_llike(collections, Eigen::MatrixXd::Zero(p, t));
  
  if (ncovz > 0) {
    G_out.resize(ncovz, collections);
  }
  
  std::vector<Eigen::MatrixXd> YPRED_out(collections, Eigen::MatrixXd::Zero(p, t));
  Eigen::MatrixXd RMSE_(1, collections);
  Eigen::MatrixXd MAE_(1, collections);
  Eigen::VectorXd chi_sq_pred_ = Eigen::VectorXd::Zero(collections);
  Eigen::VectorXd chi_sq_obs_ = Eigen::VectorXd::Zero(collections);
  
  
  // Averaging structures
  Eigen::MatrixXd Ypred_mean = Eigen::MatrixXd::Zero(p, t);
  Eigen::MatrixXd Ypred2_mean = Eigen::MatrixXd::Zero(p, t); // Sum of squares for variance
  Eigen::MatrixXd Eta_tilde_mean = Eigen::MatrixXd::Zero(p, t);
  Eigen::MatrixXd thetay_mean = Eigen::MatrixXd::Zero(p, t);
  Eigen::MatrixXd meanZmean = Eigen::MatrixXd::Zero(p, t);
  Eigen::MatrixXd meanY1mean = Eigen::MatrixXd::Zero(p, t);
  
  std::vector<Eigen::MatrixXd> Btime_postmean(ncovx, Eigen::MatrixXd::Zero(1, t));
  std::vector<Eigen::MatrixXd> Btime2_postmean(ncovx, Eigen::MatrixXd::Zero(1, t));
  
	
  // Diagnostics
  Eigen::VectorXd pvalue_ResgrRespred_sum = Eigen::VectorXd::Zero(valid_values.size());
  Eigen::VectorXd pvalue_YgrYhat_sum = Eigen::VectorXd::Zero(valid_values.size());
  double chisq_count = 0.0;
  double quant = 0.95;
  Eigen::VectorXd p95_obs_ = cpp_prctile(Y, quant);
  Eigen::VectorXd percentile95_sum = Eigen::VectorXd::Zero(p);
  double store_CRPS_1_sum = 0.0;
  double store_CRPS_2_sum = 0.0;
  
  
  
  // --- MCMC Loop ---
  Rcpp::Rcout << "Starting MCMC (" << MCMC_samples << " iterations)..." << std::endl;
  
  for (int irep = 0; irep < MCMC_samples; ++irep) {
    // Print progress
    if ((irep + 1) % print_interval == 0) {
      Rcpp::Rcout << "Iteration: " << irep + 1 << " / " << MCMC_samples << std::endl;
      Rcpp::checkUserInterrupt(); // Allow user to interrupt
    }
    
    
    // Step I.a: Sampling T-beta (time effects)
    // Compute residual y_k
    Eigen::VectorXd y2_vec = eta_tilde_vec - meanZ_vec;
    for (int k = 0; k < ncovx; ++k) {
      Eigen::VectorXd y_k = y2_vec;
      for (int k2 = 0; k2 < ncovx; ++k2) {
        if (k == k2) continue;
        y_k -= bigG[k2] * Eigen::Map<const Eigen::VectorXd>(Btimedraw[k2].data(), t);
      }
      
      // Construct prior precision matrix K = H^T * invS * H
      Eigen::SparseMatrix<double> invS_diag(Tq, Tq);
      std::vector<Eigen::Triplet<double>> invS_triplets;
      for (int i = 0; i < Tq; ++i) {
        invS_triplets.emplace_back(i, i, Q1invdraw_time(k));
      }
      invS_triplets[0] = Eigen::Triplet<double>(0, 0, 1.0 / V_beta_0);  // diffuse prior
      invS_diag.setFromTriplets(invS_triplets.begin(), invS_triplets.end());
      
      Eigen::SparseMatrix<double> K = H.transpose() * invS_diag * H;
      
      // Observation precision
      Eigen::SparseMatrix<double> GinvOmega11 = bigG[k].transpose() * (1.0 / s2_err_mis_draw);
      Eigen::SparseMatrix<double> GinvOmega11G = GinvOmega11 * bigG[k];
      
      // Posterior precision invP = K + G'SigmaInvG
      Eigen::SparseMatrix<double> invP = K + GinvOmega11G;
      invP = 0.5 * (invP + Eigen::SparseMatrix<double>(invP.transpose()));  // Ensure symmetry
      
      // Posterior mean
      Eigen::VectorXd tmp = GinvOmega11 * y_k;
      
      // Sample from posterior
      Eigen::MatrixXd bb = rmvnorm_eigen_sparse(1, invP, tmp);  // T x 1
      
      // Reshape
      Eigen::MatrixXd Bdrawc = Eigen::Map<Eigen::MatrixXd>(bb.data(), 1, t);
      Btimedraw[k] = Bdrawc;
      
    }
    
    
    
    
    // Step I.d: Sampling constant coefficients (C-beta and Z-gamma)
    // Compute residual y2 = eta_tilde - T*beta_t - S*beta_s - ST*beta_st
    Eigen::VectorXd y2_reg = eta_tilde_vec;
    for (int k = 0; k < ncovx; ++k) {
      y2_reg -= bigG[k] * Btimedraw[k].transpose();
    }
		
    // Construct combined regressor matrix [Z, X]
    Eigen::MatrixXd X_reshaped(p * t, ncovx);
    for (int k = 0; k < ncovx; ++k) {
      X_reshaped.col(k) = Eigen::Map<const Eigen::VectorXd>(X[k].data(), p * t);
    }
    
    if (ncovz > 0) {    
			// Prior precision (diagonal, small values for diffuse prior)
			Eigen::MatrixXd K_reg = Eigen::MatrixXd::Identity(ncovz, ncovz) * (1.0 / V_gamma);
			
			// Posterior precision and mean
			Eigen::MatrixXd GinvOmega11_reg = Z_mat.transpose() * (1.0 / s2_err_mis_draw);
			Eigen::MatrixXd GinvOmega11G_reg = GinvOmega11_reg * Z_mat;
			Eigen::MatrixXd invP_reg = K_reg + GinvOmega11G_reg;
			invP_reg = 0.5 * (invP_reg + invP_reg.transpose());  // Ensure symmetry
			
			Eigen::VectorXd tmp_reg = GinvOmega11_reg * y2_reg;
    
			// Sample from posterior
			Eigen::MatrixXd P_reg = invP_reg.inverse();
			Eigen::LLT<Eigen::MatrixXd> llt(P_reg);
			Eigen::MatrixXd L_reg = llt.matrixL();
			gamma_draw = P_reg * tmp_reg;
			Eigen::VectorXd normal_std_reg = Eigen::VectorXd::Zero(ncovz);
			for (int i = 0; i < ncovz; ++i) {
					normal_std_reg(i) += R::rnorm(0.0, 1.0);
				}
			gamma_draw += L_reg * normal_std_reg;
      Eigen::VectorXd Zgamma = Z_mat * gamma_draw;
      meanZ = Eigen::Map<const Eigen::MatrixXd>(Zgamma.data(), p, t);
      meanZ_vec = Zgamma;
    }
    
    
    // Step I.e: Computing the current mean component from X
    Eigen::MatrixXd current_meanY1 = Eigen::MatrixXd::Zero(p, t);
    
    for (int k = 0; k < ncovx; ++k) {
      Eigen::MatrixXd B_t_full = Btimedraw[k].colwise().replicate(p);  // 1 x t → p x t
      Eigen::ArrayXXd appoX_B = X[k].array() * B_t_full.array();
      current_meanY1 += appoX_B.matrix();
    }
    
    
    // Step II: Sampling variances of state vectors
    for (int k = 0; k < ncovx; ++k) {
      // --- TIME (T-beta) ---
      Eigen::MatrixXd Btime_k = Btimedraw[k];  // 1 x t
      Eigen::VectorXd e2 = Btime_k.row(0).segment(1, t - 1) - Btime_k.row(0).segment(0, t - 1);
      double newnu2 = a1 + static_cast<double>(t - 1) / 2.0;
      double newS2 = b1 + 0.5 * e2.squaredNorm();
      Q1invdraw_time(k) = R::rgamma(newnu2, 1.0 / newS2);
      
      
    }
    
    
    // Step III: Sampling the measurement error variance
    Eigen::MatrixXd thetay_draw = current_meanY1 + meanZ + offset;
    Eigen::MatrixXd yhat = eta_tilde - thetay_draw;
    
    double sse_2 = yhat.array().square().sum();
    double g1_pos = s2_a + static_cast<double>(p * t) / 2.0;
    double g2_pos = s2_b + 0.5 * sse_2;
    
    double precision_draw = R::rgamma(g1_pos, 1.0 / g2_pos);
    s2_err_mis_draw = 1.0 / precision_draw;
    
    
    // Step IV: Handling missing data (imputation)
    // Generate predicted Y ~ N(thetay_draw, s2_err_mis_draw * I)
    double current_sd = std::sqrt(s2_err_mis_draw);
    auto norm_mcmc = [&] (double) {return R::rnorm(0.0, current_sd);};
    Eigen::MatrixXd Ypred = Eigen::MatrixXd::NullaryExpr(p, t, norm_mcmc);
    Ypred += thetay_draw;
    
    // Impute missing values in eta_tilde
		if (nan_indices.size() > 0) {
			for (const auto& [i, j] : nan_indices) {
				eta_tilde(i, j) = Ypred(i, j);
			}
		}
    eta_tilde_vec = Eigen::Map<const Eigen::VectorXd>(eta_tilde.data(), eta_tilde.size());
    
    
    // --- Store results post-burn-in and thinning ---
    if (irep >= nburn && (irep - nburn + 1) % thin == 0) {
      // Store posterior samples
      S2_err_mis_(0, collect_count) = s2_err_mis_draw;
      
      for (int k = 0; k < ncovx; ++k) {
        Q1inv_time_out(k, collect_count) = Q1invdraw_time(k);
      }
      
      if (ncovz > 0) {
        for (int j = 0; j < ncovz; ++j) {
          G_out(j, collect_count) = gamma_draw(j);
        }
      }
      
      // Log-likelihood
      /*double log_det_Sigma = static_cast<double>(p) * std::log(s2_err_mis_draw);
      double term1 = -0.5 * sse_2 / s2_err_mis_draw;
      double term2 = -0.5 * static_cast<double>(p * t) * std::log(2.0 * M_PI);
      double term3 = -0.5 * log_det_Sigma;
      store_llike(collect_count) = term1 + term2 + term3;*/
			// Calculate point-wise log-likelihood
			Eigen::MatrixXd term1 = -0.5 * yhat.array().square().matrix() / s2_err_mis_draw;
			Eigen::MatrixXd term2 = Eigen::MatrixXd::Constant(p, t, -0.5 * std::log(2.0 * M_PI));
			Eigen::MatrixXd term3 = Eigen::MatrixXd::Constant(p, t, -0.5 * std::log(s2_err_mis_draw));
			store_llike[collect_count] = term1 + term2 + term3;

      Eigen::MatrixXd Ypred2 = Eigen::MatrixXd::NullaryExpr(p, t, norm_mcmc);
			Ypred2 += thetay_draw;
      YPRED_out[collect_count] = Ypred;
      
      // Diagnostics
      if (n_obs>0) {
        Eigen::VectorXd Y_obs(n_obs), Ypred_obs(n_obs), Ypred_2_obs(n_obs), theta_obs(n_obs);
        int idx = 0;
        
        for (int i = 0; i < Y.rows(); ++i) {
          for (int j = 0; j < Y.cols(); ++j) {
            if (std::isfinite(Y(i, j))) {
              Y_obs(idx) = Y(i, j);
              Ypred_obs(idx) = Ypred(i, j);
              Ypred_2_obs(idx) = Ypred2(i, j);
              theta_obs(idx) = thetay_draw(i, j);
              ++idx;
            }
          }
        }
        
        Eigen::VectorXd V_hat_vec = Eigen::VectorXd::Constant(n_obs, s2_err_mis_draw);
        Eigen::VectorXd pearson_res = (Y_obs - theta_obs).array() / V_hat_vec.array().sqrt();
        Eigen::VectorXd pearson_res_pred = (Ypred_obs - theta_obs).array() / V_hat_vec.array().sqrt();
        
        chi_sq_obs_(collect_count) = pearson_res.squaredNorm();
        chi_sq_pred_(collect_count) = pearson_res_pred.squaredNorm();
				
				if (std::isfinite(chi_sq_obs_(collect_count)) && std::isfinite(chi_sq_pred_(collect_count))) {
					if(chi_sq_obs_(collect_count) >= chi_sq_pred_(collect_count)) {
						chisq_count += 1.0;
					}
				}
				
				for (int i = 0; i < n_obs; ++i) {
					if (std::pow(pearson_res(i), 2.0) >= std::pow(pearson_res_pred(i), 2.0)) {
						pvalue_ResgrRespred_sum(i) += 1.0;
					}
					if (Y_obs(i) >= Ypred_obs(i)) {
						pvalue_YgrYhat_sum(i) += 1.0;
					}
				}

				// CRPS components
				store_CRPS_1_sum += (Ypred_obs - Ypred_2_obs).array().abs().sum();
				store_CRPS_2_sum += (Ypred_obs - Y_obs).array().abs().sum();

        
        RMSE_(0, collect_count) = std::sqrt((Y_obs - Ypred_obs).array().square().mean());
        MAE_(0, collect_count) = (Y_obs - Ypred_obs).array().abs().mean();
      } else {
        RMSE_(0, collect_count) = std::numeric_limits<double>::quiet_NaN();
        MAE_(0, collect_count) = std::numeric_limits<double>::quiet_NaN();
        chi_sq_obs_(collect_count) = std::numeric_limits<double>::quiet_NaN();
        chi_sq_pred_(collect_count) = std::numeric_limits<double>::quiet_NaN();
      }
      
      Eigen::VectorXd p95_pred = cpp_prctile(Ypred, quant);
      
      // Element-wise comparison and accumulation
      for (int i = 0; i < p; ++i) {
        if (p95_obs_(i) >= p95_pred(i)) {
          percentile95_sum(i) += 1.0;
        }
      }
      
      
      // Update averages ('ave' structure)
      Ypred_mean += Ypred;
      Ypred2_mean += Ypred.array().square().matrix();
      Eta_tilde_mean += eta_tilde;
      thetay_mean += thetay_draw;
      meanZmean += meanZ;
      meanY1mean += current_meanY1;
      
      for (int k = 0; k < ncovx; ++k) {        
        Btime_postmean[k] += Btimedraw[k];
        Btime2_postmean[k] += Btimedraw[k].array().square().matrix();
        
      }
      
      
      
      collect_count++;
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
  ave_results["meanZmean"] = meanZmean / n_samples_collected;
  ave_results["meanY1mean"] = meanY1mean / n_samples_collected;
  
  // Convert vector of matrices to 3D arrays for R
  ave_results["Btime_postmean"] = cube_to_array(Btime_postmean, n_samples_collected);
  ave_results["Btime2_postmean"] = cube_to_array(Btime2_postmean, n_samples_collected);
  
  // DIC calculation
	Eigen::VectorXd store_llike_vec(collections);
	for (int i = 0; i < collect_count; ++i) {
		store_llike_vec(i) = store_llike[i].sum();
	}
  double mean_llike = store_llike_vec.mean();
  double D_bar = -2.0 * mean_llike;
  
  // Extract posterior means
  Eigen::MatrixXd theta_hat = thetay_mean / n_samples_collected;;
  Eigen::MatrixXd eta_tilde_hat = Eta_tilde_mean / n_samples_collected;
  double s2_err_hat = S2_err_mis_.row(0).mean();
  
  Eigen::MatrixXd resid_hat = eta_tilde_hat - theta_hat;
  double sse_hat = resid_hat.array().square().sum();
  double log_det_Sigma_hat = static_cast<double>(p) * std::log(s2_err_hat);
  
  double term1_hat = -0.5 * sse_hat / s2_err_hat;
  double term2_hat = -0.5 * static_cast<double>(p * t) * std::log(2.0 * M_PI);
  double term3_hat = -0.5 * log_det_Sigma_hat;
  double log_lik_hat = term1_hat + term2_hat + term3_hat;
  double D_hat = -2.0 * log_lik_hat;
  
	ave_results["Dbar"] = D_bar;
	ave_results["pD"] = D_bar - D_hat;
  ave_results["DIC"] = D_bar + (D_bar - D_hat);  // DIC = D_bar + pD = 2*D_bar - D_hat
  
  
  // Calculate other summary stats for 'ave'
  if (n_obs>0) {
    ave_results["pvalue_ResgrRespred"] = pvalue_ResgrRespred_sum.mean() / n_samples_collected;
    ave_results["pvalue_YgrYhat"] = pvalue_YgrYhat_sum.mean() / n_samples_collected;
    ave_results["chisq_pvalue"] = chisq_count / n_samples_collected;
    ave_results["CRPS"] = (0.5 * store_CRPS_1_sum - store_CRPS_2_sum) / (n_samples_collected * n_obs);
  } else {
    ave_results["pvalue_ResgrRespred"] = NA_REAL;
    ave_results["pvalue_YgrYhat"] = NA_REAL;
    ave_results["chisq_pvalue"] = NA_REAL;
    ave_results["CRPS"] = NA_REAL;
  }
  
  // Percentile p-value vector
  ave_results["percentile95_pvalue"] = Rcpp::wrap((percentile95_sum / n_samples_collected).eval());
  
  // PMCC: Posterior Mean + Variance (using only non-missing Y)
  Eigen::MatrixXd Ypred_mean_ave = Ypred_mean / n_samples_collected;
  Eigen::MatrixXd Ypred2_mean_ave = Ypred2_mean / n_samples_collected;
  Eigen::MatrixXd varYPRED = Ypred2_mean_ave - Ypred_mean_ave.array().square().matrix();
  
  double pmcc_sum = 0.0;
  if (n_obs > 0) {
    // for (const auto& [i, j] : nan_indices) continue;  // skip missing
    for (int i = 0; i < Y.rows(); ++i) {
      for (int j = 0; j < Y.cols(); ++j) {
        if (std::isfinite(Y(i, j))) {
          double diff = Y(i, j) - Ypred_mean_ave(i, j);
          pmcc_sum += diff * diff + varYPRED(i, j);
        }
      }
    }
  }
  ave_results["PMCC"] = pmcc_sum;
  
  // Finalize 'out' list
  out_results["S2_err_mis_"] = S2_err_mis_;
  out_results["Q1inv_time"] = Q1inv_time_out;
  if (ncovz > 0) {
    out_results["G_"] = G_out;
  }
  out_results["store_llike"] = cube_to_array(store_llike, 1.0);
  out_results["YPRED_"] = cube_to_array(YPRED_out, 1.0);
	
  out_results["RMSE_"] = RMSE_;
  out_results["MAE_"] = MAE_;
  out_results["chi_sq_pred_"] = chi_sq_pred_;
  out_results["chi_sq_obs_"] = chi_sq_obs_;
  
  // Store final state for potential restart
  out_results["eta_tilde"] = eta_tilde;
  out_results["Btimedraw"] = cube_to_array(Btimedraw, 1.0);
  
  
  
  // Return both lists
  return Rcpp::List::create(Rcpp::Named("ave") = ave_results,
                            Rcpp::Named("out") = out_results);
  
}

