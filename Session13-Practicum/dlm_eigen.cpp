// 2025
// Code prepared by Carlo Zaccardi

#include <RcppEigen.h>
#include <algorithm>
#include <variant>
#include <vector>
#include <simulation_smoother_eigen.h>
#include <spatial_car_eigen.h>
#include <spatial_exp_eigen.h>
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
    const Eigen::MatrixXd& W_dense,
    const bool& point_referenced,
    const int& nrep,
    const int& nburn,
    const int& thin,
	const int& print_interval,
    const double& V_beta_0,		  // Prior variance of initial state
    const double& V_gamma,		  // Prior variance of constant coefficients
    const double& a1,
    const double& b1,
    const double& s2_a,
    const double& s2_b,
    const bool keepYpred,
    Rcpp::Nullable<Rcpp::List> out_prev_nullable
) {
  const int p = X[0].rows();
  const int t = X[0].cols();
  const int ncovx = X.size();
  
  Eigen::SparseMatrix<double> W;
  double min_rho2, max_rho2;
  if (!point_referenced) {
    W = W_dense.sparseView();
    min_rho2 = 0.1;
    max_rho2 = 1.0;
  } else {
    Rcpp::List min_and_max = minMaxUpperTriangular(W_dense);
    min_rho2 = Rcpp::as<double>(min_and_max["min_nz"])/3.0;
    max_rho2 = Rcpp::as<double>(min_and_max["max"])/3.0;
  }
  
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
  
  Eigen::SparseMatrix<double> D(p, p);
  if (!point_referenced) {
    // Adjacency matrix processing
    Eigen::VectorXd rowSums = Eigen::VectorXd::Zero(p);
    // Compute row sums of W
    for (int k = 0; k < W.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(W, k); it; ++it) {
        rowSums(it.row()) += it.value();
      }
    }
    // Set diagonal of D
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(p);
    for (int i = 0; i < p; ++i) {
      triplets.emplace_back(i, i, rowSums(i));
    }
    D.setFromTriplets(triplets.begin(), triplets.end());
  }
  
  
  // ST indicator
  Eigen::VectorXd ST = Eigen::VectorXd::Ones(ncovx);
  
  // --- Initialization ---
  bool restart = !out_prev_nullable.isNull();
  
  // MCMC parameters
  Eigen::VectorXd rho1_space_draw = Eigen::VectorXd::Constant(ncovx, 0.01);
  Eigen::VectorXd rho1_space_time_draw = Eigen::VectorXd::Constant(ncovx, 0.01);
  Eigen::VectorXd rho2_space_draw = Eigen::VectorXd::Constant(ncovx, std::max(0.5*max_rho2, min_rho2));
  Eigen::VectorXd rho2_space_time_draw = Eigen::VectorXd::Constant(ncovx, std::max(0.5*max_rho2, min_rho2));
  Eigen::VectorXd Q1invdraw_time = Eigen::VectorXd::Constant(ncovx, 10.0); // Precision for T-beta
  double s2_err_mis_draw = 0.01; // Measurement error variance
  
  // Effects: each slice is a matrix
  Eigen::VectorXd Bdraw_vec;
  std::vector<Eigen::MatrixXd> Bdraw(ncovx, Eigen::MatrixXd::Zero(1, 1));       // Constant effect
  std::vector<Eigen::MatrixXd> Btimedraw(ncovx, Eigen::MatrixXd::Zero(1, t));   // Time effect T-beta
  std::vector<Eigen::MatrixXd> Bspacedraw(ncovx, Eigen::MatrixXd::Zero(p, 1));  // Space effect S-beta
  std::vector<Eigen::MatrixXd> Bspacetimedraw(ncovx, Eigen::MatrixXd::Zero(p, t)); // Space-time effect ST-beta
  
  // Z coefficients and covariate contributions
  Eigen::VectorXd gamma_draw; // Will be resized later if ncovz > 0
  Eigen::MatrixXd meanZ = Eigen::MatrixXd::Zero(p, t); // Contribution of Z covariates
  Eigen::MatrixXd offset = Eigen::MatrixXd::Zero(p, t); // Offset term (unused in Gibbs steps)
  
  Eigen::VectorXd eta_tilde_vec = Eigen::Map<const Eigen::VectorXd>(eta_tilde.data(), eta_tilde.size());
  Eigen::VectorXd meanZ_vec = Eigen::Map<const Eigen::VectorXd>(meanZ.data(), meanZ.size());
  if (restart) {
    Rcpp::Rcout << "Restarting MCMC from previous state..." << std::endl;
    Rcpp::List out_prev = Rcpp::as<Rcpp::List>(out_prev_nullable);
    
    // Load Btimedraw, Bspacedraw, Bspacetimedraw
    {
      Rcpp::NumericVector Btime_array = out_prev["Btimedraw"];
      Rcpp::NumericVector Bspace_array = out_prev["Bspacedraw"];
      Rcpp::NumericVector Bst_array = out_prev["Bspacetimedraw"];
      
      Rcpp::IntegerVector Btime_dims = Btime_array.attr("dim");
      Rcpp::IntegerVector Bspace_dims = Bspace_array.attr("dim");
      Rcpp::IntegerVector Bst_dims = Bst_array.attr("dim");
      
      for (int k = 0; k < ncovx; ++k) {
        Btimedraw[k] = Eigen::Map<const Eigen::MatrixXd>(
          &Btime_array[0] + k * Btime_dims[0] * Btime_dims[1],
                                                          Btime_dims[0], Btime_dims[1]);
        
        Bspacedraw[k] = Eigen::Map<const Eigen::MatrixXd>(
          &Bspace_array[0] + k * Bspace_dims[0] * Bspace_dims[1],
                                                             Bspace_dims[0], Bspace_dims[1]);
        
        Bspacetimedraw[k] = Eigen::Map<const Eigen::MatrixXd>(
          &Bst_array[0] + k * Bst_dims[0] * Bst_dims[1],
                                                    Bst_dims[0], Bst_dims[1]);
      }
    }
    
    // Load rho and Q1inv vectors
    rho1_space_draw = extract_last_col(out_prev["RHO1_space_"]);
    rho1_space_time_draw = extract_last_col(out_prev["RHO1_space_time_"]);
    rho2_space_draw = extract_last_col(out_prev["RHO2_space_"]);
    rho2_space_time_draw = extract_last_col(out_prev["RHO2_space_time_"]);
    Q1invdraw_time = extract_last_col(out_prev["Q1inv_time"]);
    
    // Load s2_err_mis_draw
    {
      Rcpp::NumericMatrix s2_mat = out_prev["S2_err_mis_"];
      s2_err_mis_draw = s2_mat(0, s2_mat.ncol() - 1);
    }
    
    // Load eta_tilde
    {
      Rcpp::NumericMatrix eta_mat = out_prev["eta_tilde"];
      eta_tilde = Rcpp::as<Eigen::MatrixXd>(eta_mat);
      eta_tilde_vec = Eigen::Map<const Eigen::VectorXd>(eta_tilde.data(), eta_tilde.size());
    }
    
    // Load gamma_draw and meanZ if ncovz > 0
    if (ncovz > 0) {
      Rcpp::NumericMatrix G_mat = out_prev["G_"];
      gamma_draw = extract_last_col(G_mat);
      
      Eigen::VectorXd Zgamma = Z_mat * gamma_draw;
      meanZ = Eigen::Map<const Eigen::MatrixXd>(Zgamma.data(), p, t);
      meanZ_vec = Zgamma;
    }
    
    // Load Bdraw
    {
      Rcpp::NumericVector Bdraw_array = out_prev["Bdraw"];
      Rcpp::IntegerVector Bdraw_dims = Bdraw_array.attr("dim");
	  Bdraw_vec.resize(ncovx);
      for (int k = 0; k < ncovx; ++k) {
        Bdraw[k] = Eigen::Map<const Eigen::MatrixXd>(
          &Bdraw_array[0] + k * Bdraw_dims[0] * Bdraw_dims[1],
                                                          Bdraw_dims[0], Bdraw_dims[1]);
		Bdraw_vec(k) = Bdraw[k](0,0);
      }
    }
    
  } else { // Initial LM fit if not restarting
    Eigen::MatrixXd regressors;
    
    // Reshape X into (p*t) x ncovx
    Eigen::MatrixXd X_reshaped(p * t, ncovx);
    for (int k = 0; k < ncovx; ++k) {
      Eigen::Map<const Eigen::VectorXd> x_vec(X[k].data(), X[k].size());
      X_reshaped.col(k) = x_vec;
    }
    
    if (ncovz > 0) {
      // Combine Z and X into regressors
      regressors.resize(p * t, ncovz + ncovx);
      regressors << Z_mat, X_reshaped;
    } else {
      regressors = X_reshaped;
    }
    
    // Fit linear model
    Eigen::VectorXd bml = lmfit(regressors, eta_tilde_vec);
    if (bml.size() != ncovz + ncovx) {
      Rcpp::warning("LM fit returned unexpected number of coefficients. Initializing Bdraw to zero.");
      Bdraw_vec = Eigen::VectorXd::Zero(ncovx);
    } else {
      if (ncovz > 0) {
        gamma_draw = bml.head(ncovz);
        Eigen::VectorXd Zgamma = Z_mat * gamma_draw;
        meanZ = Eigen::Map<const Eigen::MatrixXd>(Zgamma.data(), p, t);
        meanZ_vec = Zgamma;
      }
      Bdraw_vec = bml.tail(ncovx);
    }
    for (int k = 0; k < ncovx; ++k) {
      Bdraw[k](0,0) = Bdraw_vec(k);
    }
  }
  
  
  
  
  
  // Dimensions
  const int q = 1;
  const int q2 = p;
  const int Tq = t * q;
  const int Tq2 = t * q2;
  
  
  
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
  
  // Construct sparse difference matrix H2 (Tq2 x Tq2)
  Eigen::SparseMatrix<double> H2(Tq2, Tq2);
  std::vector<Eigen::Triplet<double>> triplets_H2;
  for (int i = 0; i < Tq2; ++i) {
    triplets_H2.emplace_back(i, i, 1.0);
    if (i >= q2) {
      triplets_H2.emplace_back(i, i - q2, -1.0);
    }
  }
  H2.setFromTriplets(triplets_H2.begin(), triplets_H2.end());
  
  
  std::vector<Eigen::SparseMatrix<double>> bigG(ncovx);   // pt x t
  std::vector<Eigen::SparseMatrix<double>> bigG2(ncovx);  // pt x pt
  std::vector<Eigen::SparseMatrix<double>> bigG3(ncovx);  // pt x p
  
  for (int k = 0; k < ncovx; ++k) {
    std::vector<Eigen::Triplet<double>> triplets_G;
    std::vector<Eigen::Triplet<double>> triplets_G2;
    std::vector<Eigen::Triplet<double>> triplets_G3;
    
    for (int i = 0; i < t; ++i) {
      const Eigen::VectorXd& x_col = X[k].col(i);  // p x 1
      
      // bigG: pt x t (each column is stacked x_col)
      for (int j = 0; j < p; ++j) {
        triplets_G.emplace_back(i * p + j, i, x_col(j));
      }
      
      // bigG2: pt x pt (block diagonal of diag(x_col))
      for (int j = 0; j < p; ++j) {
        int idx = i * p + j;
        triplets_G2.emplace_back(idx, idx, x_col(j));
      }
      
      // bigG3: pt x p (each block is diag(x_col))
      for (int j = 0; j < p; ++j) {
        triplets_G3.emplace_back(i * p + j, j, x_col(j));
      }
    }
    
    Eigen::SparseMatrix<double> G(p * t, t);
    G.setFromTriplets(triplets_G.begin(), triplets_G.end());
    bigG[k] = G;
    
    Eigen::SparseMatrix<double> G2(p * t, p * t);
    G2.setFromTriplets(triplets_G2.begin(), triplets_G2.end());
    bigG2[k] = G2;
    
    Eigen::SparseMatrix<double> G3(p * t, p);
    G3.setFromTriplets(triplets_G3.begin(), triplets_G3.end());
    bigG3[k] = G3;
  }
  
  std::vector<Eigen::SparseMatrix<double>> Q1invdraw_space(ncovx);
  std::vector<Eigen::SparseMatrix<double>> Q1invdraw_spacetime(ncovx);
  if (!point_referenced) {
    
    for (int k = 0; k < ncovx; ++k) {
      if (rho1_space_draw(k) == 0.0) {
        Rcpp::stop("Initial rho1_space_draw is zero.");
      }
      if (rho1_space_time_draw(k) == 0.0) {
        Rcpp::stop("Initial rho1_space_time_draw is zero.");
      }
      
      Eigen::SparseMatrix<double> Q_space = D - rho2_space_draw(k) * W;
      Q1invdraw_space[k] = (1.0 / rho1_space_draw(k)) * Q_space;
      
      Eigen::SparseMatrix<double> Q_spacetime = D - rho2_space_time_draw(k) * W;
      Q1invdraw_spacetime[k] = (1.0 / rho1_space_time_draw(k)) * Q_spacetime;
    }
  } else {
    for (int k = 0; k < ncovx; ++k) {
      if (rho1_space_draw(k) == 0.0) {
        Rcpp::stop("Initial rho1_space_draw is zero.");
      }
      if (rho1_space_time_draw(k) == 0.0) {
        Rcpp::stop("Initial rho1_space_time_draw is zero.");
      }
      
      Eigen::SparseMatrix<double> Q_space = (-W_dense.array() * (1.0 / rho2_space_draw(k))).exp().matrix().inverse().sparseView();
      Q1invdraw_space[k] = (1.0 / rho1_space_draw(k)) * Q_space;
      
      Eigen::SparseMatrix<double> Q_spacetime = (-W_dense.array() * (1.0 / rho2_space_time_draw(k))).exp().matrix().inverse().sparseView();
      Q1invdraw_spacetime[k] = (1.0 / rho1_space_time_draw(k)) * Q_spacetime;
    }
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
  Eigen::MatrixXd RHO1_space_(ncovx, collections);
  Eigen::MatrixXd RHO1_space_time_(ncovx, collections);
  Eigen::MatrixXd RHO2_space_(ncovx, collections);
  Eigen::MatrixXd RHO2_space_time_(ncovx, collections);
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
  
  std::vector<Eigen::MatrixXd> B_postmean(ncovx, Eigen::MatrixXd::Zero(1, 1));
  std::vector<Eigen::MatrixXd> B2_postmean(ncovx, Eigen::MatrixXd::Zero(1, 1));
  std::vector<Eigen::MatrixXd> Btime_postmean(ncovx, Eigen::MatrixXd::Zero(1, t));
  std::vector<Eigen::MatrixXd> Btime2_postmean(ncovx, Eigen::MatrixXd::Zero(1, t));
  std::vector<Eigen::MatrixXd> Bspace_postmean(ncovx, Eigen::MatrixXd::Zero(p, 1));
  std::vector<Eigen::MatrixXd> Bspace2_postmean(ncovx, Eigen::MatrixXd::Zero(p, 1));
  std::vector<Eigen::MatrixXd> Bspacetime_postmean(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> Bspacetime2_postmean(ncovx, Eigen::MatrixXd::Zero(p, t));
  
  std::vector<Eigen::MatrixXd> B2_c_t(ncovx, Eigen::MatrixXd::Zero(1, t));
  std::vector<Eigen::MatrixXd> B2_c_s(ncovx, Eigen::MatrixXd::Zero(p, 1));
  std::vector<Eigen::MatrixXd> B2_c_t_s_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_c_t_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_t_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_t_s(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_t_s_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_c_t_s(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_s_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> B2_c_s_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  
  std::vector<Eigen::MatrixXd> E_t_s(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> E_t_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  std::vector<Eigen::MatrixXd> E_s_st(ncovx, Eigen::MatrixXd::Zero(p, t));
  
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
      y2_vec -= bigG2[k] * Eigen::Map<Eigen::VectorXd>(Bspacetimedraw[k].data(), Bspacetimedraw[k].size());
      y2_vec -= bigG3[k] * Bspacedraw[k].col(0);
	  Eigen::MatrixXd xxxx = X[k] * Bdraw_vec(k);
      y2_vec -= Eigen::Map<Eigen::VectorXd>(xxxx.data(), xxxx.size());
    }
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
      
      // Reshape and center
      Eigen::MatrixXd Bdrawc = Eigen::Map<Eigen::MatrixXd>(bb.data(), 1, t);
      Btimedraw[k] = Bdrawc.row(0).array() - Bdrawc.row(0).mean();
      
    }
    
    
    
    // Step I.b: Sampling ST-beta (space-time effects)
    // Compute residual y_k
    y2_vec = eta_tilde_vec - meanZ_vec;
    for (int k = 0; k < ncovx; ++k) {
      y2_vec -= bigG[k] * Eigen::Map<Eigen::VectorXd>(Btimedraw[k].row(0).data(), Btimedraw[k].cols());
      y2_vec -= bigG3[k] * Bspacedraw[k].col(0);
      Eigen::MatrixXd xxxx = X[k] * Bdraw_vec(k);
      y2_vec -= Eigen::Map<Eigen::VectorXd>(xxxx.data(), xxxx.size());
    }
    
    for (int k = 0; k < ncovx; ++k) {
      Eigen::VectorXd y_k = y2_vec;
      for (int k2 = 0; k2 < ncovx; ++k2) {
        if (k == k2) continue;
        y_k -= bigG2[k2] * Eigen::Map<const Eigen::VectorXd>(Bspacetimedraw[k2].data(), p * t);
      }
      
      // Construct block-diagonal prior precision matrix invS_diag
      Eigen::SparseMatrix<double> invS_diag(p * t, p * t);
      std::vector<Eigen::Triplet<double>> triplets;
      for (int i = 0; i < p; ++i) {
        triplets.emplace_back(i, i, 1.0 / V_beta_0);  // prior for beta_0
      }
      for (int i = 1; i < t; ++i) {
        for (int j = 0; j < Q1invdraw_spacetime[k].outerSize(); ++j) {
          for (Eigen::SparseMatrix<double>::InnerIterator it(Q1invdraw_spacetime[k], j); it; ++it) {
            int row = i * p + it.row();
            int col = i * p + it.col();
            triplets.emplace_back(row, col, it.value());
          }
        }
      }
      invS_diag.setFromTriplets(triplets.begin(), triplets.end());
      
      // Total prior precision
      Eigen::SparseMatrix<double> K = H2.transpose() * invS_diag * H2;
      
      // Observation precision
      Eigen::SparseMatrix<double> GinvOmega11 = bigG2[k].transpose() * (1.0 / s2_err_mis_draw);
      Eigen::SparseMatrix<double> GinvOmega11G = GinvOmega11 * bigG2[k];
      
      // Posterior precision
      Eigen::SparseMatrix<double> invP = K + GinvOmega11G;
      invP = 0.5 * (invP + Eigen::SparseMatrix<double>(invP.transpose()));  // ensure symmetry
      
      // Posterior mean
      Eigen::VectorXd tmp = GinvOmega11 * y_k;
      
      // Sample from posterior
      Eigen::MatrixXd bb = rmvnorm_eigen_sparse(1, invP, tmp);  // pt x 1
      // Reshape and apply sum-to-zero constraints
      Eigen::MatrixXd Bdrawc = Eigen::Map<Eigen::MatrixXd>(bb.data(), p, t);
      Bdrawc = Bdrawc.colwise() - Bdrawc.rowwise().mean();  // center columns (space)
      Bdrawc = Bdrawc.rowwise() - Bdrawc.colwise().mean();  // center rows (time)
      Bspacetimedraw[k] = Bdrawc * ST(k);  // apply ST indicator
    }
    
    // Step I.c: Sampling S-beta (spatial effects)
    y2_vec = eta_tilde_vec - meanZ_vec;
    for (int k = 0; k < ncovx; ++k) {
      y2_vec -= bigG[k] * Eigen::Map<Eigen::VectorXd>(Btimedraw[k].row(0).data(), Btimedraw[k].cols());
      y2_vec -= bigG2[k] * Eigen::Map<Eigen::VectorXd>(Bspacetimedraw[k].data(), Bspacetimedraw[k].size());
      Eigen::MatrixXd xxxx = X[k] * Bdraw_vec(k);
      y2_vec -= Eigen::Map<Eigen::VectorXd>(xxxx.data(), xxxx.size());
    }
    
    for (int k = 0; k < ncovx; ++k) {
      // Compute residual y_k
      Eigen::VectorXd y_k = y2_vec;
      for (int k2 = 0; k2 < ncovx; ++k2) {
        if (k == k2) continue;
        y_k -= bigG3[k2] * Bspacedraw[k2];
      }
      
      // Prior precision
      const Eigen::SparseMatrix<double>& K = Q1invdraw_space[k];
      
      // Observation precision
      Eigen::SparseMatrix<double> GinvOmega11 = bigG3[k].transpose() * (1.0 / s2_err_mis_draw);
      Eigen::SparseMatrix<double> GinvOmega11G = GinvOmega11 * bigG3[k];
      
      // Posterior precision
      Eigen::SparseMatrix<double> invP = K + GinvOmega11G;
      invP = 0.5 * (invP + Eigen::SparseMatrix<double>(invP.transpose()));  // ensure symmetry
      
      // Posterior mean
      Eigen::VectorXd tmp = GinvOmega11 * y_k;
      
      // Sample from posterior
      Eigen::MatrixXd bb = rmvnorm_eigen_sparse(1, invP, tmp);  // p x 1
      
      // Center (sum-to-zero constraint)
      Eigen::VectorXd centered = bb.col(0).array() - bb.col(0).mean();
      Bspacedraw[k] = centered;
    }
    
    
    // Step I.d: Sampling constant coefficients (C-beta and Z-gamma)
    // Compute residual y2 = eta_tilde - T*beta_t - S*beta_s - ST*beta_st
    Eigen::VectorXd y2_reg = eta_tilde_vec;
    for (int k = 0; k < ncovx; ++k) {
      y2_reg -= bigG[k] * Btimedraw[k].transpose();
      y2_reg -= bigG2[k] * Eigen::Map<const Eigen::VectorXd>(Bspacetimedraw[k].data(), p * t);
      y2_reg -= bigG3[k] * Bspacedraw[k];
    }
    
    // Construct combined regressor matrix [Z, X]
    Eigen::MatrixXd X_reshaped(p * t, ncovx);
    for (int k = 0; k < ncovx; ++k) {
      X_reshaped.col(k) = Eigen::Map<const Eigen::VectorXd>(X[k].data(), p * t);
    }
    
    Eigen::MatrixXd X_reg;
    if (ncovz > 0) {
      X_reg.resize(p * t, ncovz + ncovx);
      X_reg << Z_mat, X_reshaped;
    } else {
      X_reg = X_reshaped;
    }
    
    // Prior precision (diagonal, small values for diffuse prior)
    Eigen::MatrixXd K_reg = Eigen::MatrixXd::Identity(ncovz + ncovx, ncovz + ncovx) * (1.0 / V_gamma);
    
    // Posterior precision and mean
    Eigen::MatrixXd GinvOmega11_reg = X_reg.transpose() * (1.0 / s2_err_mis_draw);
    Eigen::MatrixXd GinvOmega11G_reg = GinvOmega11_reg * X_reg;
    Eigen::MatrixXd invP_reg = K_reg + GinvOmega11G_reg;
    invP_reg = 0.5 * (invP_reg + invP_reg.transpose());  // Ensure symmetry
    
    Eigen::VectorXd tmp_reg = GinvOmega11_reg * y2_reg;
    
    // Sample from posterior
	Eigen::MatrixXd P_reg = invP_reg.inverse();
	Eigen::LLT<Eigen::MatrixXd> llt(P_reg);
	Eigen::MatrixXd L_reg = llt.matrixL();
    Eigen::VectorXd gamma_beta_draw = P_reg * tmp_reg;
	Eigen::VectorXd normal_std_reg = Eigen::VectorXd::Zero(gamma_beta_draw.size());
	for (int i = 0; i < gamma_beta_draw.size(); ++i) {
      normal_std_reg(i) += R::rnorm(0.0, 1.0);
    }
	gamma_beta_draw += L_reg * normal_std_reg;
	
    
    
    // Update gamma_draw and Bdraw
    if (ncovz > 0) {
      gamma_draw = gamma_beta_draw.head(ncovz);
      Eigen::VectorXd Zgamma = Z_mat * gamma_draw;
      meanZ = Eigen::Map<const Eigen::MatrixXd>(Zgamma.data(), p, t);
      meanZ_vec = Zgamma;
    }
    Eigen::VectorXd Bdraw_vec = gamma_beta_draw.tail(ncovx);
    for (int k = 0; k < ncovx; ++k) {
      Bdraw[k](0,0) = Bdraw_vec(k);
    }
    
    
    // Step I.e: Computing the current mean component from X
    Eigen::MatrixXd current_meanY1 = Eigen::MatrixXd::Zero(p, t);
    
    for (int k = 0; k < ncovx; ++k) {
      Eigen::MatrixXd B_c_full = Eigen::MatrixXd::Constant(p, t, Bdraw_vec(k));
      Eigen::MatrixXd B_t_full = Btimedraw[k].colwise().replicate(p);  // 1 x t → p x t
      Eigen::MatrixXd B_s_full = Bspacedraw[k].rowwise().replicate(t); // p x 1 → p x t
      const Eigen::MatrixXd& B_st_full = Bspacetimedraw[k];
      Eigen::ArrayXXd appoX_B = X[k].array() * (B_c_full + B_t_full + B_s_full + B_st_full).array();
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
      
      // --- SPACE-TIME (ST-beta) ---
      Eigen::MatrixXd Bst_k = Bspacetimedraw[k];  // p x t
      Eigen::MatrixXd Btemp = Eigen::MatrixXd::Zero(p, t);
      Btemp.block(0, 1, p, t - 1) = Bst_k.block(0, 1, p, t - 1) - Bst_k.block(0, 0, p, t - 1);
      
      if (ST(k) == 1.0) {
        if (!point_referenced) {
          rho2_space_time_draw(k) = MH_spatial_correlation_CAR(W, Btemp, rho1_space_time_draw(k), 0.1, 1.0);
          int pstar = (rho2_space_time_draw(k) == 1.0) ? (p - 1) : p;
          Eigen::SparseMatrix<double> Q_inv_st = D - rho2_space_time_draw(k) * W;
          rho1_space_time_draw(k) = posterior_conditional_variance(Btemp, Q_inv_st, 0.01, 0.01, pstar, t - 1);
          Q1invdraw_spacetime[k] = (1.0 / rho1_space_time_draw(k)) * Q_inv_st;
        } else {
          rho2_space_time_draw(k) = MH_spatial_correlation_EXP(W_dense, Btemp, rho1_space_time_draw(k), min_rho2, max_rho2);
          Eigen::SparseMatrix<double> Q_inv_st = (-W_dense.array() * (1.0 / rho2_space_time_draw(k))).exp().matrix().inverse().sparseView();
          rho1_space_time_draw(k) = posterior_conditional_variance(Btemp, Q_inv_st, 0.01, 0.01, p, t - 1);
          Q1invdraw_spacetime[k] = (1.0 / rho1_space_time_draw(k)) * Q_inv_st;
        }
      } else {
        rho1_space_time_draw(k) = 1e-9;
        Q1invdraw_spacetime[k].resize(p, p);
        Q1invdraw_spacetime[k].setZero();
      }
      
      // --- SPACE (S-beta) ---
      Eigen::MatrixXd Bs_k = Bspacedraw[k];  // p x 1
      if (!point_referenced) {
        rho2_space_draw(k) = MH_spatial_correlation_CAR(W, Bs_k, rho1_space_draw(k), 0.1, 1.0);
        int pstar_s = (rho2_space_draw(k) == 1.0) ? (p - 1) : p;
        Eigen::SparseMatrix<double> Q_inv_s = D - rho2_space_draw(k) * W;
        rho1_space_draw(k) = posterior_conditional_variance(Bs_k, Q_inv_s, 0.01, 0.01, pstar_s, 1);
        Q1invdraw_space[k] = (1.0 / rho1_space_draw(k)) * Q_inv_s;
      } else {
        rho2_space_draw(k) = MH_spatial_correlation_EXP(W_dense, Bs_k, rho1_space_draw(k), min_rho2, max_rho2);
        Eigen::SparseMatrix<double> Q_inv_s = (-W_dense.array() * (1.0 / rho2_space_draw(k))).exp().matrix().inverse().sparseView();
        rho1_space_draw(k) = posterior_conditional_variance(Bs_k, Q_inv_s, 0.01, 0.01, p, 1);
        Q1invdraw_space[k] = (1.0 / rho1_space_draw(k)) * Q_inv_s;
      }
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
        RHO1_space_(k, collect_count) = rho1_space_draw(k);
        RHO1_space_time_(k, collect_count) = rho1_space_time_draw(k);
        RHO2_space_(k, collect_count) = rho2_space_draw(k);
        RHO2_space_time_(k, collect_count) = rho2_space_time_draw(k);
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
      double term3 = -0.5 * static_cast<double>(t) * log_det_Sigma;
      store_llike(collect_count) = term1 + term2 + term3;*/
			// Calculate point-wise log-likelihood
			Eigen::MatrixXd term1 = -0.5 * yhat.array().square() / s2_err_mis_draw;
			Eigen::MatrixXd term2 = Eigen::MatrixXd::Constant(p, t, -0.5 * std::log(2.0 * M_PI));
			Eigen::MatrixXd term3 = Eigen::MatrixXd::Constant(p, t, -0.5 * std::log(s2_err_mis_draw));
			store_llike[collect_count] = term1 + term2 + term3;

      Eigen::MatrixXd Ypred2 = Eigen::MatrixXd::NullaryExpr(p, t, norm_mcmc);
			Ypred2 += thetay_draw;
      if (keepYpred) {
        YPRED_out[collect_count] = Ypred;
      }
      
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
        B_postmean[k] += Bdraw[k];
        B2_postmean[k] += Bdraw[k].array().square().matrix();
        
        Btime_postmean[k] += Btimedraw[k];
        Btime2_postmean[k] += Btimedraw[k].array().square().matrix();
        
        Bspace_postmean[k] += Bspacedraw[k];
        Bspace2_postmean[k] += Bspacedraw[k].array().square().matrix();
        
        Bspacetime_postmean[k] += Bspacetimedraw[k];
        Bspacetime2_postmean[k] += Bspacetimedraw[k].array().square().matrix();
        
        double Bdrawk = Bdraw_vec(k);
        // B_c_t: 1 x t
        Eigen::MatrixXd B_c_t = Eigen::MatrixXd::Ones(1, t) * Bdrawk + Btimedraw[k];
        B2_c_t[k] += B_c_t.array().square().matrix();
        
        // B_c_s: p x 1
        Eigen::MatrixXd B_c_s = Eigen::MatrixXd::Ones(p, 1) * Bdrawk + Bspacedraw[k];
        B2_c_s[k] += B_c_s.array().square().matrix();
        
        // B_c_t_s_st: p x t
        Eigen::MatrixXd B_c_t_s_st = Eigen::MatrixXd::Ones(p, t) * Bdrawk
					+ Btimedraw[k].colwise().replicate(p)
          + Bspacedraw[k].rowwise().replicate(t)
          + Bspacetimedraw[k];
        B2_c_t_s_st[k] += B_c_t_s_st.array().square().matrix();
          
				// B_c_t_st: p x t
				Eigen::MatrixXd B_c_t_st = Eigen::MatrixXd::Ones(p, t) * Bdraw[k] + Btimedraw[k].colwise().replicate(p)
				+ Bspacetimedraw[k];
				B2_c_t_st[k] += B_c_t_st.array().square().matrix();
						
				// B_t_st: p x t
				Eigen::MatrixXd B_t_st = Btimedraw[k].colwise().replicate(p) + Bspacetimedraw[k];
				B2_t_st[k] += B_t_st.array().square().matrix();
          
				// B_t_s: p x t
				Eigen::MatrixXd B_t_s = Btimedraw[k].colwise().replicate(p) + Bspacedraw[k].rowwise().replicate(t);
				B2_t_s[k] += B_t_s.array().square().matrix();
          
				// B_t_s_st: p x t
				Eigen::MatrixXd B_t_s_st = Btimedraw[k].colwise().replicate(p)
					+ Bspacedraw[k].rowwise().replicate(t)
					+ Bspacetimedraw[k];
				B2_t_s_st[k] += B_t_s_st.array().square().matrix();
					
				// B_c_t_s: p x t
				Eigen::MatrixXd B_c_t_s = Eigen::MatrixXd::Ones(p, t) * Bdrawk
					+ Btimedraw[k].colwise().replicate(p)
					+ Bspacedraw[k].rowwise().replicate(t);
				B2_c_t_s[k] += B_c_t_s.array().square().matrix();
					
				// B_s_st: p x t
				Eigen::MatrixXd B_s_st = Bspacedraw[k].rowwise().replicate(t) + Bspacetimedraw[k];
				B2_s_st[k] += B_s_st.array().square().matrix();
					
				// B_c_s_st: p x t
				Eigen::MatrixXd B_c_s_st = Eigen::MatrixXd::Ones(p, t) * Bdrawk
				+ Bspacedraw[k].rowwise().replicate(t)
					+ Bspacetimedraw[k];
				B2_c_s_st[k] += B_c_s_st.array().square().matrix();
					
				// Interactions
				E_t_s[k] += (Btimedraw[k].colwise().replicate(p).array() * Bspacedraw[k].rowwise().replicate(t).array()).matrix();
				E_t_st[k] += (Btimedraw[k].colwise().replicate(p).array() * Bspacetimedraw[k].array()).matrix();
				E_s_st[k] += (Bspacedraw[k].rowwise().replicate(t).array() * Bspacetimedraw[k].array()).matrix();
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
  ave_results["B_postmean"] = cube_to_array(B_postmean, n_samples_collected);
  ave_results["B2_postmean"] = cube_to_array(B2_postmean, n_samples_collected);
  ave_results["Btime_postmean"] = cube_to_array(Btime_postmean, n_samples_collected);
  ave_results["Btime2_postmean"] = cube_to_array(Btime2_postmean, n_samples_collected);
  ave_results["Bspace_postmean"] = cube_to_array(Bspace_postmean, n_samples_collected);
  ave_results["Bspace2_postmean"] = cube_to_array(Bspace2_postmean, n_samples_collected);
  ave_results["Bspacetime_postmean"] = cube_to_array(Bspacetime_postmean, n_samples_collected);
  ave_results["Bspacetime2_postmean"] = cube_to_array(Bspacetime2_postmean, n_samples_collected);
  ave_results["B2_c_t_s_st"] = cube_to_array(B2_c_t_s_st, n_samples_collected);
  
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
  
  double term1_hat = -0.5 * sse_hat / s2_err_hat;
  double term2_hat = -0.5 * static_cast<double>(p * t) * std::log(2.0 * M_PI);
  double term3_hat = -0.5 * static_cast<double>(p * t) * std::log(s2_err_hat);
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
  out_results["RHO1_space_"] = RHO1_space_;
  out_results["RHO1_space_time_"] = RHO1_space_time_;
  out_results["RHO2_space_"] = RHO2_space_;
  out_results["RHO2_space_time_"] = RHO2_space_time_;
  out_results["Q1inv_time"] = Q1inv_time_out;
  if (ncovz > 0) {
    out_results["G_"] = G_out;
  }
  out_results["store_llike"] = cube_to_array(store_llike, 1.0);
  if (keepYpred) {
    out_results["YPRED_"] = cube_to_array(YPRED_out, 1.0);
  }
  out_results["RMSE_"] = RMSE_;
  out_results["MAE_"] = MAE_;
  out_results["chi_sq_pred_"] = chi_sq_pred_;
  out_results["chi_sq_obs_"] = chi_sq_obs_;
  
  // Store final state for potential restart
  out_results["eta_tilde"] = eta_tilde;
  out_results["Bdraw"] = cube_to_array(Bdraw, 1.0);
  out_results["Btimedraw"] = cube_to_array(Btimedraw, 1.0);
  out_results["Bspacedraw"] = cube_to_array(Bspacedraw, 1.0);
  out_results["Bspacetimedraw"] = cube_to_array(Bspacetimedraw, 1.0);
  
  
  
  // Return both lists
  return Rcpp::List::create(Rcpp::Named("ave") = ave_results,
                            Rcpp::Named("out") = out_results);
  
}

