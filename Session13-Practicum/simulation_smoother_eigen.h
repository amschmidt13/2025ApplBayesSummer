#ifndef SIMULATION_SMOOTHER
#define SIMULATION_SMOOTHER


#include <RcppEigen.h>

Eigen::MatrixXd rmvnorm_eigen_sparse(int n, const Eigen::SparseMatrix<double>& precision, const Eigen::VectorXd& location) {
  const int T = precision.rows();

  // Generate standard normal noise
  auto norm = [&] (double) {return R::rnorm(0.0, 1.0);};
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::NullaryExpr(T, n,  norm);

  // Repeat location vector across columns
  Eigen::MatrixXd location_matrix = location.replicate(1, n);

  // Cholesky decomposition of precision matrix (upper triangular)
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper, Eigen::NaturalOrdering<int>> chol_solver;
  chol_solver.compute(precision);
  if (chol_solver.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed.");
  }

  // Solve for mean-adjusted samples
  Eigen::MatrixXd U_mu_mat = chol_solver.matrixU().transpose().triangularView<Eigen::Lower>().solve(location_matrix);
  Eigen::MatrixXd draw = chol_solver.matrixU().triangularView<Eigen::Upper>().solve(U_mu_mat + epsilon);

  return draw;
}


#endif