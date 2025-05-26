sourceCpp("dlm_eigen_uv.cpp") # Compile the C++ code

DLM = function(y, X, Z, nrep, nburn, thin, print.interval, prior_list) {
  
  V_beta_0 = prior_list$V_beta_0
  V_gamma = prior_list$V_gamma
  a1 = prior_list$a1
  b1 = prior_list$b1
  s2_a = prior_list$s2_a
  s2_b = prior_list$s2_b
  
  ncovx = dim(X)[2]
  X = lapply(1:ncovx, \(i) rbind(X[,i]))
  if (!is.null(Z)) {
    ncovz = dim(Z)[2]
    Z = lapply(1:ncovz, \(i) rbind(Z[,i]))
  }
  
  re = dlm_cpp(rbind(y), X, Z, nrep, nburn, thin, print.interval,
               V_beta_0, V_gamma, a1, b1, s2_a, s2_b)
  return(re)
}