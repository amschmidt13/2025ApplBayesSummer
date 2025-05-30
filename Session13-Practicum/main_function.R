sourceCpp("dlm_eigen.cpp") # Compile the C++ code

DLM.st = function(y, X, Z, D, point.referenced, nrep, nburn, thin, 
                  print.interval, prior_list, out_prev = NULL) {
  
  stopifnot("Please specify nburn>1" = nburn>1)
  stopifnot("Please specify nrep>1" = nrep>1)
  
  V_beta_0 = prior_list$V_beta_0
  V_gamma = prior_list$V_gamma
  a1 = prior_list$a1
  b1 = prior_list$b1
  s2_a = prior_list$s2_a
  s2_b = prior_list$s2_b
  
  ncovx = dim(X)[3]
  X = lapply(1:ncovx, \(i) X[,,i])
  if (!is.null(Z)) {
    ncovz = dim(Z)[3]
    Z = lapply(1:ncovz, \(i) Z[,,i])
  }
  
  re = dlm_cpp(y, X, Z, D, point.referenced, nrep, nburn, thin, print.interval,
               V_beta_0, V_gamma, 
               a1, b1, s2_a, s2_b, TRUE, out_prev)
  return(re)
}
