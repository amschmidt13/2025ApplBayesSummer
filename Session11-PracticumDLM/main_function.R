sourceCpp("dlm_uv.cpp") # Compile the C++ code

DLM = function(y, X, Z, nrep, nburn, thin, prior_list) {
  dim(X) = c(1, dim(X))
  if (!is.null(Z)) {dim(Z) = c(1, dim(Z))}
  V_beta_0 = prior_list$V_beta_0
  V_gamma = prior_list$V_gamma
  a1 = prior_list$a1
  b1 = prior_list$b1
  s2_a = prior_list$s2_a
  s2_b = prior_list$s2_b
  
  re = dlm_cpp(rbind(y), X, Z, nrep, nburn, thin, V_beta_0, V_gamma, a1, b1, s2_a, s2_b)
  return(re)
}