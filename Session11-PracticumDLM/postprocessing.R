plot.tvc = function(dlmod, beta_all=NULL) {
  require(ggplot2)
  ave <- dlmod$ave
  
  gg = list()
  for (k in 1:dim(ave$Btime_postmean)[3]) {
    Bmean = ave$Btime_postmean[1,,k]
    Bsd = sqrt(ave$Btime2_postmean[1,,k] - ave$Btime_postmean[1,,k]^2)
    q = qnorm(.975)
    df = data.frame(k=k, Time=1:length(Bmean), Bmean=Bmean,
                    B025 = Bmean - q*Bsd, B975 = Bmean + q*Bsd)
    if (! is.null(beta_all)) df$Btrue = beta_all[1,,k]
    gg[[k]] = ggplot(df, aes(x=Time)) +
      geom_line(aes(y= Bmean)) +
      geom_line(aes(y= B025), linetype='dashed') +
      geom_line(aes(y= B975), linetype='dashed') +
      ylab(paste('beta', k-1))
    if (! is.null(beta_all)) 
      gg[[k]] = gg[[k]] + geom_line(aes(y=Btrue), col='red')
  }
  return(gg)
}


plot.fitted = function(dlmod, Y) {
  require(ggplot2)
  ave <- dlmod$ave
  
  df = data.frame(Time = 1:length(Y),
                  Y = Y,
                  Fitted = ave$Ypred_mean[1,])
  df = pivot_longer(df, Y:Fitted, names_to = "Outcome", values_to = 'Value')
  
  gg = ggplot(df, aes(x=Time)) +
    geom_line(aes(y=Value, col=Outcome)) +
    ggtitle('Observed vs. Fitted')
  return(gg)
}


posterior.chains = function(dlmod) {
  require(coda)
  out = dlmod$out
  nrep = length(out$S2_err_mis_)
  
  param = cbind(out$S2_err_mis_[1,])
  param = cbind(param, t(out$Q1inv_time^-1))
  nc = NCOL(param)
  cn = c('sigma2epsilon', paste0('sig2beta', 2:nc-2))
  colnames(param) = cn
  if ("G_" %in% names(out)) {
    param = cbind(param, t(out$G_))
    nr = NROW(out$G_)
    colnames(param) = c(cn, paste0('gamma', 1:nr))
  }
  
  param = mcmc(param)
  return(param)
  
}
