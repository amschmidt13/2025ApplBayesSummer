plot.tvc = function(dlmod) {
  require(ggplot2)
  require(ggpubr)
  ave <- dlmod$ave
  ncov = dim(ave$Btime_postmean)[3]
  q = qnorm(.975)
  
  gg = list()
  for (k in 1:ncov) {
    Bmean = ave$Btime_postmean[1,,k]
    Bsd = sqrt(ave$Btime2_postmean[1,,k] - ave$Btime_postmean[1,,k]^2)
    df = data.frame(k=k, Time=1:length(Bmean), Bmean=Bmean,
                    B025 = Bmean - q*Bsd, B975 = Bmean + q*Bsd)
    
    gg[[k]] = ggplot(df, aes(x=Time)) +
      geom_line(aes(y= Bmean)) +
      geom_line(aes(y= B025), linetype='dashed') +
      geom_line(aes(y= B975), linetype='dashed') +
      labs(y=paste('beta', k-1), title = 'Temporal Effect')
    
  }
  return(gg)
}

plot.svc = function(dlmod, Coo_sf, region) {
  require(ggplot2)
  require(ggpubr)
  ave <- dlmod$ave
  ncov = dim(ave$Btime_postmean)[3]
  q = qnorm(.975)
  
  gg = list()
  for (k in 1:ncov) {
    Bmean = ave$Bspace_postmean[,1,k]
    Bsd = sqrt(ave$Bspace2_postmean[,1,k] - ave$Bspace_postmean[,1,k]^2)
    Coo_sf$Bmean = Bmean
    Coo_sf$Bsd = Bsd
    Coo_sf$B025 = Bmean - q*Bsd
    Coo_sf$B975 = Bmean + q*Bsd
    
    
    gg1 = ggplot(Coo_sf) +
      geom_sf(data = region, fill='white') +
      geom_sf(aes(col=Bmean), size=4) +
      labs(title = paste('Spatial Effect of beta', k-1), color = 'Mean')
    gg2 = ggplot(Coo_sf) +
      geom_sf(data = region, fill='white') +
      geom_sf(aes(col=Bsd), size=4) +
      labs(title = paste('Spatial Effect of beta', k-1), color = 'Std. Dev')
    gg[[k]] = ggarrange(gg1, gg2, nrow = 1, ncol = 2)
  }
  return(gg)
}


plot.stvc = function(dlmod, ids) {
  require(ggplot2)
  require(ggpubr)
  ave <- dlmod$ave
  ncov = dim(ave$Bspacetime_postmean)[3]
  q = qnorm(.975)
  nn = ceiling(sqrt(length(ids)))
  
  gg = list()
  for (k in 1:ncov) {
    Bmean = ave$Bspacetime_postmean[,,k]
    Bsd = sqrt(ave$Bspacetime2_postmean[,,k] - ave$Bspacetime_postmean[,,k]^2)
    df = expand.grid(k=k, Space=1:NROW(Bmean), Time=1:NCOL(Bmean))
    df$Bmean=as.vector(Bmean)
    df$B025 = as.vector(Bmean - q*Bsd)
    df$B975 = as.vector(Bmean + q*Bsd)
    df = df %>% filter(Space %in% ids)
    
    gg[[k]] = ggplot(df, aes(x=Time)) +
      geom_line(aes(y= Bmean)) +
      geom_line(aes(y= B025), linetype='dashed') +
      geom_line(aes(y= B975), linetype='dashed') +
      labs(title=paste('Spatio-Temporal Effect of beta', k-1)) +
      facet_wrap(~Space, nn, nn, labeller = 'label_both')
    
  }
  return(gg)
}


plot.fitted = function(dlmod, Y, id) {
  require(ggplot2)
  ave <- dlmod$ave
  
  df = data.frame(Time = 1:NCOL(Y),
                  Y = Y[id,],
                  Fitted = ave$Ypred_mean[id,])
  df = pivot_longer(df, Y:Fitted, names_to = "Outcome", values_to = 'Value')
  
  gg = ggplot(df, aes(x=Time)) +
    geom_line(aes(y=Value, col=Outcome)) +
    ggtitle('Observed vs. Fitted')
  return(gg)
}


