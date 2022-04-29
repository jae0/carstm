carstm_posterior_simulations_summary = function( X, D=1, probs=c(0.025, 0.975) ) {
  nd = length(dim(X))
  # D == dim to summarize along 
  if (nd==3) X = colSums( X, na.rm=TRUE )  # need to aggregate along rows (domain)

  out = data.frame( year = as.numeric(rownames(X)) )
  out$mean = apply( simplify2array(X), D, mean, na.rm=TRUE )
  out$sd = apply( simplify2array(X), D, sd, na.rm=TRUE )
  out$median = apply( simplify2array(X), D, median, na.rm=TRUE )
  out$lb025 = apply( simplify2array(X), D, quantile, probs=probs[1], na.rm=TRUE )
  out$ub975 = apply( simplify2array(X), D, quantile, probs=probs[2], na.rm=TRUE )

  return(out)
}
