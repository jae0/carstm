
posterior_summary = function(X, pivot = NULL, id=TRUE) {
    if (is.vector(X)) X = t(X)
    
    nd = length(dim(X))
    if (is.null(pivot)) pivot = 1:(nd-1)
    if (nd==2) {
        out =  data.frame( cbind( 
            mean = apply( X, pivot, mean, na.rm=TRUE ),
            sd = apply( X, pivot, sd, na.rm=TRUE ),
            quant0.025 = apply(X, pivot, quantile, probs=0.025, na.rm=TRUE),
            quant0.5   = apply(X, pivot, median, na.rm=TRUE),
            quant0.975 = apply(X, pivot, quantile, probs=0.975, na.rm=TRUE)
        ))
        row.names(out) = row.names(X)

    } else {
        out =  data.frame( abind( 
            mean = apply( X, pivot, mean, na.rm=TRUE ),
            sd = apply( X, pivot, sd, na.rm=TRUE ),
            quant0.025 = apply(X, pivot, quantile, probs=0.025, na.rm=TRUE),
            quant0.5   = apply( X, pivot, median, na.rm=TRUE),
            quant0.975 = apply(X, pivot, quantile, probs=0.975, na.rm=TRUE),
            hier.names=TRUE
        ))
        row.names(out) = row.names(X)
    }
    if (id) out$ID = row.names(out)
    return( out )
}
