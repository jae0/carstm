

format_results = function(X, labels=NULL ){
    out = data.table( t( apply( X, 1, c) ) )
    if (!is.null(labels)) row.names(out) = labels
    return(out)
}

