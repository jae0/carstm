split_vector2matrix = function( u, vnames=c("AUID", "yr_factor"), sep=":", matchto=NULL ) {
  # could be array .. alter matrix to array ..
  un = names(u)
  if(is.null(un)) un=rownames(u)
  if(is.null(un)) stop( "Data vector requires names or rownames as an attribute" )
  m = 1:length(un)
  for (j in 1:length(vnames)) m = intersect(m, grep(vnames[j], un, fixed=TRUE ) )
  w = matrix( unlist(strsplit(un[m], sep)), ncol=length(vnames), byrow=TRUE)
  for (i in 1:length(vnames)) w = gsub(vnames[i], "", w)
  colnames(w) = vnames
  rownames(w) = un[m]

  if (is.null(matchto)) return (w)

  #require matching indices
  ids = list()
  for( i in 1:length(vnames)) ids[[vnames[i]]] = match( w[,vnames[i]], as.character(matchto[[i]] ) )
  out = matrix( NA, ncol=length(matchto[[2]]), nrow=length(matchto[[1]]) )
  rownames(out) = as.character(matchto[[1]] )
  colnames(out) = as.character(matchto[[2]] )
  out[ cbind(ids[[vnames[1]]], ids[[vnames[2]]])] = u[m]

  # reshape
  mm = matrix(NA, nrow=length(matchto[[1]]), ncol=length(matchto[[2]]), dimnames=matchto )
  out_row = match( rownames(out), as.character( matchto[[1]] ) )
  out_col = match( colnames(out), as.character( matchto[[2]] ) )
  mm[ out_row, out_col  ] = out[]
  out = mm
  return(out)
}
