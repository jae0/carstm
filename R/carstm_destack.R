
carstm_destack = function( inputdata=NULL, tag="predictions", dimensionality=NULL, AUID=NULL, year=NULL, dyear=NULL ) {

  res = list( dimensionality = dimensionality )

  if ( dimensionality == "space") {
    ii = which(
      M$tag==tag &
      M$AUID %in% AUID
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    res$AUID = AUID
    res$matchfrom = list( AUID=M$AUID[ii] )
    res$matchto   = list( AUID=res$AUID )
  }

  if ( dimensionality == "space-year") {
    ii = which(
      M$tag==tag &
      M$AUID %in% AUID &
      M$year %in% year
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    res$AUID = AUID
    res$year = year
    res$matchfrom = list( AUID=M$AUID[ii], year=M$year[ii] )
    res$matchto   = list( AUID=res$AUID,   year=res$year )
  }

  if ( dimensionality == "space-year-season") {
    ii = which(
      M$tag==tag &
      M$AUID %in% AUID &
      M$year %in% year
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    res$AUID = AUID
    res$year = year
    res$dyear = dyear
    res$matchfrom = list( AUID=M$AUID[ii], year=M$year[ii], dyear=M$dyear[ii] )
    res$matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
  }

  res$ii = ii

  return(res)
}
