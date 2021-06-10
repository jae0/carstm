
  carstm_model_glm = function( p, M, fn_fit, file_compress_method=FALSE ) {
    

    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (class(M)=="character") assign("M", eval(parse(text=M) ) )


    vn = p$variabletomodel

    if (exists("data_transformation", p)) M[, vn]  = p$data_transformation$forward( M[, vn] ) # make all positive

    fit = try( glm( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )

    if (is.null(fit)) warning("model fit error")
    if ("try-error" %in% class(fit) ) warning("model fit error")

    message( "Saving carstm fit: ", fn_fit )

    save( fit, file=fn_fit, compress=file_compress_method )

    # do the computations here as fit can be massive ... best not to copy, etc ..
    message( "Computing summaries ..." )
 
    O = list()

    O[["summary"]] = summary(fit)

  # row indices for predictions
    S = as.character( p[["AUID"]] )

    nAUID = length(S)
 
    if ( p$aegis_dimensionality == "space") {
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      MM = paste( M$AUID, M$year, sep=":")
      ipred = match(withsolutions, MM)
      mfrom = list( S=M$AUID[ipred] )
      mto   = list( S=S )
    }

    if ( p$aegis_dimensionality == "space-year") {
      T = as.character( p$yrs )
      M$year = as.character(M$year)
      # filter by AUID and years in case additional data in other areas and times are used in the input data
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( M$AUID, M$year, sep=":")
      ipred = ipred[ match(withsolutions, MM) ]
      mfrom = list( S=M$AUID[ipred], T=M$year[ipred]  )
      mto   = list( S=S, T=T )
    }

    if ( p$aegis_dimensionality == "space-year-season") {
      T = as.character( p$yrs )
      U = as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) )

      M$year = as.character(M$year)
      M$dyear = as.character( discretize_data( M$dyear, p$discretization[["dyear"]] ) )

      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("dyear", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( M$AUID, M$year, M$dyear, sep=":")
      ipred = match(withsolutions, MM)
      mfrom = list( S=M$AUID[ipred], T=M$year[ipred], U=M$dyear[ipred] )
      mto   = list( S=S, T=T, U=U )
    }

    preds = predict( fit, newdata=M[ipred,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
    
    O = list()

    O[["summary"]] = summary(fit)

    vn =  paste(vn, "predicted", sep=".")
    input = preds$fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", p) ) O[[vn]] = data_transformation$backward( O[[vn]] ) # make all positive

    vn =  paste(vn, "predicted_se", sep=".")
    input = preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )

    vn =  paste(vn, "predicted_lb", sep=".")
    input = preds$fit - 1.96*preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", p) ) O[[vn]] = data_transformation$backward( O[[vn]] ) # make all positive

    vn =  paste(vn, "predicted_ub", sep=".")
    input = preds$fit + 1.96*preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", p) ) O[[vn]] = data_transformation$backward( O[[vn]] ) # make all positive

    return(O)

  }