
  carstm_model_gam = function( O, data, fn_fit=tempfile(pattern="fit_", fileext=".RDS"), 
    compress="gzip", compression_level=1,  redo_fit=TRUE , ... ) {
    
    # TODO:: assumes a fixed name convention .. look at carstm_model_inla to disconnect and use O$vn$*

    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (class(data)=="character") assign("data", eval(parse(text=data) ) )

    fn_res=tempfile(pattern="res_", fileext=".RDS")

    vn = O$variabletomodel

    if (exists("data_transformation", O)) data[, vn]  = O$data_transformation$forward( data[, vn] ) # make all positive

    fit = NULL

    if (redo_fit) {
      fit = try( gam( O$formula , data=data, family=O$family ) )

      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")

      message( "Saving carstm fit: ", fn_fit )

      read_write_fast( data=fit, file=fn_fit, compress=compress, compression_level=compression_level )

    }

    if (is.null(fit)) load( fn_fit )

    # do the computations here as fit can be massive ... best not to copy, etc ..
    message( "Computing summaries ..." )


    # O = list()

    O[["summary"]] = summary(fit)

  # row indices for predictions
    S = as.character( O[["AUID"]] )
    nAUID = length(S)

 
    if ( O$dimensionality == "space") {
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      MM = paste( data$AUID, data$year, sep=":")
      ipred = match(withsolutions, MM)
      mfrom = list( S=data$AUID[ipred] )
      mto   = list( S=S )
    }

    if ( O$dimensionality == "space-time") {
      T = as.character( O$yrs )
      data$year = as.character(data$year)
      # filter by AUID and years in case additional data in other areas and times are used in the input data
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( data$AUID, data$year, sep=":")
      ipred = ipred[ match(withsolutions, MM) ]
      mfrom = list( S=data$AUID[ipred], T=data$year[ipred]  )
      mto   = list( S=S, T=T )
    }

    if ( O$dimensionality == "space-time-cyclic") {
      T = as.character( O$yrs )
      U = as.character( discretize_data( (O$dyears + diff(O$dyears)[1]/2), O$discretization[["dyear"]] ) )

      data$year = as.character(data$year)
      data$dyear = as.character( discretize_data( data$dyear, O$discretization[["dyear"]] ) )

      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("dyear", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( data$AUID, data$year, data$dyear, sep=":")
      ipred = match(withsolutions, MM)
      mfrom = list( S=data$AUID[ipred], T=data$year[ipred], U=data$dyear[ipred] )
      mto   = list( S=S, T=T, U=U )
    }


    preds = predict( fit, newdata=data[ipred,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
    vn =  paste(O$variabletomodel, "predicted", sep=".")
    input = preds$fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", O$family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", O$family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", O) ) O[[vn]] = O$data_transformation$backward( O[[vn]] ) # make all positive

    vn =  paste(O$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )

    vn =  paste(O$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - 1.96*preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", O$family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", O$family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", O) ) O[[vn]] = O$data_transformation$backward( O[[vn]] ) # make all positive

    vn =  paste(O$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + 1.96*preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", O$family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", O$family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", O) ) O[[vn]] = O$data_transformation$backward( O[[vn]] ) # make all positive

    read_write_fast( data=O, file=fn_res, compress=compress, compression_level=compression_level )

    message( "carstm summary saved as: ", fn_res )

    return(O)

  }