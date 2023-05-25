
  carstm_model_glm = function( O, data, fn_fit=tempfile(pattern="fit_", fileext=".rdata"), fn_res=tempfile(pattern="res_", fileext=".rdata"), compress=TRUE,   redo_fit=TRUE , ...  ) {
    
    # TODO:: assumes a fixed name convention .. look at carstm_model_inla to disconnect and use O$vn$*

    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (class(data)=="character") assign("data", eval(parse(text=data) ) )


    vn = O$variabletomodel

    if (exists("data_transformation", O)) data[, vn]  = O$data_transformation$forward( data[, vn] ) # make all positive

    fit = NULL

    if (redo_fit) {
      fit = try( glm( O$formula , data=data, family=O$family ) )

      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")

      message( "Saving carstm fit: ", fn_fit )

      save( fit, file=fn_fit, compress=compress )

    }

    if (is.null(fit)) load( fn_fit )

    # do the computations here as fit can be massive ... best not to copy, etc ..
    message( "Computing summaries ..." )
 
    # O = list()


    O[["summary"]] = summary(fit)

  # row indices for predictions
    S = as.character( O[["AUID"]] )

    nAUID = length(S)
 
    # parsing formula is better as it makes it independent of storage dim (aegis_dimensionality); see carstm_model_inla
    if (O$aegis_dimensionality=="space") O$dimensionality = "space"  
    if (O$aegis_dimensionality=="space-year") O$dimensionality = "space-time"  
    if (O$aegis_dimensionality=="space-year-season") O$dimensionality = "space-time-cyclic"  
  
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
    
    O = list()

    O[["summary"]] = summary(fit)

    vn =  paste(vn, "predicted", sep=".")
    input = preds$fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", O$family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", O$family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", O) ) O[[vn]] = data_transformation$backward( O[[vn]] ) # make all positive

    vn =  paste(vn, "predicted_se", sep=".")
    input = preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )

    vn =  paste(vn, "predicted_lb", sep=".")
    input = preds$fit - 1.96*preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", O$family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", O$family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", O) ) O[[vn]] = data_transformation$backward( O[[vn]] ) # make all positive

    vn =  paste(vn, "predicted_ub", sep=".")
    input = preds$fit + 1.96*preds$se.fit
    O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
    if ( grepl( ".*poisson", O$family)) O[[vn]] = exp(O[[vn]])
    if ( grepl( ".*lognormal", O$family)) O[[vn]] = exp(O[[vn]])
    if (exists("data_transformation", O) ) O[[vn]] = data_transformation$backward( O[[vn]] ) # make all positive

    save( O, file=fn_res, compress=compress )

    message( "carstm summary saved as: ", fn_res )

    return(O)

  }