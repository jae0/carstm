
carstm_summary = function( p=NULL, fit=NULL, M=NULL, sppoly=NULL, operation="load", mrange=NULL, improve.hyperparam.estimates=FALSE, ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  if (is.null(sppoly))  sppoly = areal_units( p=p )  # will redo if not found
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

  fn_res = carstm_filenames( p=p, returntype="carstm_modelled_results", areal_units_fn=areal_units_fn )
  outputdir = dirname(fn_res)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  res = NULL
  if (operation=="load") {  # carstm_model.*carstm_modelled
    if (file.exists(fn_res)) load( fn_res)
    if (is.null(res)) message("carstm_summary: Modelled surface not found, this needs to be created or alternate sppoly specified.")
    return( res )
  }

  # to improve hyper param estimates..
  if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.25, diff.logdens=18 )  # get improved estimates for the hyperparameters


  # results go here
  res = list( M=M, dimensionality = p$aegis_dimensionality )

  # row indices for predictions
  if ( p$aegis_dimensionality == "space") {
    AUID = sppoly[["AUID"]]
    res$i_preds = which(
      M$tag=="predictions" &
      M$AUID %in% AUID
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    res$AUID = AUID
    res$matchfrom = list( AUID=M$AUID[res$i_preds] )
    res$matchto   = list( AUID=res$AUID )
  }

  if ( p$aegis_dimensionality == "space-year") {
    AUID = sppoly[["AUID"]]
    year = as.character( p$yrs )
    M$year = as.character(M$year)

    res$i_preds = which(
      M$tag=="predictions" &
      M$AUID %in% AUID &
      M$year %in% year
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data

    res$AUID = AUID
    res$year = year
    res$matchfrom = list( AUID=M$AUID[res$i_preds], year=M$year[res$i_preds] )
    res$matchto   = list( AUID=res$AUID,   year=res$year )
  }

  if ( p$aegis_dimensionality == "space-year-season") {
    AUID = sppoly[["AUID"]]
    year = as.character( p$yrs )
    dyear = as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) )
    M$year = as.character(M$year)
    M$dyear = as.character( discretize_data( M$dyear, p$discretization[["dyear"]] ) )

    res$i_preds = which(
      M$tag=="predictions" &
      M$AUID %in% AUID &
      M$year %in% year
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    res$AUID = AUID
    res$year = year
    res$dyear = dyear
    res$matchfrom = list( AUID=M$AUID[res$i_preds], year=M$year[res$i_preds], dyear=M$dyear[res$i_preds] )
    res$matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
  }


  nAUID = length(res$AUID)


  if ( grepl("glm", p$carstm_modelengine) |  grepl("gam", p$carstm_modelengine) ) {

    if ( p$aegis_dimensionality == "space") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      MM = paste( M$AUID, M$year, sep=":")[res$i_preds]
      res$i_preds = match(withsolutions, MM)
      res$matchfrom = list( AUID=M$AUID[res$i_preds] )
    }

    if ( p$aegis_dimensionality == "space-year") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( M$AUID, M$year, sep=":")[res$i_preds]
      res$i_preds = res$i_preds[ match(withsolutions, MM) ]
      res$matchfrom = list( AUID=M$AUID[res$i_preds], year=M$year[res$i_preds]  )
    }

    if ( p$aegis_dimensionality == "space-year-dyear") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("dyear", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( M$AUID, M$year, M$dyear, sep=":")[res$i_preds]
      res$i_preds = match(withsolutions, MM)
      res$matchfrom = list( AUID=M$AUID[res$i_preds], year=M$year[res$i_preds], dyear=M$dyear[res$i_preds] )
    }

    if ( grepl("glm", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[res$i_preds,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*poisson", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*poisson", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*poisson", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
    }


    if ( grepl("gam", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[res$i_preds,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*poisson", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*poisson", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*poisson", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
    }
  }


  if ( grepl("inla", p$carstm_modelengine) ) {
    # row indices for spatial effects (if any)
    if (exists("summary.random", fit)) {

      if (exists("auid", fit$summary.random)) {

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single nonspatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          res$i_nonspatial = which(resout$type=="nonspatial")
          res$ns_matchfrom = list( AUID=resout$AUID[res$i_nonspatial]  )
          res$ns_matchto   = list( AUID=res$AUID  )

        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          #  nonspatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          res$i_nonspatial = which(resout$type=="nonspatial")
          res$ns_matchfrom = list( AUID=resout$AUID[res$i_nonspatial], year=resout$year[res$i_nonspatial] )
          res$ns_matchto   = list( AUID=res$AUID,   year=res$year  )
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # nonspatial at all time slices
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          res$i_nonspatial = which(resout$type=="nonspatial")
          res$ns_matchfrom = list( AUID=resout$AUID[res$i_nonspatial], year=as.character(resout$year[res$i_nonspatial]), dyear=as.character( discretize_data( resout$dyear[res$i_nonspatial], p$discretization[["dyear"]] ) ) )
          res$ns_matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
        }

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single spatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          res$i_spatial = which(resout$type=="spatial")
          res$sp_matchfrom = list( AUID=resout$AUID[res$i_spatial]  )
          res$sp_matchto   = list( AUID=res$AUID  )

        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          # spatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          res$i_spatial = which(resout$type=="spatial")
          res$sp_matchfrom = list( AUID=resout$AUID[res$i_spatial], year=resout$year[res$i_spatial] )
          res$sp_matchto   = list( AUID=res$AUID,   year=res$year  )
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # at every time slice
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          res$i_spatial = which(resout$type=="spatial")
          res$sp_matchfrom = list( AUID=resout$AUID[res$i_spatial], year=as.character(resout$year[res$i_spatial]), dyear=as.character( discretize_data( resout$dyear[res$i_spatial], p$discretization[["dyear"]] ) ) )
          res$sp_matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
        }

      }
    }


    # predictions

    if (exists("summary.fitted.values", fit)) {

      vn = paste( p$variabletomodel, "predicted", sep=".")
      input = fit$summary.fitted.values[ res$i_preds, "mean" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )

      NA_mask = NULL
      if ( exists("carstm_predict_force_range", p)) {
        if (p$carstm_predict_force_range) {
          if (!is.null(mrange)) {
            # only an issue for factorial models ... missing locations get predicted to absurd values as there is no information .. filter out
            NA_mask = which(res[[vn]]  > max(mrange, na.rm=TRUE) )
          }
        }
      }
      if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive


      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ res$i_preds, "0.025quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ res$i_preds, "0.975quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      if ( grepl( "family.*=.*lognormal", p$carstm_model_call)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_se", sep=".")
      input = fit$summary.fitted.values[ res$i_preds, "sd" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA

    }

    ## --------- predictions complete ------



    ## --------- start random effects -------

    # match conditions for random effects .. i_preds are locations of predictions in "fit"
    # random effects results ..
    if (exists("summary.random", fit)) {

      if (exists("iid_error", fit$summary.random)) {
        # IID random effects
        vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[res$i_preds, "mean" ]
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      }

      if (exists("auid", fit$summary.random)) {

        input = fit$summary.random$auid[ res$i_nonspatial, "mean" ]
        vn = paste( p$variabletomodel, "random_auid_nonspatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$ns_matchfrom, matchto=res$ns_matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

        input = fit$summary.random$auid[ res$i_spatial, "mean" ]  # offset structure due to bym2
        vn = paste( p$variabletomodel, "random_auid_spatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$sp_matchfrom, matchto=res$sp_matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

      }
    }
  }

  save( res, file=fn_res, compress=TRUE)

  if (operation!="load") res=fn_res  # if computing as file was not found in load step, then return data, else, file name

  return (res)

}



