
carstm_model = function( p, M=NULL, DS="redo", improve.hyperparam.estimates=FALSE, file_compress_method=FALSE, ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  sppoly = areal_units( p=p )  # required by car fit
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

  fn_fit = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
  fn_res = carstm_filenames( p=p, returntype="carstm_modelled_summary", areal_units_fn=areal_units_fn )
  outputdir = dirname(fn_fit)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
  
  fit = NULL
  if (DS=="carstm_modelled_fit") {
    if (file.exists(fn_fit)) load( fn_fit )
    if (is.null(fit)) message("carstm modelled fit not found.")
    return( fit )
  }

  res = NULL
  if (DS=="carstm_modelled_summary") {  # carstm_model.*carstm_modelled
    if (file.exists(fn_res)) load( fn_res)
    if (is.null(res)) message("carstm summary not found.")
    return( res )
  }
 
  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )

  # INLA does not like duplicates ... causes optimizer to crash frequently
  eps = exp( log( .Machine$double.eps ) / 2)  # ~ 1.5e-8

  vn = p$variabletomodel
  M[, vn]  = M[, vn]  + runif( nrow(M), -eps, eps )

  if (exists("data_transformation", p)) M[, vn]  = p$data_transformation$forward( M[, vn] ) # make all positive


  mrange = NULL
  # get hyper param scalings
  if ( grepl("inla", p$carstm_modelengine) ) {

    # hyperparms
    j = which( is.finite(M[,vn]) )
    mrange = range( M[ j, vn ]  )  # on data scale not internal

    if (!exists("p$carstm_model_family", p )) p$carstm_model_family = "normal"

    if ( p$carstm_model_family == "lognormal" ) {
      m = log( M[ j, vn ])
    } else if ( grepl( ".*poisson", p$carstm_model_family)) {
      m = log( M[ j, vn ] / M[ j, "data_offset" ]  )
      mrange = range( M[ j, vn ]/ M[ j, "data_offset" ]  )  # on data scale not internal
      mrange = mrange * median(M[ M$tag=="predictions", "data_offset" ] )
    } else if ( p$carstm_model_family =="binomial" )  {
      m = M[ j, vn ]
    } else {
      m = M[ j, vn ]
    }

    H = carstm_hyperparameters( sd(m), alpha=0.5, median(m) )
    m = NULL


    p = parameters_add_without_overwriting( p,
      options.control.family = inla.set.control.family.default()
    )



    gc()

    fit  = NULL

      if (!exists("options.control.inla", p )) p$options.control.inla = list(
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive"), # default h=0.02
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.05, cmin=0, tolerance=1e-9),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.1, cmin=0),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.001, cmin=0), # default h=0.02 ?or 0.01
        list( optimise.strategy="smart", stupid.search=TRUE, strategy="adaptive", h=0.0001, cmin=0), # default h=0.02
        list( optimise.strategy="smart", h=0.2 ),
        list( optimise.strategy="smart", h=0.4 ),
        list( stupid.search=TRUE, fast=FALSE, step.factor=0.1),
        list( stupid.search=TRUE, cmin=0)
      )


      for ( civ in 1:length(p$options.control.inla)) {
        fit = try( inla( p$carstm_model_formula , data=M, family=p$carstm_model_family,
          control.compute=list(dic=TRUE, waic=FALSE, cpo=FALSE, config=FALSE),
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list(compute=FALSE, link=1 ),
          control.fixed= list(mean.intercept=0, prec.intercept=0.001, mean=0, prec=0.001),
          control.family = p$options.control.family,
          control.inla   = p$options.control.inla[[civ]],
          verbose=TRUE
        ))
        if (!inherits(fit, "try-error" )) break()
      }

      if (inherits(fit, "try-error" )) {
        message("If you are using MSWindows and you get a popup complaining about inla stopped working,")
        message("try setting the flag in the following link to 1, using regedit. Be careful.")
        message("e.g., see: https://monitormyweb.com/guides/how-to-disable-stopped-working-message-in-windows")
        stop( "solution did not converge")
      }

      # to improve hyper param estimates..
      if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.25, diff.logdens=18 )  # get improved estimates for the hyperparameters

  }

  if ( grepl("glm", p$carstm_modelengine) ) {
    fit = try( glm( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )
  }

  if ( grepl("gam", p$carstm_modelengine) ) {
    fit = try( gam( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )
  }

  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")

  message( "Saving carstm fit (this can be slow): ", fn_fit )

  save( fit, file=fn_fit, compress=file_compress_method )


  # ----------------
  # summarize
  # do the computations here as fit can be massive ... best not to copy, etc ..
  message( "Computing summaries ..." )


  # results go here
  res = list( M=M, dimensionality=p$aegis_dimensionality, summary=summary(fit), 
    sppoly=sppoly, fn_res=fn_res 
  )
  
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
      if ( grepl( ".*poisson", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( ".*poisson", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( ".*poisson", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
    }


    if ( grepl("gam", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[res$i_preds,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( ".*poisson", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( ".*poisson", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + 1.96*preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( ".*poisson", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
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
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive


      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ res$i_preds, "0.025quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ res$i_preds, "0.975quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      if ( grepl( ".*lognormal", p$carstm_model_family)) res[[vn]] = exp(res[[vn]])
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
        # carstm_map( res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )

        input = fit$summary.random$auid[ res$i_spatial, "mean" ]  # offset structure due to bym2
        vn = paste( p$variabletomodel, "random_auid_spatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$sp_matchfrom, matchto=res$sp_matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
        # carstm_map( res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )

      }
    }
  }

  save( res, file=fn_res, compress=file_compress_method )

  message( "carstm summary saved as: ", fn_res )

  print( res$summary )

  # if ( grepl("inla", p$carstm_modelengine) ) {

  #   # a few plots
  #   # copied from INLA::plot.inla and streamlined to reduce access to the fit object which can be huge
  #   print( plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=TRUE ) )

  # } else {
  #   print( plot(fit) )
  # }

  return( fit )
}
