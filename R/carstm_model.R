
carstm_model = function( p, M=NULL, DS="redo" ) {

  auids = p$auid
  if (exists( "inputdata_spatial_discretization_planar_km", p )) auids = paste( auids, p$inputdata_spatial_discretization_planar_km ,   sep="_" )
  if (exists( "inputdata_temporal_discretization_yr", p )) auids = paste( auids, round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )

  auids_suffix = paste( auids, p$variabletomodel, p$carstm_modelengine,  "rdata", sep="." )
  fn = file.path( p$modeldir, paste("carstm_modelled_results", auids_suffix, sep="." ) )
  fn_fit = file.path( p$modeldir, paste( "carstm_modelled_fit", auids_suffix, sep=".") )

  if (DS %in% c("carstm_modelled_fit", "carstm_modelled"))  {
    if (DS=="carstm_modelled") {
      if (file.exists(fn)) {
        load( fn)
        return( res )
      }
    }
    if (DS=="carstm_modelled_fit") {
      if (file.exists(fn_fit)) {
        load( fn_fit )
        return( fit )
      }
    }
  }

  # prediction surface
  sppoly = areal_units( p=p )  # will redo if not found

  # init results list
  res = list(
    StrataID = sppoly[["StrataID"]],
    strata   = as.numeric(sppoly[["StrataID"]] )  )

  if (exists("data_transformation", p)) M[, p$variabletomodel]  = p$data_transformation$forward( M[, p$variabletomodel] ) # make all positive

  if ( grepl("inla", p$carstm_modelengine) ) {
    # hyperparms
    H = carstm_hyperparameters( sd(M[,p$variabletomodel], na.rm=TRUE), alpha=0.5, median( M[,p$variabletomodel], na.rm=TRUE) )
    M$strata  = as.numeric( M$StrataID)
    M$zi = discretize_data( M$z, p$discretization$z )

    M$iid_error = 1:nrow(M) # for inla indexing for set level variation
    if ( p$aegis_dimensionality == "space-year") {
      M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
      M$year = floor(M$tiyr)
    }
    if ( p$aegis_dimensionality == "space-year-season") {
      M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
      M$year = floor(M$tiyr)
      M$dyear  =  factor( as.character( trunc(  (M$tiyr - M$year )/ p$tres )*p$tres), levels=p$dyears)
    }

  }

  fit  = NULL
  assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")
  save( fit, file=fn_fit, compress=TRUE )


  # match conditions for predictions .. ii are locations of predictions in "fit"
  if ( p$aegis_dimensionality == "space") {
    ii = which(
      M$tag=="predictions" &
      M$strata %in% res$strata
    )  # filter by strata and years in case additional data in other areas and times are used in the input data
    matchfrom = list( strata=M$strata[ii] )
    matchto   = list( strata=res$strata )
  }

  if ( p$aegis_dimensionality == "space-year") {
    ii = which(
      M$tag=="predictions" &
      M$strata %in% res$strata &
      M$year %in% p$yrs
    )  # filter by strata and years in case additional data in other areas and times are used in the input data
    matchfrom = list( strata=M$strata[ii], year=as.character(M$year[ii]) )
    matchto   = list( strata=res$strata, year=as.character(p$yrs)  )
  }

  if ( p$aegis_dimensionality == "space-year-season") {
    ii = which(
      M$tag=="predictions" &
      M$strata %in% res$strata &
      M$year %in% p$yrs
    )  # filter by strata and years in case additional data in other areas and times are used in the input data
    matchfrom = list( strata=M$strata[ii], year=as.character(M$year[ii]), dyear=M$dyear[ii] )
    matchto   = list( strata=res$strata, year=as.character(p$yrs), dyear=factor(p$dyears) )
 }


  if ( grepl("glm", p$carstm_modelengine) ) {

    preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

    vn =  paste(p$variabletomodel, "predicted", sep=".")
    input = preds$fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )

    vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

  }


  if ( grepl("gam", p$carstm_modelengine) ) {

    preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

    vn =  paste(p$variabletomodel, "predicted", sep=".")
    input = preds$fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )

    vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
  }


  if ( grepl("inla", p$carstm_modelengine) ) {

    if (exists("summary.fitted.values", fit)) {
      vn = paste( p$variabletomodel, "predicted", sep=".")
      input = fit$summary.fitted.values[ ii, "mean" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_se", sep=".")
      input = fit$summary.fitted.values[ ii, "sd" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ ii, "0.025quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ ii, "0.975quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
    }
  }


  ## --------- predictions complete ------



  ## --------- start random effects -------

  if ( grepl("inla", p$carstm_modelengine) ) {

    nstrata = length(res$StrataID)

    # match conditions for random effects .. ii are locations of predictions in "fit"
    # random effects results ..
    if (exists("summary.random", fit)) {

      if (exists("iid_error", fit$summary.random)) {
        # IID random effects
        if ( p$aegis_dimensionality == "space") {
          ii = which(
            M$tag=="predictions" &
            M$strata %in% res$strata
          )  # filter by strata and years in case additional data in other areas and times are used in the input data
          matchfrom = list( strata=M$strata[ii] )
          matchto   = list( strata=res$strata )
        }

        if ( p$aegis_dimensionality == "space-year") {
          ii = which(
            M$tag=="predictions" &
            M$strata %in% res$strata &
            M$year %in% p$yrs
          )  # filter by strata and years in case additional data in other areas and times are used in the input data
          matchfrom = list( strata=M$strata[ii], year=as.character(M$year[ii]) )
          matchto   = list( strata=res$strata, year=as.character(p$yrs)  )
        }

        if ( p$aegis_dimensionality == "space-year-season") {
          ii = which(
            M$tag=="predictions" &
            M$strata %in% res$strata &
            M$year %in% p$yrs
          )  # filter by strata and years in case additional data in other areas and times are used in the input data
          matchfrom = list( strata=M$strata[ii], year=as.character(M$year[ii]), dyear=M$dyear[ii] )
          matchto   = list( strata=res$strata, year=as.character(p$yrs), dyear=factor(p$dyears) )
        }
        vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[ii, "mean" ]
        res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      }

      if (exists("strata", fit$summary.random)) {

        if (nrow(fit$summary.random$strata) == nstrata*2) {
          # a single nonspatial effect (no grouping across time)
          kk = 1:nstrata
          matchfrom = list( strata=fit$summary.random$strata$ID[kk]  )
          matchto   = list( strata=res$strata  )
        } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$ny ) {
          #  nonspatial effects grouped by year
          kk = ii
          matchfrom = list( strata=M$strata[kk], year=M$year[kk] )
          matchto   = list( strata=res$strata, year=p$yrs )
        } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$nt ) {
          # nonspatial at all time slices
          kk = ii
          matchfrom = list( strata=M$strata[kk], year=M$year[kk], dyear=M$dyear[kk] )
          matchto   = list( strata=res$strata, year=p$yrs, dyear=factor(p$dyears) )
        }
        vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
        input = fit$summary.random$strata[ kk, "mean" ]
        res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )


        if (nrow(fit$summary.random$strata) == nstrata*2) {
          # a single spatial effect (no grouping across time)
          kk = 1:nstrata
          matchfrom = list( strata=fit$summary.random$strata$ID[kk]  )
          matchto   = list( strata=res$strata  )
        } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$ny ) {
          # spatial effects grouped by year
          kk = ii
          matchfrom = list( strata=M$strata[kk], year=M$year[kk] )
          matchto   = list( strata=res$strata, year=p$yrs )
        } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$nt ) {
          # at every time slice
          kk = ii
          matchfrom = list( strata=M$strata[kk], year=M$year[kk], dyear=M$dyear[kk] )
          matchto   = list( strata=res$strata, year=p$yrs, dyear=factor(p$dyears) )
        }

        vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
        input = fit$summary.random$strata[ kk+max(kk), "mean" ]  # offset structure due to bym2
        res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

      }
    }
  }

  save( res, file=fn, compress=TRUE )
  return(res)
}
