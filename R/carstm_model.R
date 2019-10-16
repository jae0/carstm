
carstm_model = function( p, M=NULL, DS="redo" ) {

  auids = paste(  p$auid, p$inputdata_spatial_discretization_planar_km,
    round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )
  auids_suffix = paste( p$variabletomodel, p$carstm_modelengine, auids, "rdata", sep="." )
  fn = file.path( p$modeldir, paste("carstm_modelled", auids_suffix, sep="." ) )
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
  res = list(StrataID = sppoly[["StrataID"]])  # init results list
  res$strata = as.numeric(res$StrataID)


  if ( grepl("glm", p$carstm_modelengine) ) {
    #nothing to do, in case input dat needs tweaking
  }

  if ( grepl("gam", p$carstm_modelengine) ) {
    #nothing to do, in case input dat needs tweaking
  }

  if ( grepl("inla", p$carstm_modelengine) ) {
    # hyperparms
    H = carstm_hyperparameters( sd(M[,p$variabletomodel], na.rm=TRUE), alpha=0.5, median( M[,p$variabletomodel], na.rm=TRUE) )
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


  # predictions .. ii are locations of predictions in "fit"
  if ( p$aegis_dimensionality == "space") {
    ii = which(
      M$tag=="predictions" &
      M$strata %in% res$strata
    )  # filter by strata and years in case additional data in other areas and times are used in the input data
    matchfrom_preds = list( strata=M$strata[ii] )
    matchto_preds   = list( strata=res$strata )
  }

  if ( p$aegis_dimensionality == "space-year") {
    ii = which(
      M$tag=="predictions" &
      M$strata %in% res$strata &
      M$year %in% p$yrs
    )  # filter by strata and years in case additional data in other areas and times are used in the input data
    matchfrom_preds = list( strata=M$strata[ii], year=as.character(M$year[ii]) )
    matchto_preds   = list( strata=res$strata, year=as.character(p$yrs)  )
  }

  if ( p$aegis_dimensionality == "space-year-season") {
    ii = which(
      M$tag=="predictions" &
      M$strata %in% res$strata &
      M$year %in% p$yrs
    )  # filter by strata and years in case additional data in other areas and times are used in the input data
    matchfrom_preds = list( strata=M$strata[ii], year=as.character(M$year[ii]), dyear=M$dyear[ii] )
    matchto_preds   = list( strata=res$strata, year=as.character(p$yrs), dyear=factor(p$dyears) )
 }


  if ( grepl("glm", p$carstm_modelengine) ) {

    preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

    vn =  paste(p$variabletomodel, "predicted", sep=".")
    input = preds$fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )
    if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )

    vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - preds$se.fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )
    if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + preds$se.fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )
    if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

  }


  if ( grepl("gam", p$carstm_modelengine) ) {

    preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

    vn =  paste(p$variabletomodel, "predicted", sep=".")
    input = preds$fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )
    if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )

    vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - preds$se.fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )
    if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + preds$se.fit
    res[vn] = reformat_to_array( input =input, matchfrom=matchfrom_pred, matchto=matchto_pred )
    if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive
  }


  if ( grepl("inla", p$carstm_modelengine) ) {
    nstrata = length(res$StrataID)

    if (exists("summary.fitted.values", fit)) {
      vn = paste( p$variabletomodel, "predicted", sep=".")
      input = fit$summary.fitted.values[ ii, "0.025quant" ]
      res[vn] = reformat_to_array( input=input, matchfrom=matchfrom_pred, matchto=matchto_pred )
      if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_se", sep=".")
      input = fit$summary.fitted.values[ ii, "se" ]
      res[vn] = reformat_to_array( input=input, matchfrom=matchfrom_pred, matchto=matchto_pred )

      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ ii, "mean" ]
      res[vn] = reformat_to_array( input=input, matchfrom=matchfrom_pred, matchto=matchto_pred )
      if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ ii, "0.975quant" ]
      res[vn] = reformat_to_array( input=input, matchfrom=matchfrom_pred, matchto=matchto_pred )
      if (exists("data_transformation", p) ) res[vn] = p$data_transformation$backward( res[vn] ) # make all positive
    }


    # random effects results ..
    if (exists("summary.random", fit)) {

      if (exists("iid_error", fit$summary.random)) {
        # IID random effects
        matchfrom = list( strata=M$strata[ii], year=M$year[ii], dyear=M$dyear[ii] )
        matchto   = list( strata=res$strata, year=p$yrs, dyear=factor(p$dyears) )

        vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[ ii, "mean" ]
        res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="1950", dyear="0") )
      }

      if (exists("strata", fit$summary.random)) {

        if (nrow(fit$summary.random$strata) == nstrata*2) {
          # CAR random effects can be of variable length depending upon model construct:

          # a single spatial and nonspatial effect (no grouping across time)
          jj = 1:nstrata
          matchfrom = list( strata=fit$summary.random$strata$ID[jj]  )
          matchto   = list( strata=res$strata  )

          vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
          input = fit$summary.random$strata[ jj, "mean" ]
          res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

          vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
          input = fit$summary.random$strata[ jj+max(jj), "mean" ]
          res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

          # carstm_plot( p=p, res=res, vn=vn )

        } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$ny ) {

          # spatial and nonspatial effects grouped by year
          matchfrom = list( strata=M$strata[ii], year=M$year[ii] )
          matchto   = list( strata=res$strata, year=p$yrs )

          vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
          input = fit$summary.random$strata[ ii, "mean" ]
          res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

          vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
          input = fit$summary.random$strata[ ii+max(ii), "mean" ]
          res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

          # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000" ) )

        } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$nt ) {

          # need to test/fix ...
          matchfrom = list( StrataID=M$StrataID[ii], year=M$year[ii], dyear=M$dyear[ii] )
          matchto   = list( StrataID=res$StrataID, year=p$yrs, dyear=factor(p$dyears) )

          vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
          input = fit$summary.random$strata[ ii, "mean" ]
          res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

          vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
          input = fit$summary.random$strata[ ii+max(ii), "mean" ]
          res[vn] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

          # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )
        }
      }
    }
  }

  save( res, file=fn, compress=TRUE )
  return(res)
}
