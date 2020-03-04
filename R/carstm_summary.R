
carstm_summary = function( p=NULL, operation="compute", extrapolation_limit=NA, extrapolation_replacement=NA, ... ) {

  # require areal_units_fn,

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ) )
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  required.vars = c("areal_units_fn", "inputdata_spatial_discretization_planar_km", "inputdata_temporal_discretization_yr", "variabletomodel",
    "carstm_modelengine", "modeldir", "carstm_model_label", "yrs", "variabletomodel" )

  for (i in required.vars) {
    if (!exists(i, p)) {
      message( "Missing parameter" )
      message( i )
      stop()
    }
  }

  # same file naming as in carstm ..
  outputdir = file.path(p$modeldir, p$carstm_model_label)
  areal_units_fns = p$areal_units_fn

  if (exists( "inputdata_spatial_discretization_planar_km", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_spatial_discretization_planar_km, 6),   sep="_" )
  if (exists( "inputdata_temporal_discretization_yr", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  areal_units_fns_suffix = paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine,  "rdata", sep="." )

  fn_res = file.path( outputdir, paste("carstm_modelled_results", areal_units_fns_suffix, sep="." ) )

  if (operation=="load") {  # carstm_model.*carstm_modelled
    if (file.exists(fn_res)) {
      load( fn_res)
      return( res )
    }
  }


  # construct meanweights matrix used to convert number to weight
  sppoly = areal_units( p=p )

  M = snowcrab_carstm( p=p, DS="carstm_inputs" )
  M$yr = M$year  # req for meanweights
  # mean weight by auidxyear


  # init results list
  year = NULL
  if ( p$aegis_dimensionality %in% c("space-year", "space-year-dyear") ) {
    year = as.character( p$yrs )
    M$year = as.character(M$year)
  }

  dyear = NULL
  if ( p$aegis_dimensionality %in% c("space-year-dyear") ){
    dyear = as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) )
    M$dyear = as.character( discretize_data( M$dyear, p$discretization[["dyear"]] ) )
  }


  fit = carstm_model( p=p, DS="carstm_modelled_fit", carstm_model_label=p$carstm_model_label ) # to load currently saved res

  # initialize results in array format
  res = NULL
  res = carstm_inla_destack( inputdata=M, dimensionality=p$aegis_dimensionality, AUID=sppoly[["AUID"]], year=year, dyear=dyear )

  # match conditions for predictions .. ii are locations of predictions in "fit"
  if ( grepl("glm", p$carstm_modelengine) |  grepl("gam", p$carstm_modelengine) ) {

    if ( p$aegis_dimensionality == "space") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = data.frame( matrix( unlist(strsplit( withsolutions, ":", fixed=TRUE)), ncol=1, byrow=TRUE), stringsAsFactors=FALSE )
      ws_matchfrom = list( AUID=as.character(withsolutions[,1])  )
      NA_mask = reformat_to_array( input =rep(1, nrow(withsolutions)), matchfrom=ws_matchfrom, matchto=res$matchto )
    }
    if ( p$aegis_dimensionality == "space-year") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      withsolutions = data.frame( matrix( unlist(strsplit( withsolutions, ":", fixed=TRUE)), ncol=2, byrow=TRUE), stringsAsFactors=FALSE )
      ws_matchfrom = list( AUID=as.character(withsolutions[,1]), year=as.character(withsolutions[,2]) )
      NA_mask = reformat_to_array( input =rep(1, nrow(withsolutions)), matchfrom=ws_matchfrom, matchto=res$matchto )
    }
    if ( p$aegis_dimensionality == "space-year-season") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("dyear", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      withsolutions = data.frame( matrix( unlist(strsplit( withsolutions, ":", fixed=TRUE)), ncol=3, byrow=TRUE), stringsAsFactors=FALSE )
      ws_matchfrom = list( AUID=as.character(withsolutions[,1]), year=as.character(withsolutions[,2]), dyear=as.character(withsolutions[,3])  )
      NA_mask = reformat_to_array( input =rep(1, nrow(withsolutions)), matchfrom=ws_matchfrom, matchto=res$matchto )
    }


    if ( grepl("glm", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[res$ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    }


    if ( grepl("gam", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[res$ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + preds$se.fit
      res[[vn]] = reformat_to_array( input =input, matchfrom=res$matchfrom, matchto=res$matchto ) * NA_mask[]
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
    }

  }


  if ( grepl("inla", p$carstm_modelengine) ) {

    fit = inla.hyperpar(fit, dz=0.2, diff.logdens=20 )  # get improved estimates for the hyperparameters

    if (exists("summary.fitted.values", fit)) {

      vn = paste( p$variabletomodel, "predicted", sep=".")
      input = fit$summary.fitted.values[ res$ii, "mean" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ res$ii, "0.025quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ res$ii, "0.975quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_se", sep=".")
      input = fit$summary.fitted.values[ res$ii, "sd" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )

    }
  }


  ## --------- predictions complete ------



  ## --------- start random effects -------

  if ( grepl("inla", p$carstm_modelengine) ) {

    nAUID = length(res$AUID)

    # match conditions for random effects .. ii are locations of predictions in "fit"
    # random effects results ..
    if (exists("summary.random", fit)) {

      if (exists("iid_error", fit$summary.random)) {
        # IID random effects
        vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[res$ii, "mean" ]
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
      }

      if (exists("auid", fit$summary.random)) {

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single nonspatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          kk = which(resout$type=="nonspatial")
          ns_matchfrom = list( AUID=resout$AUID[kk]  )
          ns_matchto   = list( AUID=res$AUID  )
          input = fit$summary.random$auid[ kk, "mean" ]
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          #  nonspatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          kk = which(resout$type=="nonspatial")
          ns_matchfrom = list( AUID=resout$AUID[kk], year=resout$year[kk] )
          ns_matchto   = list( AUID=res$AUID, year=as.character(p$yrs) )
          input = fit$summary.random$auid[ kk, "mean" ]
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # nonspatial at all time slices
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          kk = which(resout$type=="nonspatial")
          ns_matchfrom = list( AUID=resout$AUID[kk], year=as.character(resout$year[kk]), dyear=as.character( discretize_data( resout$dyear[kk], p$discretization[["dyear"]] ) ) )
          ns_matchto   = list( AUID=res$AUID,   year=as.character(p$yrs),      dyear=as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) ) )
          input = fit$summary.random$auid[ kk, "mean" ]
        }

        vn = paste( p$variabletomodel, "random_auid_nonspatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=ns_matchfrom, matchto=ns_matchto )
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single spatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          kk = which(resout$type=="spatial")
          sp_matchfrom = list( AUID=resout$AUID[kk]  )
          sp_matchto   = list( AUID=res$AUID  )
          input = fit$summary.random$auid[ kk, "mean" ]  # offset structure due to bym2
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          # spatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          kk = which(resout$type=="spatial")
          sp_matchfrom = list( AUID=resout$AUID[kk], year=resout$year[kk] )
          sp_matchto   = list( AUID=res$AUID, year=p$yrs )
          input = fit$summary.random$auid[ kk, "mean" ]  # offset structure due to bym2
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # at every time slice
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          kk = which(resout$type=="spatial")
          sp_matchfrom = list( AUID=resout$AUID[kk], year=as.character(resout$year[kk]), dyear=as.character( discretize_data( resout$dyear[kk], p$discretization[["dyear"]] ) ) )
          sp_matchto   = list( AUID=res$AUID,   year=as.character(p$yrs),      dyear=as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) ) )
          input = fit$summary.random$auid[ kk, "mean" ]  # offset structure due to bym2
        }
        vn = paste( p$variabletomodel, "random_auid_spatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=sp_matchfrom, matchto=sp_matchto )
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

      }
    }
  }

  save( res, file=fn_res, compress=TRUE)

  return (res)

}



