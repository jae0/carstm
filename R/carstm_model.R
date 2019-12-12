
carstm_model = function( p, M=NULL, DS="redo", ... ) {

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ) )
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  areal_units_fns = p$areal_units_fn
  if (exists( "inputdata_spatial_discretization_planar_km", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_spatial_discretization_planar_km, 6),   sep="_" )
  if (exists( "inputdata_temporal_discretization_yr", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )

  areal_units_fns_suffix = paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine,  "rdata", sep="." )
  outputdir = file.path(p$modeldir, p$carstm_model_label)

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  fn = file.path( outputdir, paste("carstm_modelled_results", areal_units_fns_suffix, sep="." ) )
  fn_fit = file.path( outputdir, paste( "carstm_modelled_fit", areal_units_fns_suffix, sep=".") )

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
  res = list( AUID = sppoly[["AUID"]]  )

  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )

  if (exists("data_transformation", p)) M[, p$variabletomodel]  = p$data_transformation$forward( M[, p$variabletomodel] ) # make all positive

  if ( grepl("inla", p$carstm_modelengine) ) {
    # hyperparms
    j = which( is.finite(M[,p$variabletomodel]) )
    if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) {
      m = log( M[ j, p$variabletomodel ])
    } else if ( grepl( "family.*=.*poisson", p$carstm_modelcall)) {
      m = log( M[ j, p$variabletomodel ] / M[ j, "data_offset" ]  )
    } else {
      m = M[,p$variabletomodel]
    }

    H = carstm_hyperparameters( sd(m), alpha=0.5, median(m) )
    m = NULL

    # adjust based upon RAM requirements and ncores
    inla.setOption(num.threads= p$inla_num.threads)
    inla.setOption(blas.num.threads=p$inla_blas.num.threads)

  }

  gc()

  fit  = NULL
  assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")
  save( fit, file=fn_fit, compress=TRUE )


  # match conditions for predictions .. ii are locations of predictions in "fit"
  if ( p$aegis_dimensionality == "space") {
    ii = which(
      M$tag=="predictions" &
      M$AUID %in% res$AUID
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    matchfrom = list( AUID=M$AUID[ii] )
    matchto   = list( AUID=res$AUID )
    if ( grepl("glm", p$carstm_modelengine) |  grepl("gam", p$carstm_modelengine) ) {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ - grep("Intercept", withsolutions) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = data.frame( matrix( unlist(strsplit( withsolutions, ":", fixed=TRUE)), ncol=1, byrow=TRUE), stringsAsFactors=FALSE )
      ws_matchfrom = list( AUID=as.character(withsolutions[,1])  )
      NA_mask = reformat_to_array( input =rep(1, nrow(withsolutions)), matchfrom=ws_matchfrom, matchto=matchto )
    }
  }

  if ( p$aegis_dimensionality == "space-year") {
    ii = which(
      M$tag=="predictions" &
      M$AUID %in% res$AUID &
      M$year %in% p$yrs
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    matchfrom = list( AUID=M$AUID[ii], year=as.character(M$year[ii]) )
    matchto   = list( AUID=res$AUID,   year=as.character(p$yrs)  )

    if ( grepl("glm", p$carstm_modelengine) |  grepl("gam", p$carstm_modelengine) ) {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ - grep("Intercept", withsolutions) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      withsolutions = data.frame( matrix( unlist(strsplit( withsolutions, ":", fixed=TRUE)), ncol=2, byrow=TRUE), stringsAsFactors=FALSE )
      ws_matchfrom = list( AUID=as.character(withsolutions[,1]), year=as.character(withsolutions[,2]) )
      NA_mask = reformat_to_array( input =rep(1, nrow(withsolutions)), matchfrom=ws_matchfrom, matchto=matchto )
    }
  }

  if ( p$aegis_dimensionality == "space-year-season") {
    ii = which(
      M$tag=="predictions" &
      M$AUID %in% res$AUID &
      M$year %in% p$yrs
    )  # filter by AUID and years in case additional data in other areas and times are used in the input data
    matchfrom = list( AUID=M$AUID[ii], year=as.character(M$year[ii]), dyear=as.character( discretize_data( M$dyear[ii], p$discretization[["dyear"]] ) ) )
    matchto   = list( AUID=res$AUID,   year=as.character(p$yrs),      dyear=as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) ) )

    if ( grepl("glm", p$carstm_modelengine) |  grepl("gam", p$carstm_modelengine) ) {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ - grep("Intercept", withsolutions) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("dyear", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      withsolutions = data.frame( matrix( unlist(strsplit( withsolutions, ":", fixed=TRUE)), ncol=3, byrow=TRUE), stringsAsFactors=FALSE )
      ws_matchfrom = list( AUID=as.character(withsolutions[,1]), year=as.character(withsolutions[,2]), dyear=as.character(withsolutions[,3])  )
      NA_mask = reformat_to_array( input =rep(1, nrow(withsolutions)), matchfrom=ws_matchfrom, matchto=matchto )
    }
  }


  if ( grepl("glm", p$carstm_modelengine) ) {

    preds = predict( fit, newdata=M[ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2


    vn =  paste(p$variabletomodel, "predicted", sep=".")
    input = preds$fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto ) * NA_mask[]
    if ( grepl( "link.*=.*log", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto ) * NA_mask[]

    vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto ) * NA_mask[]
    if ( grepl( "link.*=.*log", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto ) * NA_mask[]
    if ( grepl( "link.*=.*log", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

  }


  if ( grepl("gam", p$carstm_modelengine) ) {

    preds = predict( fit, newdata=M[ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

    vn =  paste(p$variabletomodel, "predicted", sep=".")
    input = preds$fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if ( grepl( "link.*=.*log", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_se", sep=".")
    input = preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )

    vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
    input = preds$fit - preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if ( grepl( "link.*=.*log", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

    vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
    input = preds$fit + preds$se.fit
    res[[vn]] = reformat_to_array( input =input, matchfrom=matchfrom, matchto=matchto )
    if ( grepl( "link.*=.*log", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
    if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive
  }


  if ( grepl("inla", p$carstm_modelengine) ) {

    if (exists("summary.fitted.values", fit)) {

      vn = paste( p$variabletomodel, "predicted", sep=".")
      input = fit$summary.fitted.values[ ii, "mean" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ ii, "0.025quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ ii, "0.975quant" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) res[[vn]] = exp(res[[vn]])
      if (exists("data_transformation", p) ) res[[vn]] = p$data_transformation$backward( res[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_se", sep=".")
      input = fit$summary.fitted.values[ ii, "sd" ]
      res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )

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
        if ( p$aegis_dimensionality == "space") {
          ii = which(
            M$tag=="predictions" &
            M$AUID %in% res$AUID
          )  # filter by AUID and years in case additional data in other areas and times are used in the input data
          matchfrom = list( AUID=M$AUID[ii] )
          matchto   = list( AUID=res$AUID )
        }

        if ( p$aegis_dimensionality == "space-year") {
          ii = which(
            M$tag=="predictions" &
            M$AUID %in% res$AUID &
            M$year %in% p$yrs
          )  # filter by AUID and years in case additional data in other areas and times are used in the input data
          matchfrom = list( AUID=M$AUID[ii], year=as.character(M$year[ii]) )
          matchto   = list( AUID=res$AUID, year=as.character(p$yrs)  )
        }

        if ( p$aegis_dimensionality == "space-year-season") {
          ii = which(
            M$tag=="predictions" &
            M$AUID %in% res$AUID &
            M$year %in% p$yrs
          )  # filter by AUID and years in case additional data in other areas and times are used in the input data
          matchfrom = list( AUID=M$AUID[ii], year=as.character(M$year[ii]), dyear=as.character( discretize_data( M$dyear[ii], p$discretization[["dyear"]] ) ) )
          matchto   = list( AUID=res$AUID,   year=as.character(p$yrs),      dyear=as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) ) )
        }
        vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[ii, "mean" ]
        res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
      }

      if (exists("auid", fit$summary.random)) {

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single nonspatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          kk = which(resout$type=="nonspatial")
          matchfrom = list( AUID=resout$AUID[kk]  )
          matchto   = list( AUID=res$AUID  )
          input = fit$summary.random$auid[ kk, "mean" ]
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          #  nonspatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          kk = which(resout$type=="nonspatial")
          matchfrom = list( AUID=resout$AUID[kk], year=resout$year[kk] )
          matchto   = list( AUID=res$AUID, year=as.character(p$yrs) )
          input = fit$summary.random$auid[ kk, "mean" ]
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # nonspatial at all time slices
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          kk = which(resout$type=="nonspatial")
          matchfrom = list( AUID=resout$AUID[kk], year=as.character(resout$year[kk]), dyear=as.character( discretize_data( resout$dyear[kk], p$discretization[["dyear"]] ) ) )
          matchto   = list( AUID=res$AUID,   year=as.character(p$yrs),      dyear=as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) ) )
          input = fit$summary.random$auid[ kk, "mean" ]
        }

        vn = paste( p$variabletomodel, "random_auid_nonspatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single spatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          kk = which(resout$type=="spatial")
          matchfrom = list( AUID=resout$AUID[kk]  )
          matchto   = list( AUID=res$AUID  )
          input = fit$summary.random$auid[ kk, "mean" ]  # offset structure due to bym2
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          # spatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          kk = which(resout$type=="spatial")
          matchfrom = list( AUID=resout$AUID[kk], year=resout$year[kk] )
          matchto   = list( AUID=res$AUID, year=p$yrs )
          input = fit$summary.random$auid[ kk, "mean" ]  # offset structure due to bym2
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # at every time slice
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          kk = which(resout$type=="spatial")
          matchfrom = list( AUID=resout$AUID[kk], year=as.character(resout$year[kk]), dyear=as.character( discretize_data( resout$dyear[kk], p$discretization[["dyear"]] ) ) )
          matchto   = list( AUID=res$AUID,   year=as.character(p$yrs),      dyear=as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) ) )
          input = fit$summary.random$auid[ kk, "mean" ]  # offset structure due to bym2
        }
        vn = paste( p$variabletomodel, "random_auid_spatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=matchfrom, matchto=matchto )
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

      }
    }
  }

  save( res, file=fn, compress=TRUE )
  return(res)
}
