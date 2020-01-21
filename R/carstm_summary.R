
carstm_summary = function( p=NULL, operation="load_RES", wgts=1, extrapolation_limit=NA, extrapolation_replacement=NA, ... ) {

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
  fn     = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "aggregated_timeseries",  "rdata", sep="." ), sep="." ) )
  fn_no = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "space_timeseries_number",  "rdata", sep="." ), sep="." ) )
  fn_bio = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "space_timeseries_biomass",  "rdata", sep="." ), sep="." ) )
  fn_pa = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "space_timeseries_pa",  "rdata", sep="." ), sep="." ) )


  if ( operation=="load_timeseries" ) {
    RES = NA
    if (file.exists(fn)) load( fn)
    return( RES )
  }

  if ( operation=="load_spacetime_number" ) {
    nums = NA
    if (file.exists(fn_no)) load( fn_no )
    return( nums )
  }

  if ( operation=="load_spacetime_biomass" ) {
    biom = NA
    if (file.exists(fn_bio)) load( fn_bio )
    return( biom )
  }

  if ( operation %in% c("load_spacetime_pa", "load_spacetime") ) {
    pa = NA
    if (file.exists(fn_pa)) load( fn_pa )
    return( pa )
  }



  if (operation=="compute") {
    # construct meanweights matrix used to convert number to weight
    M = snowcrab_carstm( p=p, DS="carstm_inputs" )
    M$yr = M$year  # req for meanweights

    sppoly = areal_units( p=p )

    res = carstm_model( p=p, DS="carstm_modelled", carstm_model_label=p$carstm_model_label ) # to load currently saved res

    X = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
    X[!is.finite(X)] = NA

    if (p$selection$type %in% c("presence_absence") ) {
      pa = X
      save( pa, file=fn_pa, compress=TRUE )
      RES = data.frame( yrs = p$yrs )
      # sa weighted average prob habitat
      RES$cfaall    = colSums( pa * sppoly$au_sa_km2/ sum(sppoly$au_sa_km2), na.rm=TRUE )
      RES$cfanorth  = colSums( pa * sppoly$cfanorth_surfacearea/ sum(sppoly$au_sa_km2), na.rm=TRUE )
      RES$cfasouth  = colSums( pa * sppoly$cfasouth_surfacearea/ sum(sppoly$au_sa_km2), na.rm=TRUE )
      RES$cfa23     = colSums( pa * sppoly$cfa23_surfacearea/ sum(sppoly$au_sa_km2), na.rm=TRUE )
      RES$cfa24     = colSums( pa * sppoly$cfa24_surfacearea/ sum(sppoly$au_sa_km2), na.rm=TRUE )
      RES$cfa4x     = colSums( pa * sppoly$cfa4x_surfacearea/ sum(sppoly$au_sa_km2), na.rm=TRUE )
      save( RES, file=fn, compress=TRUE )
    }

    if (p$selection$type %in% c("biomass", "number") ) {
      if (p$selection$type == "biomass") {
        biom = X
        if (is.na(extrapolation_limit)) extrapolation_limit = max(M$totwgt/M$data_offset, na.rm=T) # 28921.8426
        uu = which( biom > extrapolation_limit )
        if (length(uu) > 0 ) {
          if (is.character(extrapolation_replacement)) if (extrapolation_replacement=="extrapolation_limit" ) extrapolation_replacement = extrapolation_limit
          biom[ uu] = extrapolation_replacement
          warning("\n Extreme-valued predictions were found, capping them to max observed rates .. \n you might want to have more informed priors, or otherwise set extrapolation_replacement=NA to replacement value \n")
        }
        qnt = quantile( biom, probs=0.99, na.rm=TRUE)
        biom[biom > qnt] = qnt
        nums = biom / wgts
      }


      if (p$selection$type == "number") {
        nums = X
        if (is.na(extrapolation_limit)) extrapolation_limit = max(M$totno/M$data_offset, na.rm=T) # 17301.5199
        uu = which( nums > extrapolation_limit )
        if (length(uu) > 0 ) {
          if (is.character(extrapolation_replacement)) if (extrapolation_replacement=="extrapolation_limit" ) extrapolation_replacement = extrapolation_limit
          nums[ uu] = extrapolation_replacement
          warning("\n Extreme-valued predictions were found, capping them to max observed rates .. \n you might want to have more informed priors, or otherwise set extrapolation=NA to replacement value \n")
        }
        qnt = quantile( nums, probs=0.99, na.rm=TRUE)
        nums[nums > qnt] = qnt
        biom = nums * wgts
      }

      save( biom, file=fn_bio, compress=TRUE )
      save( nums, file=fn_no, compress=TRUE )

      # {no, kg} /km^2 -> {kn, kt}/ km^2
      RES = data.frame( yrs = p$yrs )
      RES$cfaall    = colSums( biom * sppoly$au_sa_km2/ 10^6, na.rm=TRUE )
      RES$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea/ 10^6, na.rm=TRUE )
      RES$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea/ 10^6, na.rm=TRUE )
      RES$cfa23     = colSums( biom * sppoly$cfa23_surfacearea/ 10^6, na.rm=TRUE )
      RES$cfa24     = colSums( biom * sppoly$cfa24_surfacearea/ 10^6, na.rm=TRUE )
      RES$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea/ 10^6, na.rm=TRUE )
      save( RES, file=fn, compress=TRUE )
    }

    return( fn )
  }

}

