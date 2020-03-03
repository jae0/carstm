
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

  fn_fit = file.path( outputdir, paste( "carstm_modelled_fit", areal_units_fns_suffix, sep=".") )
  if (DS=="carstm_modelled_fit") {
    if (file.exists(fn_fit)) {
      load( fn_fit )
      return( fit )
    }
  }


  # prediction surface
  sppoly = areal_units( p=p )  # will redo if not found

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
    } else if ( grepl( "family.*=.*binomial", p$carstm_modelcall)) {
      m = M[ j, p$variabletomodel ]
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

  return(fit)

}
