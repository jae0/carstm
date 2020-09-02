
carstm_model = function( p, M=NULL, DS="redo", ... ) {

  p = parameters_control(p, list(...), control="add") # add passed args to parameter list, priority to args

  sppoly = areal_units( p=p )  # required by car fit
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

  aufns = carstm_filenames( p=p, projecttype="carstm_outputs", areal_units_fn=areal_units_fn )

  outputdir = file.path(p$modeldir, p$carstm_model_label)

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  fn_fit = file.path( outputdir, paste( "carstm_modelled_fit", aufns, sep=".") )
  if (DS=="carstm_modelled_fit") {
    if (file.exists(fn_fit)) {
      load( fn_fit )
      return( fit )
    }
  }

  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )

  if (exists("data_transformation", p)) M[, p$variabletomodel]  = p$data_transformation$forward( M[, p$variabletomodel] ) # make all positive

  mrange = NULL
  # get hyper param scalings
  if ( grepl("inla", p$carstm_modelengine) ) {
    # hyperparms
    j = which( is.finite(M[,p$variabletomodel]) )
    mrange = range( M[ j, p$variabletomodel ]  )  # on data scale not internal
    if ( grepl( "family.*=.*lognormal", p$carstm_modelcall)) {
      m = log( M[ j, p$variabletomodel ])
    } else if ( grepl( "family.*=.*poisson", p$carstm_modelcall)) {
      m = log( M[ j, p$variabletomodel ] / M[ j, "data_offset" ]  )
      mrange = range( M[ j, p$variabletomodel ]/ M[ j, "data_offset" ]  )  # on data scale not internal
      mrange = mrange * median(M[ M$tag=="predictions", "data_offset" ] )
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

  res = carstm_summary( p=p, operation="compute", fit=fit, M=M, sppoly=sppoly, mrange=mrange )

  if (DS!="carstm_modelled_fit") fit = fn_fit

  return(fit)
}
