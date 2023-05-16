
carstm_model = function( p=list(), data=NULL, dimensionality=NULL,  
  sppoly =NULL, areal_units_fn=NULL, DS="redo",  
  compress=TRUE, fn_fit=NULL, fn_res=NULL, 
   ... ) {

     if (0) {
      data=NULL
      E=NULL
      sppoly =NULL
      areal_units_fn=NULL
      dimensionality=NULL
      DS="redo"
      improve.hyperparam.estimates=FALSE
      carstm_modelengine = "inla"
      compress=TRUE
      fn_fit=NULL
      fn_res=NULL
     }
     
  # compute and extract in one go esp as inla data files are too large, otherwise

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  p = parameters_add_without_overwriting( p,
    areal_units_type = "lattice", #
    areal_units_overlay = "none", #
    carstm_model_label = "default",
    carstm_modelengine ="inla" # glm and gam also possible ... though not very useful
  )

  if (!exists("dimensionality", p)) {
    if (p$aegis_dimensionality=="space") p$dimensionality = "space"  
    if (p$aegis_dimensionality=="space-year") p$dimensionality = "space-time"  
    if (p$aegis_dimensionality=="space-year-season") p$dimensionality = "space-time-cyclic"  
  }
  if (!exists("dimensionality", p)) stop("dimensionality needs to be specified")

  if (is.null(sppoly)) {
    message("sppoly not specified, loading from file using parameters provided ..." )
    sppoly = try(  areal_units( p=p ) , silent = TRUE ) # required by car fit
    message("\n")
  }

  if (!inherits( sppoly, "try-error" )) {
    if (is.null(areal_units_fn))  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  }
  
  if (exists("carstm_modelengine", p)) carstm_modelengine = p$carstm_modelengine

  if (!is.null(areal_units_fn)) {
    if (is.null(fn_fit)) fn_fit = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
    if (is.null(fn_res)) fn_res = carstm_filenames( p=p, returntype="carstm_modelled_summary", areal_units_fn=areal_units_fn )
  }

  fit = NULL
  if (DS=="carstm_modelled_fit") {
    if (!is.null(fn_fit)) {
      message("Loading carstm fit: ", fn_fit )
      if (file.exists(fn_fit)) {
        if (grepl("\\.RDS$", fn_fit)) {
          fit = readRDS(fn_fit)
        } else {
          load( fn_fit )
        }
      }
      if (is.null(fit)) message("carstm modelled fit not found.")
    }
    return( fit )
  }

  O = NULL
  if (DS=="carstm_modelled_summary") {  # carstm_model.*carstm_modelled
    if (!is.null(fn_res)) {
      # message("Loading  data summary:  ", fn_res )
      O = NULL
      if (file.exists(fn_res)) {
        if (grepl("\\.RDS$", fn_res)) {
          O = readRDS(fn_res)
        } else {
          load( fn_res )
        }
      }
      if (is.null(O)) message(" summary not found.")
      return( O )
    } else {
      fit  = NULL
      # message("Loading results from fit: ", fn_fit )
      if (file.exists(fn_fit)) {
        if (grepl("\\.RDS$", fn_fit)) {
          fit = readRDS(fn_fit)
        } else {
          load( fn_fit )
        }
      }
      if (!is.null(fit)) {
        if (exists( "results", fit)) {
          return( fit$results )
        }
      }
      message("modelled results not found. .. trying to run extraction: '' ")
      out = carstm_model_inla( O=p, data=data, sppoly=sppoly, fn_fit=fn_fit, fn_res=fn_res, compress=compress, redo_fit=FALSE, ... )
      return(out)
    }
  }


  outputdir = dirname(fn_fit)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
   
  if ( grepl("glm", carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    out = carstm_model_glm( O=p, data=data, fn_fit=fn_fit,  fn_res=fn_res, compress=compress, ... ) 
  }

  if ( grepl("gam", carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    out = carstm_model_gam( O=p, data=data, fn_fit=fn_fit,  fn_res=fn_res, compress=compress, ... ) 
  }

  if (grepl("bayesx", carstm_modelengine) ) {
    # might be useful..
  }

  if (grepl("diseasemapping", carstm_modelengine) ) {
    # might be useful..
  }


  if ( grepl("inla_old_method", carstm_modelengine) ) {
    out = carstm_model_inla_working_copy( O=p, data=data, sppoly=sppoly, fn_fit=fn_fit, fn_res=fn_res, compress=compress, ... ) 
  }

  if ( grepl("inla", carstm_modelengine) ) {
    out = carstm_model_inla( O=p, data=data, sppoly=sppoly, fn_fit=fn_fit, fn_res=fn_res, compress=compress, ... )
  }

  return( out )
}
