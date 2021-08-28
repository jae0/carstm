
carstm_model = function( p=list(), data=NULL, dimensionality=NULL,  
  sppoly =NULL, space.id = NULL, time.id = NULL, cyclic.id=NULL, areal_units_fn=NULL, DS="redo", 
  improve.hyperparam.estimates=FALSE, compress=TRUE, fn_fit=NULL, fn_res=NULL, 
   ... ) {

     if (0) {
      data=NULL
      E=NULL
      sppoly =NULL
      space.id = NULL
      time.id = NULL,
      cyclic.id=NULL,
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
    message("Loading sppoly from file ..." )
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
      if (file.exists(fn_fit)) load( fn_fit )
      if (is.null(fit)) message("carstm modelled fit not found.")
    }
    return( fit )
  }

  O = NULL
  if (DS=="carstm_modelled_summary") {  # carstm_model.*carstm_modelled
    if (!is.null(fn_res)) {
      message("Loading carstm data summary:  ", fn_res )
      if (file.exists(fn_res)) load( fn_res)
      if (is.null(O)) message("carstm summary not found.")
      return( O )
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


  if ( grepl("inla", carstm_modelengine) ) {
    out = carstm_model_inla( O=p, data=data, sppoly=sppoly, space.id=space.id, time.id = time.id,
      cyclic.id=cyclic.id, fn_fit=fn_fit, fn_res=fn_res, compress=compress, ... ) 
  }
   
  return( out )
}
