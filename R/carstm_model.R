
carstm_model = function( p, M=NULL, E=NULL, sppoly=NULL, region.id=NULL, areal_units_fn=NULL, DS="redo", improve.hyperparam.estimates=FALSE, compression_level=6, fn_fit=NULL, fn_res=NULL,
   ... ) {

  # compute and extract in one go esp as inla data files are too large, otherwise

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  if (is.null(sppoly)) {
    message("Loading sppoly from file ..." )
    sppoly = try(  areal_units( p=p ) , silent = TRUE ) # required by car fit
    message("\n")
  }

  if (!inherits( sppoly, "try-error" )) {
    if (is.null(areal_units_fn))  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  }
  
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
  
   
  if ( grepl("glm", p$carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    O = carstm_model_glm( p=p, M=M, fn_fit=fn_fit,  fn_res=fn_res, compression_level=compression_level, ... ) 
  }

  if ( grepl("gam", p$carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    O = carstm_model_gam( p=p, M=M, fn_fit=fn_fit,  fn_res=fn_res, compression_level=compression_level, ... ) 
  }

  if (grepl("bayesx", p$carstm_modelengine) ) {
    # might be useful..
  }


  if ( grepl("inla", p$carstm_modelengine) ) {
    O = carstm_model_inla( p=p, M=M, E=E, sppoly=sppoly, region.id=region.id, fn_fit=fn_fit, fn_res=fn_res, compression_level=compression_level, ... ) 
      #  print(O[["summary"]])
  }
   
  return( O )
}
