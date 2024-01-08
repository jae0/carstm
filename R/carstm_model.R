
carstm_model = function( p=list(), data=NULL, sppoly =NULL, areal_units_fn=NULL, DS="redo", 
    space_id=NULL, time_id=NULL, cyclic_id=NULL, theta=NULL, carstm_directory=NULL, 
    nposteriors=NULL, posterior_simulations_to_retain=NULL,
    compress="gzip", compression_level=1, fn_fit=NULL, fn_res=NULL, debug=FALSE, 
    ... ) {

     if (0) {
      data=NULL
      sppoly =NULL
      areal_units_fn=NULL
      DS="redo"
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

 
  if (is.null(areal_units_fn))  areal_units_fn = areal_units_filename(p)
  
  if (exists("carstm_modelengine", p)) carstm_modelengine = p$carstm_modelengine

  if (!is.null(areal_units_fn)) {
    if (is.null(fn_fit)) fn_fit = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
    if (is.null(fn_res)) fn_res = carstm_filenames( p=p, returntype="carstm_modelled_summary", areal_units_fn=areal_units_fn )
  }

  if (exists("carstm_directory", p)) {
    # override default load/save locations:
    fn_fit = file.path( p$carstm_directory, basename( fn_fit) )
    fn_res = file.path( p$carstm_directory, basename( fn_res) )
  }

  if (!is.null(carstm_directory)) {
    # override default load/save locations:
    fn_fit = file.path( carstm_directory, basename( fn_fit) )
    fn_res = file.path( carstm_directory, basename( fn_res) )
  }

  print(fn_res)
  
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
    out = carstm_model_glm( O=p, data=data, fn_fit=fn_fit,  fn_res=fn_res, compress=compress, compression_level=compression_level, ... ) 
  }

  if ( grepl("gam", carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    out = carstm_model_gam( O=p, data=data, fn_fit=fn_fit,  fn_res=fn_res, compress=compress, compression_level=compression_level, ... ) 
  }


  if (grepl("julia", carstm_modelengine) ) {
    # fast ...  can call via JuliaR
  }

  if (grepl("brms", carstm_modelengine) ) {
    # fast ... 
  }


  if (grepl("stmv", carstm_modelengine) ) {
    # via carstm/stmv hybrid mode   ... 
  }

  if ( grepl("inla", carstm_modelengine) ) {
    out = carstm_model_inla( O=p, data=data, sppoly=sppoly, fn_fit=fn_fit, fn_res=fn_res, 
      compress=compress, compression_level=compression_level, 
      nposteriors=nposteriors, posterior_simulations_to_retain=posterior_simulations_to_retain,
      space_id=space_id, time_id=time_id, cyclic_id=cyclic_id, theta=theta, debug=debug, ... )
  }

  return( out )
}
