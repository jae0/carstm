
carstm_model = function( p=list(), data=NULL, sppoly =NULL, areal_units_fn=NULL, DS="redo", 
    space_id=NULL, time_id=NULL, cyclic_id=NULL, theta=NULL, carstm_directory=NULL, 
    toget=NULL, nposteriors=NULL, posterior_simulations_to_retain=NULL,
    compress="qs-preset", qs_preset="high", compression_level=1, fn_fit=NULL, fn_res=NULL, debug=FALSE, 
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

  # print(fn_res)
  
  fit = NULL
  if (DS=="carstm_modelled_fit") {
    if (!is.null(fn_fit)) {
      message("Loading carstm fit: ", fn_fit )
      if (file.exists(fn_fit)) {
        if (grepl("\\.RDS$", fn_fit)) {
          fit = aegis::read_write_fast(fn_fit)
        } else {
          load( fn_fit )
        }
      }
      if (is.null(fit)) message("carstm modelled fit not found.")
    }
    return( fit )
  }

  O = NULL
  if (DS=="carstm_modelled_summary") {   

    if (!is.null(fn_res)) {
      # message("Loading  data summary:  ", fn_res )
      O = NULL
      if (file.exists(fn_res)) {
        if (grepl("\\.RDS$", fn_res)) {
          O = aegis::read_write_fast(fn_res)
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
          fit = aegis::read_write_fast(fn_fit)
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


  if (DS=="carstm_modelinfo") {   
    fn_modelinfo = gsub( "_fit~", "_modelinfo~", fn_fit, fixed=TRUE )
    O = read_write_fast( file=fn_modelinfo )
    return( O )
  }

  if (DS=="carstm_marginals") {   
    fn_marginals = gsub( "_fit~", "_marginals~", fn_fit, fixed=TRUE )
    Omarginals = read_write_fast( file=fn_marginals )
    return( Omarginals )
  }

  if (DS=="carstm_randomeffects") {   
    fn_randomeffects = gsub( "_fit~", "_randomeffects~", fn_fit, fixed=TRUE )
    Orandom = read_write_fast( file=fn_randomeffects ) 
    return( Orandom )  
  }

  if (DS=="carstm_predictions") {   
    fn_preds = gsub( "_fit~", "_predictions~",  fn_fit, fixed=TRUE )
    Opredictions = read_write_fast( file=fn_preds ) 
    return( Opredictions)  
  }

  if (DS=="carstm_samples") {   
      fn_samples = gsub( "_fit~", "_samples~", fn_fit, fixed=TRUE )
      Osamples = read_write_fast(  file=fn_samples ) 
      return( Osamples )
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
      compress=compress, compression_level=compression_level, qs_preset=qs_preset,
      toget=toget,
      nposteriors=nposteriors, posterior_simulations_to_retain=posterior_simulations_to_retain,
      space_id=space_id, time_id=time_id, cyclic_id=cyclic_id, theta=theta, debug=debug, ... )
  }

  return( out )
}
