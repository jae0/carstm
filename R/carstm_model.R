
carstm_model = function( p=list(), data=NULL, sppoly =NULL, areal_units_fn=NULL, DS=NULL, 
    space_id=NULL, time_id=NULL, cyclic_id=NULL, theta=NULL, carstm_directory=NULL, 
    toget=NULL, nposteriors=NULL, posterior_simulations_to_retain=NULL,
    compress="", compression_level=9, fn_fit=NULL, debug=FALSE, 
    ... ) {

     if (0) {
      data=NULL
      sppoly =NULL
      areal_units_fn=NULL
      DS=NULL
      carstm_modelengine = "inla"
      compress=""
      fn_fit=NULL
      
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
  }

  if (exists("carstm_directory", p)) {
    # override default load/save locations:
    fn_fit = file.path( p$carstm_directory, basename( fn_fit) )
  }

  if (!is.null(carstm_directory)) {
    # override default load/save locations:
    fn_fit = file.path( carstm_directory, basename( fn_fit) )
  }
 
  outputdir = dirname(fn_fit)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
   
  if ( grepl("glm", carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    out = carstm_model_glm( O=p, DS=DS, data=data, fn_fit=fn_fit,  
      compress=compress, compression_level=compression_level, ... ) 
  }

  if ( grepl("gam", carstm_modelengine) ) {
    # not a CAR but for comparison with no spatial random effect model
    out = carstm_model_gam( O=p, DS=DS, data=data, fn_fit=fn_fit,  
      compress=compress, compression_level=compression_level, ... ) 
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

    out = carstm_model_inla( DS=DS, O=p, data=data, sppoly=sppoly, fn_fit=fn_fit, 
      compress=compress, compression_level=compression_level,       toget=toget,
      nposteriors=nposteriors, posterior_simulations_to_retain=posterior_simulations_to_retain,
      space_id=space_id, time_id=time_id, cyclic_id=cyclic_id, theta=theta, debug=debug, ... )
  }

  return( out )
}
