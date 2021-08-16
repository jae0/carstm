carstm_filenames = function( p=list(), returnvalue="full_filename", ...  ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  # returntype, areal_units_fn=NULL, variabletomodel, carstm_modelengine="inla"

  if (!exists("areal_units_fn", p)) {
    sppoly = areal_units( p=p )  # required by car fit
    p$areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  }

  if (!exists("fnroot", p)) {
    p$fnroot = paste(
      p$areal_units_fn,
      p$returntype,
      p$variabletomodel,
      p$carstm_modelengine, 
      sep="|"
    )
    if (exists( "inputdata_spatial_discretization_planar_km", p )) {
      p$fnroot = paste( p$fnroot, round(p$inputdata_spatial_discretization_planar_km*1000), sep="|"  )  # m
    }
    if (exists("inputdata_temporal_discretization_yr", p)) {
      p$fnroot = paste( p$fnroot, round(p$inputdata_temporal_discretization_yr*365),  sep="|" )  # convert to days
    }
    p$fnroot = paste( p$fnroot, "rdata", sep="." )
  }

  if (!exists("outputdir", p)) {
    p$outputdir = file.path( p$modeldir, p$carstm_model_label )
    if (exists("returntype", p)) { 
      if ( grepl( "carstm_inputs", p$returntype) ) p$outputdir = file.path( p$modeldir )  # as input data can be shared across a number of scenarios
    } 
  }

  p$fnfull = file.path( p$outputdir, p$fnroot )

  if (returnvalue=="full_filename") return(p$fnfull)
  if (returnvalue=="file_root") return(p$fnroot)
  if (returnvalue=="output_directory") return(p$outputdir)
  if (returnvalue=="areal_units_filename") return(p$areal_units_fn)

  return(p)
}
