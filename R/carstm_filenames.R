carstm_filenames = function( p,  returntype, areal_units_fn=NULL  ) {

  if (is.null(areal_units_fn)) {
    sppoly = areal_units( p=p )  # required by car fit
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  }

  fnroot = paste(
    areal_units_fn,
    returntype,
    p$variabletomodel,
    p$carstm_modelengine, 
    sep="|"
  )

  if ( grepl( "carstm_inputs", returntype) ) {
    outputdir = file.path( p$modeldir )  # as input data can be shared across a number of scenarios
  } else {
    outputdir = file.path( p$modeldir, p$carstm_model_label )
  }

  if (exists( "inputdata_spatial_discretization_planar_km", p )) {
    fnroot = paste( fnroot, round(p$inputdata_spatial_discretization_planar_km*1000), sep="|"  )  # m
  }

  if (exists("inputdata_temporal_discretization_yr", p)) {
    fnroot = paste( fnroot, round(p$inputdata_temporal_discretization_yr*365),  sep="|" )  # convert to days
  }

  fnroot = paste( fnroot, "rdata", sep="." )

  fnfull = file.path( outputdir, fnroot )

  return(fnfull)
}
