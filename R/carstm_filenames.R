carstm_filenames = function( p,  projecttype="carstm_inputs", projectname=NULL, areal_units_fn=NULL, dropextension=FALSE ) {

  if (is.null(areal_units_fn)) {
    sppoly = areal_units( p=p )  # required by car fit
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  }

  if (projecttype=="carstm_inputs") {
    fnroot = paste(
      projectname,
      projecttype,
      areal_units_fn,
      p$variabletomodel,
      p$inputdata_spatial_discretization_planar_km,
      sep="."
    )

    if (exists("inputdata_temporal_discretization_yr", p)) {
      fnroot = paste(
        fnroot,
        round(p$inputdata_temporal_discretization_yr, 6),
        sep="."
      )
    }

    if (dropextension) return (fnroot)
    fn = paste( fnroot, "rdata", sep="." )
    return(fn)
  }


  if (projecttype=="carstm_outputs") {

    if (exists( "inputdata_spatial_discretization_planar_km", p )) {
      areal_units_fn = paste(
        areal_units_fn,
        round(p$inputdata_spatial_discretization_planar_km, 6),
        sep="_"
      )
    }

    if (exists( "inputdata_temporal_discretization_yr", p )) {
      areal_units_fn = paste(
        areal_units_fn,
        round(p$inputdata_temporal_discretization_yr, 6),
        sep="_"
      )
    }

    fnroot = paste(
      areal_units_fn,
      p$variabletomodel,
      p$carstm_modelengine,
      sep="."
    )

    if (dropextension) return (fnroot)
    fn = paste( fnroot, "rdata", sep="." )
    return(fn)
  }

}
