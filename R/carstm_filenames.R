carstm_filenames = function( p=list(), returnvalue="full_filename", ...  ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  # returntype, areal_units_fn=NULL, variabletomodel, carstm_modelengine="inla"

  if (!exists("areal_units_fn", p)) {
    sppoly = areal_units( p=p )  # required by car fit
    p$areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
  }

  if (!exists("fnroot", p)) {
    p$fnroot = paste(
      gsub( ".rdata$", "", p$areal_units_fn),
      p$returntype,
      p$variabletomodel,
      p$carstm_modelengine, 
      sep="~"
    )
    p$fnroot = paste( p$fnroot, "rdata", sep="." )
  }

  if (!exists("outputdir", p)) {
    p$outputdir = file.path( p$modeldir, p$carstm_model_label )
  }

  p$fnfull = file.path( p$outputdir, p$fnroot )

  if (returnvalue=="full_filename") return(p$fnfull)
  if (returnvalue=="file_root") return(p$fnroot)
  if (returnvalue=="output_directory") return(p$outputdir)
  if (returnvalue=="areal_units_filename") return(p$areal_units_fn)

  return(p)
}
