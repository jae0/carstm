carstm_filenames = function( p=list(), returnvalue="full_filename", fn=NULL, ...  ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  # returntype, areal_units_fn=NULL, variabletomodel, carstm_modelengine="inla"
  args = list(...)

  if (!exists("areal_units_fn", p)) {
    if (exists("sppoly", args)) {
      sppoly = args[["sppoly"]]
      p$areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
    } else {
      p$areal_units_fn = areal_units_filename(p)
    }
  }  


  if (!exists("fnroot", p)) {
    
    fnrt = gsub( "\\.rdata$|\\.rdz$", "", p$areal_units_fn)

    p$fnroot = paste(
      fnrt,
      p$returntype,
      p$variabletomodel,
      p$carstm_modelengine, 
      sep="~"
    )
    p$fnroot = paste( p$fnroot, "rdz", sep="." )
  }

  if (!exists("outputdir", p)) {
    p$outputdir = file.path( p$modeldir, p$carstm_model_label )
  }

  p$fnfull = file.path( p$outputdir, p$fnroot )

  if (returnvalue=="full_filename") return(p$fnfull)
  if (returnvalue=="file_root") return(p$fnroot)
  if (returnvalue=="output_directory") return(p$outputdir)
  if (returnvalue=="areal_units_filename") return(p$areal_units_fn)
  if (returnvalue=="filename") {
      aufns = carstm_filenames( p=p, returntype="modelled_fit", areal_units_fn=p$areal_units_fn )
      # same file naming as in carstm ..
      outputdir = dirname( aufns )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
      outfn = paste( gsub("\\.rdata$|\\.rdz$", "", aufns), fn, "rdz", sep="." )
    return( outfn )   
  }  
  # so if returntype is not amongst the above, it just gets added to the fnroot and returned

  return(p)
}
