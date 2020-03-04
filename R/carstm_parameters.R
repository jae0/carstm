
carstm_parameters = function( p=NULL, redo=FALSE, ... ) {

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  p$libs = c( p$libs, RLibrary ( "sp", "spdep", "rgeos", "INLA" ))  # standard CRAN libs -- geostatistical support
  p$libs = c( p$libs, RLibrary( "aegis", "bio.taxonomy", "carstm") ) #,  ) # locally developed code
  # aegis -- support for data layers
  # carstm -- support for stan bym-car models

  if ( !exists("project_name", p)) stop("project_name is required")

  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  if (!exists("areal_units_source", p )) p$areal_units_source = "lattice" #
  if (!exists("areal_units_constraint", p )) p$areal_units_constraint = "none" #
  if (!exists("areal_units_overlay", p )) p$areal_units_overlay = "none" #
  if (!exists("areal_units_resolution_km", p )) stop( "areal_units_resolution_km should be defined ... " ) # km
  if (!exists("areal_units_proj4string_planar_km", p )) stop( "areal_units_proj4string_planar_km should be defined ... " ) # km
  if (!exists("timeperiod", p) )  p$timeperiod="default"

  if ( !exists("carstm_model_label", p) ) p$carstm_model_label = "production"  # default is to choose production results

  if (!exists("areal_units_fn", p) ) p$areal_units_fn = paste(
    p$spatial_domain,
    paste0(p$areal_units_overlay, collapse="_"),
    p$areal_units_resolution_km,
    p$areal_units_source,
    p$areal_units_constraint,
    p$timeperiod,
    sep="_"
  )

  if (!exists("nsims", p) ) p$nsims = 5000

  if (!exists("discretization", p) ) p$discretization = list()
  if (!exists("z", p$discretization) ) p$discretization$z = c(0, 10, 20, 40, 80, 100, 150, 200, 250, 300, 350, 400, 500, 1000, 2000, 5000 )  # depth cut points
  # if (!exists("z", p$discretization) ) p$discretization$z = c(2.5, 5, 10, 20, 40, 80, 160, 320, 640 )
  if (!exists("dz", p$discretization) ) p$discretization$dz = c(0.01, 0.1,  1, 2, 4, 6, 8, 10, 12 )  # slope cut points
  if (!exists("ddz", p$discretization) ) p$discretization$ddz = c(0.01, 0.1, 0.2, 0.4, 0.8, 1, 2, 4  )  # slope cut points
  if (!exists("substrate.grainsize", p$discretization) ) p$discretization$substrate.grainsize = c( 0, 1, 2, 4, 8, 12, 16, 20, 32 )

  if (!exists("pca1", p$discretization) ) p$discretization$pca1 = seq( -1, 1, by=0.1  )
  if (!exists("pca2", p$discretization) ) p$discretization$pca2 = seq( -1, 1, by=0.1  )

  if (!exists("t", p$discretization) ) p$discretization$t = seq( -4, 25, by=1  )
  if (!exists("tsd", p$discretization) ) p$discretization$tsd = seq( 0, 25, by=1  )
  if (!exists("tmin", p$discretization) ) p$discretization$tmin = seq( -4, 25, by=1  )
  if (!exists("tmax", p$discretization) ) p$discretization$tmax = seq( -4, 25, by=1  )
  if (!exists("degreedays", p$discretization) ) p$discretization$degreedays = c(10, 100, 200, 400, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000)

  if (!exists("dyear", p$discretization) ) p$discretization$dyear = seq( 0, 1, by=0.1  )

  p$n.season = length(p$discretization[["dyear"]]) - 1   # used by seasonal error in inla ... not really useful

  if (!exists("carstm_inputs_aggregated", p )) p$carstm_inputs_aggregated = FALSE  # default is to use all daat .. unless density is too thigh such as bathy and substrate
  if (!exists("inla_num.threads", p )) p$inla_num.threads= 1
  if (!exists("inla_blas.num.threads", p )) p$inla_blas.num.threads=1

  # set up default map projection
  if (!exists("boundingbox", p)) p$boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)) # bounding box for plots using spplot
  if (!exists("coastLayout", p)) p = c(p, aegis.coastline::coastline_layout( p=p ))
  if (!exists("mypalette", p)) p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")

  return(p)

}
