
carstm_parameters = function( p=NULL, redo=FALSE, ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  p$libs = c( p$libs, RLibrary ( "sp", "spdep", "rgeos", "INLA" ))  # standard CRAN libs -- geostatistical support
  p$libs = c( p$libs, RLibrary( "aegis", "bio.taxonomy", "carstm") ) #,  ) # locally developed code
  # aegis -- support for data layers
  # carstm -- support for stan bym-car models

  if (!exists("project_name", p)) stop("project_name is required")
  if (!exists("areal_units_resolution_km", p )) stop( "areal_units_resolution_km should be defined ... " ) # km
  if (!exists("areal_units_proj4string_planar_km", p )) stop( "areal_units_proj4string_planar_km should be defined ... " ) # km

  p = parameters_add_without_overwriting( p,
    areal_units_type = "lattice", #
    areal_units_overlay = "none", #
    carstm_model_label = "default",
    carstm_modelengine ="inla", # glm and gam also possible ... though not very useful
    nsims = 5000,
    boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
    mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
  )

  if ( !exists("carstm_inputdata_model_source", p))  p$carstm_inputdata_model_source = list()
  p$carstm_inputdata_model_source = parameters_add_without_overwriting( p$carstm_inputdata_model_source,
     bathymetry = "stmv",  # "stmv", "hybrid", "carstm" -- could use carstm but stmv does well with pure spatial processess
     substrate = "stmv",  # "stmv", "hybrid", "carstm" -- could use carstm but stmv does well with pure spatial processess
     temperature = "carstm",  # "stmv", "hybrid", "carstm"
     speciescomposition = "carstm" # "stmv", "hybrid", "carstm"
  )

  p = parameters_add_without_overwriting( p, discretization=list() )
  p$discretization = parameters_add_without_overwriting( p$discretization,
    z = c(0, 10, 20, 40, 80, 100, 150, 200, 250, 300, 350, 400, 500, 1000, 2000, 5000, 8000 ),
    dz = c(0.01, 0.1,  1, 2, 4, 6, 8, 10, 12 ),
    ddz = c(0.01, 0.1, 0.2, 0.4, 0.8, 1, 2, 4  ),
    substrate.grainsize = c( 0, 1, 2, 4, 8, 12, 16, 20, 32 ),
    pca1 = seq( -1, 1, by=0.1  ),
    pca2 = seq( -1, 1, by=0.1  ),
    t = seq( -4, 25, by=1  ),
    tsd = seq( 0, 25, by=1  ),
    tmin = seq( -4, 25, by=1  ),
    tmax = seq( -4, 25, by=1  ),
    degreedays = c(10, 100, 200, 400, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000),
    dyear = seq( 0, 1, by=0.1  )
  )

  p$discretization = parameters_add_without_overwriting( p$discretization,
    n.season = length(p$discretization[["dyear"]]) - 1   # used by seasonal error in inla ... not really useful
  )

  message ( p$project_name, " areal_units_resolution_km: ", p$areal_units_resolution_km)

  return(p)

}
