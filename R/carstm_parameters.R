
carstm_parameters = function( p=NULL, ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  p$libs = c( p$libs, RLibrary ( "sf",  "spdep", "rgeos", "INLA" ))  # standard CRAN libs -- geostatistical support
  p$libs = c( p$libs, RLibrary( "aegis", "carstm") ) #,  ) # locally developed code
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

  # used for lookups
  if ( !exists("carstm_inputdata_model_source", p))  p$carstm_inputdata_model_source = list()
  p$carstm_inputdata_model_source = parameters_add_without_overwriting( p$carstm_inputdata_model_source,
     bathymetry = "stmv",  # "stmv", "hybrid", "carstm" -- could use carstm but stmv does well with pure spatial processess
     substrate = "stmv",  # "stmv", "hybrid", "carstm" -- could use carstm but stmv does well with pure spatial processess
     temperature = "carstm",  # "stmv", "hybrid", "carstm"
     speciescomposition = "carstm" # "stmv", "hybrid", "carstm"
  )

  return(p)
}
