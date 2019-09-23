
carstm_parameters = function( p=NULL, redo=FALSE, ... ) {

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  p = aegis.survey::survey_parameters( p=p, yrs=p$yrs ) # basic spatial parameters

  # ---------------------
  p$libs = c( p$libs, RLibrary ( "sp", "spdep", "rgeos", "INLA" ))  # standard CRAN libs -- geostatistical support
  p$libs = c( p$libs, RLibrary( "aegis", "bio.taxonomy", "carstm") ) #,  ) # locally developed code
  # aegis -- support for data layers
  # carstm -- support for stan bym-car models

  if (!exists("polygon_source", p)) p$polygon_source = "pre2014"   # "pre2014" for older
  if (!exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  if (!exists("boundingbox", p)) p$boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)) # bounding box for plots using spplot
  if (!exists("trawlable_units", p)) p$trawlable_units = "standardtow"  # for groundfish.db

  # set up default map projection
  p = c(p, aegis.coastline::coastline_layout( p=p, redo=redo ) )
  p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
  # to plot a map: p = aegis.coastline::coastline_layout( p=p, plotmap=TRUE )

  return(p)

}
