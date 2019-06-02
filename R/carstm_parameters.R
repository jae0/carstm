
carstm_parameters = function( p=NULL, redo=FALSE, ... ) {

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  p = aegis::aegis_parameters( p=p, DS="survey", yrs=p$yrs ) # basic spatial parameters

  # ---------------------
  p$libs = c( p$libs, RLibrary ( "sp", "spdep", "rgeos", "INLA" ))  # standard CRAN libs -- geostatistical support
  p$libs = c( p$libs, RLibrary( "aegis", "bio.taxonomy", "carstm") ) #,  ) # locally developed code
  # aegis -- support for data layers
  # carstm -- support for stan bym-car models

  if (!exists("polygon_source", p)) p$polygon_source = "pre2014"   # "pre2014" for older
  if (!exists("internal.crs", p)) p$internal.crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  if (!exists("internal.crs_planar", p)) p$internal.crs_planar = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km"  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  if (!exists("boundingbox", p)) p$boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)) # bounding box for plots using spplot
  if (!exists("trawlable_units", p)) p$trawlable_units = "standardtow"  # for groundfish.db

  # set up default map projection
  p = c(p, aegis::coastline_layout( p=p, redo=redo )
  p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
  # to plot a map: p = aegis::coastline_layout( p=p, plotmap=TRUE )

  return(p)

}
