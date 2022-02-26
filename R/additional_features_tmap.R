
additional_features_tmap = function( p, area_lines=NULL, isobaths=NULL, coastline=NULL ) {
  require(tmap)
  plot_crs = p$aegis_proj4string_planar_km

  O =  NULL
  
  if (!is.null(area_lines )) {
    # area_lines = "cfa.regions"
    O = O + tm_shape( aegis.polygons::area_lines.db( DS=area_lines, returntype="sf", project_to=plot_crs ), projection=plot_crs ) + 
      tm_lines( col="slategray", alpha=0.75, lwd=2)  
  }

  if (!is.null(isobaths)) {
    # isobaths=c( seq(0, 400, by=50), 1000)
    O = O + tm_shape( aegis.bathymetry::isobath_db( depths=isobaths, project_to=plot_crs ), projection=plot_crs ) +
      tm_lines( col="slategray", alpha=0.5, lwd=0.5) 
  }

  if (!is.null(coastline)) {
    # coastline =  c("united states of america", "canada"), xlim=c(-70,-59), ylim=c(40, 50)
    O = O + tm_shape( st_transform( polygons_rnaturalearth(countries=coastline), st_crs(plot_crs) ), projection=plot_crs ) +
      tm_polygons( col="darkslategray", alpha=0.9, lwd=2)
  }

  return(O)
}
