
additional_features_tmap = function( p, area_lines=NULL, isobaths=NULL, coastline=NULL, xlim=c(-85,-35), ylim=c(35, 65) ) {
  require(tmap)
  plot_crs = p$aegis_proj4string_planar_km

  O =  NULL
  
  if (!is.null(area_lines )) {
    # area_lines = "cfa.regions"
    O = O + tm_shape( aegis.polygons::area_lines.db( DS=area_lines, returntype="sf", project_to=plot_crs ), projection=plot_crs ) + 
      tm_lines( col="slategray", alpha=0.75, lwd=2)  
  }

  if (!is.null(isobaths)) {
    # isobaths=c(  10, 100, 200, 300, 400, 800 )
    O = O + tm_shape( aegis.bathymetry::isobath_db( depths=isobaths, project_to=plot_crs ), projection=plot_crs ) +
      tm_lines( col="slategray", alpha=0.5, lwd=0.5) 
  }

  if (!is.null(coastline)) {
    # coastline =  c("unsited states of america", "canada")
    O = O + tm_shape( st_transform( polygons_rnaturalearth(countries=coastline, xlim=xlim, ylim=ylim), st_crs(plot_crs) ), projection=plot_crs ) +
      tm_borders( col="slategray", alpha=0.5, lwd=0.75)
  }

  return(O)

  if (0) {
    o = additional_features_tmap( p=p, area_lines="cfa.regions", isobaths=c(  10, 100, 200, 300,   800 ), 
      coastline =  c("unsited states of america", "canada"), xlim=c(-85,-35), ylim=c(35, 65) )
  }
}


