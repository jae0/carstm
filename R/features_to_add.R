
features_to_add = function( p, area_lines=NULL, isobaths=c(100), coastline=c("canada", "us"), xlim=c(-85,-35), ylim=c(35, 65), target="ggplot" ) {

  plot_crs = st_crs(p$aegis_proj4string_planar_km)
  
  if (target=="tmap") {

    require(tmap)
    O =  NULL
    
    if (!is.null(area_lines )) {
      data = aegis.polygons::area_lines.db( DS=area_lines, returntype="sf", project_to=plot_crs )
      O = O + tm_shape( data=data, projection=plot_crs ) + tm_lines( col="gray40", alpha=0.75, lwd=2)  
    }

    if (!is.null(isobaths)) {
      data = aegis.bathymetry::isobath_db( depths=isobaths, project_to=plot_crs )
      O = O + tm_shape( data=data, projection=plot_crs ) + tm_lines( col="gray80", alpha=0.4, lwd=0.5) 
    }

    if (!is.null(coastline)) {
      # coastline =  c("unsited states of america", "canada")
      data = st_transform( polygons_rnaturalearth(countries=coastline, xlim=xlim, ylim=ylim), st_crs(plot_crs) )
      O = O + tm_shape( data, projection=plot_crs ) +  tm_borders( col="slategray", alpha=0.75, lwd=0.75)
    }

    return(O)
  }


  if (target=="ggplot" ) {

    require(ggplot2)
    O = ggplot()  
    
    if (!is.null(coastline)) {
      # coastline =  c("us", "canada")s
      data =  st_transform( polygons_rnaturalearth(countries=coastline, xlim=xlim, ylim=ylim), st_crs(plot_crs) )  
      O = O + geom_sf( data=data, aes(alpha=0.6), colour="slategray", lwd=0.6, fill=NA ) + guides(colour = "none",lwd = "none", alpha = "none", fill="none") # + theme(plot.background=element_blank())
    }

    if (!is.null(isobaths)) {
      data = aegis.bathymetry::isobath_db( depths=isobaths, project_to=plot_crs )
      O = O + geom_sf( data=data,  aes(alpha=0.8),  colour="darkgray", lwd=0.6) + guides(colour = "none",lwd = "none", alpha = "none")
    }

    if (!is.null(area_lines )) {
      # area_lines = "cfa.regions"
      data = aegis.polygons::area_lines.db( DS=area_lines, returntype="sf", project_to=plot_crs )
      O = O + geom_sf(data=data,  aes(alpha=0.9), colour="darkgray", lwd=2.5)  + guides(colour = "none", lwd = "none", alpha = "none")
    }

    O = O + theme_void()

    return(O)

  }

  if (0) {
  
    o = features_to_add( 
      p=p, 
      area_lines="cfa.regions", 
      isobaths=c(100, 200, 300, 400), 
      coastline =  c("canada", "us"), 
      xlim=c(-85,-35), 
      ylim=c(35, 65),
      target="ggplot" 
    )
  
  }
}


