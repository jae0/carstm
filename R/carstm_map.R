
  carstm_map = function( 
    res, 
    vn, 
    xyz=NULL,
    sppoly=NULL,
    poly_match=NULL, 
    time_match=NULL, 
    spatial_domain=NULL,
    plot_crs=st_crs( projection_proj4string("lonlat_wgs84") ),
    coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ), 
    isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs  ), 
    managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs ),
    aggregate_function=mean,
    probs=c(0,0.975), 
    outformat="pdf",
    outfilename="",
    ...) {

    # carstm/aegis wrapper around spplot, will forward any additional spplot args 
    # TODO move to an sf-based plot

    require(sf)
    require(tmap)

    if (!is.null(xyz)) {
      # if xyz is given then then are point data to be aggregated into polygons
      sppoly$AUID = as.character(1:nrow(sppoly))

      if ("data.frame" %in% class(xyz)) {
        xyz = st_as_sf( xyz, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
      }

      xyz = st_transform(xyz, plot_crs )
      xyz$AUID = st_points_in_polygons( pts=xyz, polys=sppoly[, "AUID"], varname="AUID" )

      toplot = xyz[[vn]]   
      auid = xyz[["AUID"]] 
      oo = tapply( toplot, auid, aggregate_function, na.rm=TRUE )
      sppoly[[vn]] = NA
      sppoly[[vn]][ match( names(oo) , sppoly$AUID ) ] = oo


    } else {

      if (is.null(sppoly)) sppoly = res$sppoly

      # first index is spatial strata
      data_dimensionality = length( dim(res[[vn]]) )

      if (is.null(poly_match)) poly_match = match( res$AUID, sppoly[["AUID"]] )  # should match exactly but in case a subset is sent as sppoly

      if (data_dimensionality==1) {
      
        toplot = res[[vn]] [ poly_match ]  # year only
      
      } else if (data_dimensionality==2) {

          toplot = res[[vn]] [ poly_match, time_match[[1]] ]  # year only

      } else if (data_dimensionality==3) {

          toplot = res[[vn]] [ poly_match, time_match[[1]], time_match[[2]] ] # year/subyear

      }
      sppoly[,vn] = toplot
      sppoly = st_transform( sppoly, crs=st_crs(plot_crs) )

    }


    er = quantile( sppoly[[vn]], probs=probs, na.rm=TRUE )
    datarange = seq( er[1], er[2], length.out=7) 
    sppoly[[vn]][ which(sppoly[[vn]] < er[1]) ] = er[1] # set levels higher than max datarange to max datarange
    sppoly[[vn]][ which(sppoly[[vn]] > er[2]) ] = er[2] # set levels higher than max datarange to max datarange

    tmap_mode("plot")
    
    maintitle = ifelse ( exists("main"), main, vn ) 
    
    
    o = tm_shape( sppoly, projection=plot_crs ) +
      tm_polygons(
        vn,
        style = "cont",
        breaks = datarange,
        title= maintitle,
        border.col = NULL,
        colorNA = NULL,
        constrast=c(0,0.6),
        showNA=FALSE,
        lwd = 0.5, 
        palette = "YlOrRd",
        border.alpha = 0.5,
        legend.is.portrait = FALSE ) +
    tm_shape( coastline, projection=plot_crs ) +
      tm_polygons( col="grey80" ) +
    tm_shape( isobaths, projection=plot_crs ) +
      tm_lines( col="lightgray", alpha=0.5) +
    tm_shape( managementlines, projection=plot_crs ) +
      tm_lines( col="grey20", alpha=0.75, lwd=2) +

    tm_compass( position=c( "right", "top")) + 
    tm_scale_bar( position=c("left", "bottom" ) ) +
    tm_layout( frame=FALSE, legend.text.size= 0.7 )
    
    print(o)
  
    if ( outfilename !="" ) {
      if (outformat=="pdf") pdf( file=outfilename, width=9, height=7, bg='white', pointsize=12 )
      if (outformat=="svg") svg( filename=outfilename, width=9, height=7, bg='white', pointsize=12   )
      if (outformat=="png") png( filename=outfilename, width=3072, height=2304, pointsize=40, res=300 )
      print(o)

      dev.off()
      print(outfilename)
    }
    
    return(o)
  }

