
  carstm_map = function( 
    res=NULL, 
    vn=NULL, 
    space = "space",
    time= "time",
    season="season",
    vnstat="mean",
    tmap_mode="view",
    xyz=NULL,
    sppoly=NULL,
    poly_match=NULL, 
    time_match=NULL, 
    spatial_domain=NULL,
    plot_crs=NULL,
    coastline=NULL, 
    isobaths=NULL, 
    managementlines = NULL,
    additional_polygons=NULL,
    digits=3,
    aggregate_function=mean,
    probs=c(0,0.975), 
    outformat="png",
    outfilename="",
    ...) {


    ellps = list(...)

    require(sf)
    require(tmap)

    if (is.null(plot_crs)) plot_crs = st_crs( sppoly )
  
    # tmap_save options:
    width=1000
    height=800
    asp=0
    if (exists("width", ellps) )   width = ellps[["width"]]
    if (exists("height", ellps) )   height = ellps[["height"]]
    if (exists("asp", ellps) )   asp = ellps[["asp"]]

    # tmap plotting options:
    bbox = NULL
    id = row.names(sppoly)
    breaks = NULL
    style = "cont"; 
    palette =  "Spectral"
    title =  "" 
    showNA =  FALSE 
    lwd =  0.25 
    border.alpha =  0.95
    alpha =   0.9
    legend.is.portrait = FALSE   
    compass_position = c( 0.9, 0 )   
    compass_north = 0
    scale_bar_position = c( "left", "top" )
    legend_position = c(0.55, 0.8)
    legend_scale = 1.0
    legend_title.size=1.4
    legend_text.size=1.0
    legend.width=1

    if (exists("bbox", ellps) )   bbox = ellps[["bbox"]]
    if (exists("id", ellps) )   id = ellps[["id"]]
 
    if (exists("breaks", ellps) )   breaks = ellps[["breaks"]]
    if (exists("style", ellps) )   style = ellps[["style"]]
    if (exists("palette", ellps) )   palette = ellps[["palette"]]
    if (exists("title", ellps) )   title = ellps[["title"]]
    if (exists("showNA", ellps) )   showNA = ellps[["showNA"]]
    if (exists("lwd", ellps) )   lwd = ellps[["lwd"]]
    if (exists("border.alpha", ellps) )  border.alpha = ellps[["border.alpha"]]
    if (exists("alpha", ellps) )   alpha = ellps[["alpha"]]
    if (exists("legend.is.portrait", ellps) )   legend.is.portrait = ellps[["legend.is.portrait"]]
    if (exists("compass_position", ellps) )   compass_position = ellps[["compass_position"]]
    if (exists("compass_north", ellps) ) compass_north = ellps[["compass_north"]]
    if (exists("scale_bar_position", ellps) )   scale_bar_position = ellps[["scale_bar_position"]]
    if (exists("legend_position", ellps) )   legend_position = ellps[["legend_position"]]
    if (exists("legend_scale", ellps) )   legend_scale = ellps[["legend_scale"]]
    if (exists("legend_title.size", ellps) )  legend_title.size = ellps[["legend_title.size"]]
    if (exists("legend_text.size", ellps) )    legend_text.size = ellps[["legend_text.size"]]
    if (exists("legend.width", ellps) )   legend.width = ellps[["legend.width"]]
 
  
    if (class(vn))
    if (!is.null(xyz)) {
      # if xyz is given then then are point data to be aggregated into polygons
      sppoly$AUID = as.character(1:nrow(sppoly))

      if ("data.frame" %in% class(xyz)) {
        xyz = st_as_sf( xyz, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
      }

      xyz = st_transform(xyz, plot_crs )
      xyz[["AUID"]] = st_points_in_polygons( pts=xyz, polys=sppoly[, "AUID"], varname="AUID" )

      if ( length(vn) ==1) {
        toplot = xyz[[vn]]
      } else if ( length(vn) ==2) {
        toplot = xyz[[vn[1]]][[vn[2]]]   
      } else if ( length(vn) ==3) {
        toplot = xyz[[vn[1]]][[vn[2]]][[vn[3]]]   
      }
      oo = tapply( toplot, xyz[["AUID"]], aggregate_function, na.rm=TRUE )
      sppoly[[vn]] = NA
      sppoly[[vn]][ match( names(oo) , sppoly$AUID ) ] = oo


    } else {

      if (!is.null(res)) {
        if (is.null(sppoly)) sppoly = res$sppoly

        # first index is spatial strata
        if ( length(vn) ==1) {
          pdata = res[[vn]]
        } else if ( length(vn) ==2) {
          pdata = res[[vn[1]]][[vn[2]]]   
        } else if ( length(vn) ==3) {
          pdata = res[[vn[1]]][[vn[2]]][[vn[3]]]   
        }

        data_dimensionality = length( dim(pdata) )

        if (is.null(poly_match)) poly_match = match( res$space, sppoly[["AUID"]] )  # should match exactly but in case a subset is sent as sppoly

        if (data_dimensionality==1) {
          toplot = pdata[ poly_match, vnstat ]  # year only
        } else if (data_dimensionality==2) {
          toplot = pdata[ poly_match, time_match[[1]], vnstat ]  # year only
        } else if (data_dimensionality==3) {
          toplot = pdata[ poly_match, time_match[[1]], time_match[[2]], vnstat ] # year/subyear
        }
        sppoly[,vn] = toplot
      }  # else assume sppoly already has vn in it

      sppoly = st_transform( sppoly, crs=st_crs(plot_crs) )
    }


    ellps = list(...)

    if (is.null(coastline)) coastline = aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ) 
    if (is.null(isobaths)) isobaths = aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400), project_to=plot_crs  )

    if  ( exists("breaks", ellps)) {
      datarange = ellps[["breaks"]] 
      er = range(datarange)
    } else{ 
      er = quantile( sppoly[[vn]], probs=probs, na.rm=TRUE )
      datarange = signif( seq( er[1], er[2], length.out=7), 2) 
    }
    
    sppoly[[vn]][ which(sppoly[[vn]] < er[1]) ] = er[1] # set levels higher than max datarange to max datarange
    sppoly[[vn]][ which(sppoly[[vn]] > er[2]) ] = er[2] # set levels higher than max datarange to max datarange

    tmap_mode("plot")

    if ( !is.null(managementlines) ) {  
      o = tm_shape( sppoly, projection=plot_crs ) +
        tm_polygons(
          vn,
          style = ifelse ( exists("style", ellps), ellps[["style"]], "cont" ) ,
          breaks = datarange ,
          title= ifelse ( exists("main", ellps), ellps[["main"]], vn ) ,
          border.col = NULL,
          colorNA = NULL,
          showNA=FALSE,
          lwd = 0.5, 
          palette = ifelse ( exists("palette", ellps), ellps[["palette"]], "YlOrRd"),
          border.alpha = 0.5,
          legend.is.portrait = FALSE ) +
      tm_shape( coastline, projection=plot_crs ) +
        tm_polygons( col="grey80" ) +
      tm_shape( isobaths, projection=plot_crs ) +
        tm_lines( col="lightgray", alpha=0.6) +

      tm_shape( managementlines, projection=plot_crs ) +
          tm_lines( col="grey40", alpha=0.6, lwd=2)  +

      tm_compass( position=c( "right", "top")) + 
      tm_scale_bar( position=c("right", "bottom" ), width=0.2, text.size=0.7) +
      tm_legend( position=c("left", "top") ,  frame=TRUE, scale = 1 , title.size=1.5, text.size=0.80, legend.width=0.75) +
      tm_layout( frame=FALSE )
    
    } else {

      o = tm_shape( sppoly, projection=plot_crs ) +
        tm_polygons(
          vn,
          style = ifelse ( exists("style", ellps), ellps[["style"]], "cont" ) ,
          breaks = datarange ,
          title= ifelse ( exists("main", ellps), ellps[["main"]], vn ) ,
          border.col = NULL,
          colorNA = NULL,
          showNA=FALSE,
          lwd = 0.5, 
          palette = ifelse ( exists("palette", ellps), ellps[["palette"]], "YlOrRd"),
          border.alpha = 0.5,
          legend.is.portrait = FALSE ) +
      tm_shape( coastline, projection=plot_crs ) +
        tm_polygons( col="grey80" ) +
      tm_shape( isobaths, projection=plot_crs ) +
        tm_lines( col="lightgray", alpha=0.75) +
      tm_compass( position=c( "right", "top")) + 
      tm_scale_bar( position=c("right", "bottom" ), width=0.2, text.size=0.7) +
      tm_legend( position=c("left", "top") ,  frame=TRUE, scale = 1 , title.size=1.5, text.size=0.80, legend.width=0.75) +
      tm_layout( frame=FALSE )
    }

    print(o)

    if ( outfilename !="" ) {
      if (outformat=="pdf") pdf( file=outfilename, width=9, height=7, bg='white', pointsize=12 )
      if (outformat=="svg") svg( filename=outfilename, width=9, height=7, bg='white', pointsize=12   )
      if (outformat=="png") png( filename=outfilename, width=3072, height=2304, pointsize=12, res=300 )
      print(o)

      dev.off()
      print(outfilename)
    }
    
    return(o)
  }

