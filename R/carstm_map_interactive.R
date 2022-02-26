 
  carstm_map_interactive = function( 
    res=NULL, 
    toplot=NULL,
    vn=NULL, 
    vn_label=NULL,
    space = "space",
    time= "time",
    cyclic="cyclic",
    stat_var="mean",
    sppoly=NULL,
    smatch=NULL, 
    tmatch=NULL, 
    umatch=NULL, 
    plot_crs=NULL,
    plot_elements=c( "isobaths", "compass", "scale_bar" ),
    additional_features = NULL,
    background = NULL,
    aggregate_function=mean,
    probs=c(0,0.975), 
    outformat="mapview",
    outfilename="",
    map_mode="view",

    legend.title.size=1,
    legend.text.size = 0.6,

    width=9, height=7, bg='white', pointsize=12, pres=300,
    vwidth = 1600,  # in points
    vheight=1200,
    depths = c(50, 100, 200, 400),
    digits=3,
    tmap_zoom=6, 
    id = "space", 
    ...) {


    # simple front end to tmap view mode for plotting

    ellps = list(...)

    require(sf)
    require(tmap)
  
    # tmap_save options:
    width=1000
    height=800
    asp=0

    if (exists("width", ellps) )   width = ellps[["width"]]
    if (exists("height", ellps) )   height = ellps[["height"]]
    if (exists("asp", ellps) )   asp = ellps[["asp"]]

    # tmap plotting options:
    bbox = NULL
    breaks = NULL
    style = "cont"; 
    palette =  "YlOrRd"
    title =  "" 
    showNA =  FALSE 
    lwd =  0.25 
    border.alpha =  0.75
    alpha =   0.95
    legend.is.portrait = TRUE   
    compass_position = c( 0.925, 0.1 )  # (left-right, top-bottom)  
    compass_north = 0
    scale_bar_position = c( 0.75, "BOTTOM" )
    if ( map_mode=="view") scale_bar_position= c("left", "bottom" )
    legend_position = c("left", "top" )
    scale = 2
    legend.title.size=1
    legend.text.size=0.95
    legend.width = ifelse( map_mode=="plot", 0.25, 100 )  # plot_mode uses fractions .. view_mode uses pixels
    scale_bar_width = ifelse( map_mode=="plot", 0.1, 125 )

    if (exists("bbox", ellps) )   bbox = ellps[["bbox"]]
 
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
    if (exists("scale", ellps) )   scale = ellps[["scale"]]
    if (exists("legend.title.size", ellps) )  legend.title.size = ellps[["legend.title.size"]]
    if (exists("legend.text.size", ellps) )    legend.text.size = ellps[["legend.text.size"]]
    if (exists("legend.width", ellps) )   legend.width = ellps[["legend.width"]]


    # if toplot not passed, create from res if given
    if (is.null(toplot)) {
      
      if (!is.null(res)) {
        vv = 0
        toplot = carstm_results_unpack( res, vn ) 
        vv = which(dimnames(toplot)$stat == stat_var) 
        if ( exists("sppoly", res)) {
          if (is.null(sppoly)) sppoly = res[["sppoly"]]
        }
        if (exists(space, res)) {
          suid = res[[space]]
          if (is.null(smatch)) smatch = suid 
          js = match( as.character( sppoly[["AUID"]] ), smatch )  # should match exactly but in case sppoly is a subset 
        }
        if (exists(time, res)) {
          tuid = res[[time]]
          if (is.null(tmatch)) tmatch = tuid
          jt = match( tmatch, res[[time]] )  
        } 
        if (exists(cyclic, res)) {
          uuid = res[[cyclic]]
          if (is.null(umatch)) umatch = uuid
          ju = match( umatch, res[[cyclic]] )  
        }

        data_dimensionality = ifelse (is.vector(toplot), 1, length(dim(toplot) ) )
        if (data_dimensionality==2) {
          toplot = toplot[ js, vv ]  # year only
        } else if (data_dimensionality==3) {
          toplot = toplot[ js, jt, vv ]  # year only
        } else if (data_dimensionality==4) {
          toplot = toplot[ js, jt, ju, vv ] # year/subyear
        }
      } 
    }

    # prepare sppoly
    if (is.null(sppoly)) stop( "sppoly is required")

    if (is.null(plot_crs)) plot_crs = st_crs( sppoly )
    sppoly = st_transform( sppoly, crs=st_crs(plot_crs) )

    if (!exists(space, sppoly)) {
      if (exists("AUID", sppoly)) {
        sppoly[, space] = as.character(sppoly[["AUID"]])  
      } else {
        sppoly[, space] = as.character(1:nrow(sppoly))
      }
    }

    # add toplot to sppoly for final plots, but first check in case toplot is xyz data
    if (!is.null(toplot)) {
      
      ndata = ifelse ( is.vector(toplot), length(toplot), nrow(toplot) )
      if (ndata != nrow(sppoly) ) {
        # must have lon,lat in toplot
        toplot = st_as_sf( toplot, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
        toplot = st_transform(toplot, st_crs( sppoly ) )
        toplot_id = st_points_in_polygons( pts=res, polys=sppoly[, space], varname=space )
        toplot = tapply( toplot[[vn]], toplot_id, aggregate_function, na.rm=TRUE )
      }
 
      if  ( exists("breaks", ellps)) {
        breaks = ellps[["breaks"]] 
        er = range(breaks)
      } else{ 
        er = quantile( toplot, probs=probs, na.rm=TRUE )
        # breaks = signif( seq( er[1], er[2], length.out=7), 2) 
        breaks = pretty( er )
      }
      
      toplot[ which(toplot < er[1]) ] = er[1] # set levels higher than max datarange to max datarange
      toplot[ which(toplot > er[2]) ] = er[2] # set levels higher than max datarange to max datarange
      toplot = round( toplot, digits=digits)

      if (is.vector(toplot) | ( is.array(toplot) && length(dim(toplot) <2 )) ) {
        vn_label =  paste(vn, collapse="_")
        sppoly[, vn_label] = toplot
      } else {
        vn_label = colnames(toplot) = paste(paste(vn, collapse="_"), colnames(toplot), sep="_")
        sppoly = cbind(sppoly, toplot )
      }
      vn_label = gsub(" ", ".", vn_label)
    
    } else {
      # no data sent, assume it is an element of sppoly
      if (!exists(vn, sppoly)) message( paste("variable: ", vn, "not found in sppoly ..."))

      if  ( exists("breaks", ellps)) {
        breaks = ellps[["breaks"]] 
        er = range(breaks)
      } else{ 
        er = range( sppoly[[vn]],   na.rm=TRUE )
        breaks = pretty( er )
      }

      if (is.null(vn_label)) vn_label = vn  # this permits direct plotting of sppoly variables (if toplot and res are not sent)
      sppoly[, vn_label] = round( sppoly[[vn]], digits=digits)

    }

    if (exists("id", ellps) )   id = ellps[["id"]]
    # sppoly = st_cast(sppoly, "POLYGON" )
    sppoly = st_make_valid(sppoly)

    
    tmap_mode( "view" )
    
    tmout = NULL

  
    # https://leaflet-extras.github.io/leaflet-providers/preview/
    # OpenTopoMap, Stamen.Watercolor, Stamen.Terrain, Stamen.TonerLite, Esri.OceanBasemap 
    if (is.null(background)) {
      tmout = tmout +
        tm_basemap(leaflet::providers$CartoDB.PositronNoLabels, alpha=0.8) 
  #     tm_basemap(leaflet::providers$Esri.OceanBasemap, alpha=0.9) +
  #     tm_tiles(leaflet::providers$CartoDB.PositronOnlyLabels, alpha=0.8) 
    } else {
      tmout = tmout + background
    }

    tmout = tmout + 
      tm_shape( sppoly, projection=plot_crs ) +
      tm_polygons( 
        col=vn_label, 
        title= title,
        style = style,
        palette = palette,
        breaks = breaks,
        midpoint = NA ,
        border.col = "lightgray",
        colorNA = NULL,
        id = id,
        showNA =showNA,
        lwd = lwd,  
        border.alpha =border.alpha,
        alpha =alpha, 
        legend.is.portrait = legend.is.portrait ) +

    tmout = tmout + 
      tm_facets(ncol = 2, sync = TRUE) 
      # tm_facets(as.layers = TRUE) 

    if (!is.null(additional_features) ) {
      # e.g. management lines, etc
      tmout = tmout + additional_features 
    }

  
    if ("scale_bar" %in% plot_elements ) {
      tmout = tmout + 
        tm_scale_bar( position=scale_bar_position, width=scale_bar_width, text.size=0.7)  
    }

    tmout = tmout + 
      tm_view(set.view = tmap_zoom, view.legend.position=legend_position  ) +
      tm_layout( legend.text.size=legend.text.size, legend.title.size=legend.title.size, scale=scale, frame=FALSE ) 

    print( tmout ) 
    return(tmout)
  }

