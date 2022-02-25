 
  carstm_map = function( 
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
    plot_elements=c( "isobaths", "compass", "scale_bar" ),
    additional_features = NULL,
    aggregate_function=mean,
    probs=c(0,0.975), 
    outfilename=NULL,
    digits = 3,
    ...) {
 
    # thin wrapper around tmap plot mode

    ellps = list(...)

    require(sf)
    require(tmap)

    id = ifelse( exists("id", ellps), ellps[["id"]], "space" )
 
     # tmap_save options:
    outformat = ifelse( exists("outformat", ellps),  ellps[["outformat"]], "tmap" )
    pointsize = ifelse( exists("pointsize", ellps),   ellps[["pointsize"]], 12 )
    width_pts = ifelse( exists("width_pts", ellps),  ellps[["width_pts"]], 1200 )
    height_pts = ifelse( exists("height_pts", ellps),  ellps[["height_pts"]], 800 )
    pres = ifelse( exists("pres", ellps),  ellps[["pres"]], 192 )
    asp = ifelse( exists("asp", ellps),  ellps[["asp"]], 0 )
    width_in = ifelse( exists("width_in", ellps),  ellps[["width_in"]], 9 )
    height_in = ifelse( exists("height_in", ellps),  ellps[["height_in"]], 7 )
    bg = ifelse( exists("bg", ellps),  ellps[["bg"]], 'white' )

    # tmap plotting options:
    
    if ( exists("compass.position", ellps) ) {
      compass.position =  ellps[["compass.position"]] 
    } else {
      compass.position = c("right", "TOP" )
    }

    if ( exists("scale_bar.position", ellps) ) {
      scale_bar.position =  ellps[["scale_bar.position"]] 
    } else {
      scale_bar.position = c("RIGHT", "BOTTOM" )
    }

    if ( exists("legend.position", ellps) ) {
      legend.position =  ellps[["legend.position"]] 
    } else {
      legend.position = c("LEFT", "top" )
    }
 

    scale = ifelse( exists("scale", ellps),  ellps[["scale"]],  2 )

    style = ifelse( exists("style", ellps),   ellps[["style"]],  "cont" )
    palette = ifelse( exists("palette", ellps),   ellps[["palette"]],  "YlOrRd" )
    title = ifelse( exists("title", ellps),   ellps[["title"]],  "" )
    showNA = ifelse( exists("showNA", ellps),   ellps[["showNA"]],  FALSE )
    lwd = ifelse( exists("lwd", ellps),  ellps[["lwd"]],  0.04  )
    border.alpha = ifelse( exists("border.alpha", ellps),  ellps[["border.alpha"]],  0.7 )
    alpha = ifelse( exists("alpha", ellps),   ellps[["alpha"]],  0.9 )


    compass.north = ifelse( exists("compass.north", ellps),  ellps[["compass.north"]],  0 )

    scale_bar.width = ifelse( exists("scale_bar.width", ellps),  ellps[["scale_bar.width"]], 0.1 )

    legend.title.size = ifelse( exists("legend.title.size", ellps),  ellps[["legend.title.size"]],  1 )
    legend.text.size = ifelse( exists("legend.text.size", ellps),   ellps[["legend.text.size"]],  1 )
    legend.width = ifelse( exists("legend.width", ellps),  ellps[["legend.width"]], 0.75 )

    legend.is.portrait = ifelse( exists("legend.is.portrait", ellps),   ellps[["legend.is.portrait"]],  TRUE )
 

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

    # cannot use ifelse as it is not a singleton
    if( exists("plot_crs", ellps)) {
      plot_crs = ellps[["plot_crs"]] 
    } else {
      plot_crs = st_crs( sppoly )
    } 
 
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

    sppoly = st_make_valid(sppoly)

    tmout = NULL  
 
    tmout =  tmout + 
      tm_shape( sppoly, projection = plot_crs ) +
      tm_polygons( 
        col = vn_label, 
        title= title,
        style = style,
        palette = palette,
        breaks = breaks,
        midpoint = NA ,
        border.col = "gray30",
        colorNA = NULL,
        id = id,
        showNA =showNA,
        lwd = lwd,  
        border.alpha =border.alpha,
        alpha =alpha, 
        legend.is.portrait = legend.is.portrait ) 
 
    if (!is.null(additional_features) ) {
      # e.g. management lines, etc
      tmout = tmout + additional_features 
    }

    if ("compass" %in% plot_elements ) {
      tmout = tmout + 
        tm_compass( position=compass.position, size=1) 
      message( "You will need to tweak more settings such as angle of compass .." ) 
    }

    if ("scale_bar" %in% plot_elements ) {
      tmout = tmout + 
        tm_scale_bar( position=scale_bar.position, width=scale_bar.width, text.size=0.5)  
    }

    tmout = tmout + 
      tm_layout( frame=FALSE, legend.position=legend.position, scale=scale, legend.title.size=legend.title.size,
        legend.text.size =legend.text.size, legend.width=legend.width ) 


    if ( !is.null(outfilename) ) {

      if (outformat=="tmap") {
        tmap_save( tmout, outfilename, width=width_pts, height=height_pts, asp=asp, dpi=pres, scale=scale*0.4)
      }

      if (outformat %in% c("pdf", "svg", "png")){
        if (outformat=="pdf") pdf( file=outfilename, width=width_in, height=height_in, bg=bg, pointsize=pointsize )
        if (outformat=="svg") svg( filename=outfilename, width=width_in, height=height_in, bg=bg, pointsize=pointsize   )
        if (outformat=="png") png( filename=outfilename, width=width_pts, height=height_pts, pointsize=pointsize, res=pres )
          print(tmout)
        dev.off()
      }

      print(outfilename)
    }
 
    return(tmout)
  }

