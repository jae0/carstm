
  carstm_map = function(
    res=NULL,
    toplot=NULL,
    vn=NULL,
    vn_label=NULL,
    space = "space_name",
    time= "time_name",
    cyclic="cyclic_name",
    stat_var="mean",
    alpha=0.95,
    modelinfo=NULL,
    sppoly=NULL,
    smatch=NULL, 
    tmatch=NULL, 
    umatch=NULL, 
    plot_elements=c( "isobaths", "compass", "scale_bar" ),
    aggregate_function=mean,
    probs=c(0,0.975),
    outfilename=NULL,
    outscale= 1,
    digits = 3,
    transformation=NA,
    plotmethod="ggplot",
    tmapmode="view",
    colors=RColorBrewer::brewer.pal(5, "YlOrRd"), # display.brewer.all()
    ...) {
 

    ellps = list(...)


    require(sf)
    require(RColorBrewer) 

    # if toplot not passed, create from res if given
    if (is.null(toplot)) {

      if (!is.null(res)) {
        vv = 0
        toplot = carstm_results_unpack( res, vn )

        vv = which(dimnames(toplot)$stat == stat_var)
        if ( exists("sppoly", modelinfo)) {
          if (is.null(sppoly)) sppoly = modelinfo[["sppoly"]]
        }
        if (exists(space, modelinfo)) {
          suid = modelinfo[[space]]
          if (is.null(smatch)) smatch = suid
          js = match( as.character( sppoly[["AUID"]] ), smatch )  # should match exactly but in case sppoly is a subset
        }
        if (exists(time, modelinfo)) {
          tuid = modelinfo[[time]]
          if (is.null(tmatch)) tmatch = tuid
          jt = match( tmatch, modelinfo[[time]] )
        }
        if (exists(cyclic, modelinfo)) {
          uuid = modelinfo[[cyclic]]
          if (is.null(umatch)) umatch = uuid
          ju = match( umatch, modelinfo[[cyclic]] )
        }

        data_dimensionality = ifelse (is.vector(toplot), 1, length(dim(toplot) ) )
        if (data_dimensionality==2) {
          toplot = toplot[ js, vv ]  # year only
        } else if (data_dimensionality==3) {
          toplot = toplot[ js, jt, vv ]  # year only
        } else if (data_dimensionality==4) {
          toplot = toplot[ js, jt, ju, vv ] # year/subyear
        }
        if (!is.na(transformation)) toplot = transformation(toplot)
        if (exists("filter", sppoly)) toplot = toplot * sppoly[["filter"]]

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
        if (!is.na(transformation)) toplot = transformation(toplot)
        if (exists("filter", sppoly)) toplot  = toplot  * sppoly[["filter"]]
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
      if (length(vn) ==1 ) {
        if (!exists(vn, sppoly)) message( paste("variable: ", vn, "not found in sppoly ..."))

        if (!is.na(transformation)) sppoly[[vn]] = transformation(sppoly[[vn]])
        if (exists("filter", sppoly)) sppoly[[vn]] = sppoly[[vn]] * sppoly[["filter"]]

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
    }

    sppoly = st_make_valid(sppoly)
    attributes(sppoly[[vn_label]]) = NULL

 
 
    if (plotmethod=="ggplot") {
 # copy of simple .. elaborate as required ..
      require(ggplot2)

      title = ifelse( exists("title", ellps), ellps[["title"]],  "" )
      
      if ( exists("legend.position", ellps)) {
        legend.position = ellps[["legend.position"]] 
      } else {
        legend.position = c( 0.925, 0.15 ) 
      }

      pltadd = NULL 
      if ( exists("additional_features", ellps) ) {
        # e.g. management lines, etc
        pltadd = ellps[["additional_features"]]
        if (exists("ggplot", pltadd ) ) pltadd = pltadd[["ggplot"]][["layers"]]
      } 
      
      annot = NULL
      if ( exists("annotation", ellps) ) {
        annot = labs(caption = ellps[["annotation"]] )
      }

      bb = st_bbox(sppoly)
      xr = c(bb["xmin"], bb["xmax"])
      yr = c(bb["ymin"], bb["ymax"])
#, colour="gray90"
 
      plt = ggplot() +
        geom_sf( data=sppoly, aes(fill=.data[[vn_label]], alpha=0.95), lwd=0 )  +
        scale_fill_gradientn(name = vn_label, 
          limits=range(breaks),
          colors=alpha(colors, alpha=0.99), na.value=NA ) +
        guides(fill = guide_colorbar(
          # title.theme = element_blank(), 
          # title.theme = element_text(size = 20),
          label.theme = element_text(size = 16) ) ) +
        scale_alpha(range = c(0.8, 0.95), guide = "none") +
        annot +
        pltadd  +
        coord_sf(xlim =xr, ylim =yr, expand = FALSE) +
        theme(
          axis.line=element_blank(),
          # axis.text.x=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          legend.position=legend.position,
          legend.title = element_blank(),
          # panel.background=element_blank(),
          panel.background = element_rect(fill =NA),
          panel.border=element_blank(),
          # panel.grid.major=element_blank(),
          panel.grid.major = element_line(color = "grey"),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(), 
          plot.caption = element_text(hjust = 0, size = 14)
      )
        

      if ( !is.null(outfilename) ) {
          if (!file.exists( dirname(outfilename) ))  dir.create( dirname(outfilename), recursive=TRUE, showWarnings=FALSE )
          print(plt)
          ggsave( outfilename )
          # ggsave( outfilename, width=width_pts, height=height_pts, asp=asp, dpi=pres, scale=scale*outscale)
          print(outfilename)
      }
      return(plt)
    }

    

    if (plotmethod=="tmap") {

      # thin wrapper around tmap plot mode

      require(tmap)

      id = ifelse( exists("id", ellps), ellps[["id"]], "space" )
  
      # tmap_save options:
      outscale = ifelse( exists("outscale", ellps),  ellps[["outscale"]],  0.4 )  # to match screen to file outputs (not same)
      outformat = ifelse( exists("outformat", ellps),  ellps[["outformat"]], "tmap" )
      pointsize = ifelse( exists("pointsize", ellps),   ellps[["pointsize"]], 12 )
      width_pts = ifelse( exists("width_pts", ellps),  ellps[["width_pts"]], 1600 )
      height_pts = ifelse( exists("height_pts", ellps),  ellps[["height_pts"]], 1200 )
      pres = ifelse( exists("pres", ellps),  ellps[["pres"]], 600 )
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
      lwd = ifelse( exists("lwd", ellps),  ellps[["lwd"]],  0.02  )
      border.alpha = ifelse( exists("border.alpha", ellps),  ellps[["border.alpha"]],  0.75 )
      alpha = ifelse( exists("alpha", ellps),   ellps[["alpha"]],  0.975 )


      compass.north = ifelse( exists("compass.north", ellps),  ellps[["compass.north"]],  0 )

      scale_bar.width = ifelse( exists("scale_bar.width", ellps),  ellps[["scale_bar.width"]], 0.1 )

      legend.title.size = ifelse( exists("legend.title.size", ellps),  ellps[["legend.title.size"]],  1 )
      legend.text.size = ifelse( exists("legend.text.size", ellps),   ellps[["legend.text.size"]],  0.75 )
      legend.width = ifelse( exists("legend.width", ellps),  ellps[["legend.width"]], 0.75 )

      legend.is.portrait = ifelse( exists("legend.is.portrait", ellps),   ellps[["legend.is.portrait"]],  TRUE )
  


      tmap_mode( tmapmode )

      plt = NULL

      plt =  plt + 
        tm_shape( sppoly, projection = plot_crs ) +
        tm_polygons( 
          col = vn_label, 
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
          legend.is.portrait = legend.is.portrait )
  
      if ( exists("additional_features", ellps) ) {
        # e.g. management lines, etc
        pltadd = ellps[["additional_features"]] 
        if (exists("tmap", pltadd ) ) pltadd = pltadd[["tmap"]]
        plt = plt + pltadd 
      }


      if ("compass" %in% plot_elements ) {
        plt = plt + 
          tm_compass( position=compass.position, size=1.5) 
        message( "You might need to tweak  angle of compass .. send 'north=xx" ) 
      }

      if ("scale_bar" %in% plot_elements ) {
        plt = plt +
          tm_scale_bar( position=scale_bar.position, width=scale_bar.width, text.size=0.5)
      }

      plt = plt +
        tm_layout( frame=FALSE, legend.position=legend.position, scale=scale, legend.title.size=legend.title.size,
          legend.text.size =legend.text.size, legend.width=legend.width ) 

      tmap_mode(tmapmode)

      if ( !is.null(outfilename) ) {
        
        if (!file.exists( dirname(outfilename) ))  dir.create( dirname(outfilename), recursive=TRUE, showWarnings=FALSE )

        if (outformat=="tmap") {
          tmap_save( plt, outfilename, width=width_pts, height=height_pts, asp=asp, dpi=pres, scale=scale*outscale)
        }
 
        if (outformat=="ggsave") {
          print(plt)
          ggsave( outfilename )
          # ggsave( outfilename, width=width_pts, height=height_pts, asp=asp, dpi=pres, scale=scale*outscale)
        }

        if (outformat %in% c("pdf", "svg", "png")){
          if (outformat=="pdf") pdf( file=outfilename, width=width_in, height=height_in, bg=bg, pointsize=pointsize )
          if (outformat=="svg") svg( filename=outfilename, width=width_in, height=height_in, bg=bg, pointsize=pointsize   )
          if (outformat=="png") png( filename=outfilename, width=width_pts, height=height_pts, pointsize=pointsize, res=pres )
            print(plt)
          dev.off()
        }

        print(outfilename)
      }
  
      return(plt)
    }

  }

