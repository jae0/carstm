
  carstm_map = function( 
    res=NULL, 
    vn=NULL, 
    space = "space",
    time= "time",
    season="season",
    stat_var="mean",
    sppoly=NULL,
    smatch=NULL, 
    tmatch=NULL, 
    plot_crs=NULL,
    plot_elements=c( "isobaths", "compass", "scale_bar", "legend" ),
    additional_polygons = NULL,
    aggregate_function=mean,
    probs=c(0,0.975), 
    outformat="mapview",
    outfilename="",
    map_mode="view",
    width=9, height=7, bg='white', pointsize=12, pres=300,
    depths = c(50, 100, 200, 400),
    digits=3,
    tmap_zoom=7, 
    ...) {


    if (0) {
      space = "space"
      time= "time"
      season="season"
      stat_var="mean"
      sppoly=NULL
      smatch=NULL 
      tmatch=NULL 
      plot_crs=NULL
      coastline=NULL 
      isobaths=NULL 
      managementlines = NULL
      plot_elements=NULL
      aggregate_function=mean
      probs=c(0,0.975) 
      outformat="mapview"
      map_mode="view"
      width=9; height=7; bg='white'; pointsize=12; res=300
      depths = c(50, 100, 200, 400)
      digits=3

      p = aegis.bathymetry::bathymetry_parameters( project_class="carstm" )  # defaults are hard coded
      res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
      require(fields)
      
      coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=p$aegis_proj4string_planar_km )
      isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=p$aegis_proj4string_planar_km )
    
      vn = "predictions" 
      vn = c("random", "space", "combined") 
      

      carstm_map(  
          res = res,    
          vn = vn,
          breaks = pretty(p$discretization$z),
          palette = "viridis",
          main = "Bathymetry predicted",
          coastline = coastline,
          isobaths = isobaths 
      )

 
 
    }


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
    id = row.names(sppoly)
    breaks = NULL
    style = "cont"; 
    palette =  "YlOrRd"
    title =  "" 
    showNA =  FALSE 
    lwd =  0.25 
    border.alpha =  0.5
    alpha =   0.9
    legend.is.portrait = FALSE   
    compass_position = c( "right", "bottom" )   
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

    
    vv = 0
    if ( exists("sppoly", res)) {
      if (is.null(sppoly)) sppoly = res[["sppoly"]]
      toplot = carstm_results_unpack( res, vn ) 
      vv = which(dimnames(toplot)$stat == stat_var) 
      vn = paste(vn, collapse="_")
    }  else {
      # if xyz is given then then are point data to be aggregated into polygons
      if (!exists(space, sppoly)) {
        if (exists("AUID", sppoly)) {
          sppoly[, space] = as.character(sppoly[, "AUID"])  
        } else {
          sppoly[, space] = as.character(1:nrow(sppoly))
        }
      }
      res = st_as_sf( res, coords= c("lon", "lat"), crs=st_crs( projection_proj4string("lonlat_wgs84") ) )
      res = st_transform(res, st_crs( sppoly ) )
      res[[space]] = st_points_in_polygons( pts=res, polys=sppoly[, space], varname=space )
      toplot = tapply( res[[vn]], res[[space]], aggregate_function, na.rm=TRUE )
    }

    if (exists(space, res)) {
      if (is.null(smatch)) smatch = match( as.character(sppoly[["AUID"]]), res[[space]] )  # should match exactly but in case a subset 
    }
    if (exists(time, res)) {
      tuid = res[[time]]
      if (is.null(tmatch)) tmatch = tuid
    } 
    if (exists(season, res)) {
      uuid = res[[season]]
      if (is.null(umatch)) umatch = uuid
    }

  
    if (is.null(plot_crs)) plot_crs = st_crs( sppoly )

    sppoly = st_transform( sppoly, crs=st_crs(plot_crs) )

    data_dimensionality = ifelse (is.vector(res), 1, length(dim(res) ) )

    if (data_dimensionality==1) {
      toplot = toplot[ smatch ]
    } else if (data_dimensionality==2) {
      toplot = toplot[ smatch, vv ]  # year only
    } else if (data_dimensionality==3) {
      toplot = toplot[ smatch, tmatch, vv ]  # year only
    } else if (data_dimensionality==4) {
      toplot = toplot[ smatch, tmatch, umatch, vv ] # year/subyear
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
      plot_names = vn
      sppoly[, plot_names] = toplot
    } else {
      plot_names = colnames(toplot) = paste(vn, colnames(toplot), sep="_")
      sppoly = cbind(sppoly, toplot )
    }
    plot_names = gsub(" ", ".", plot_names)


    tmap_mode( map_mode )

    # https://leaflet-extras.github.io/leaflet-providers/preview/
    # OpenTopoMap, Stamen.Watercolor, Stamen.Terrain, Stamen.TonerLite, Esri.OceanBasemap 
    tmout = 
      tm_basemap(leaflet::providers$CartoDB.Positron, alpha=0.8) +
#      tm_basemap(leaflet::providers$Esri.OceanBasemap, alpha=0.9) +
      tm_tiles(leaflet::providers$CartoDB.PositronOnlyLabels, alpha=0.8) +
      tm_shape( sppoly, projection=plot_crs ) +
      tm_polygons( 
        col=plot_names, 
        title= ifelse ( exists("main", ellps), ellps[["main"]], vn ) ,
        style = ifelse ( exists("style", ellps), ellps[["style"]], "cont" ) ,
        palette = ifelse ( exists("palette", ellps), ellps[["palette"]], "YlOrRd"),
        breaks = breaks,
        border.col = "lightgray",
        colorNA = NULL,
        id = id,
        showNA =showNA,
        lwd = lwd,  
        border.alpha =border.alpha,
        alpha =alpha, 
        legend.is.portrait = legend.is.portrait ) +
      tm_facets(ncol = 2, sync = TRUE) 
      # tm_facets(as.layers = TRUE) 

    if (!is.null(additional_polygons) ) {
      # e.g. management lines, etc
      for (poly in additional_polygons) {
        if (any( st_is(poly, c("MULTILINESTRING", "LINESTRING") ) ) ) {
          tmout = tmout + 
            tm_shape( poly, projection=plot_crs ) +
            tm_lines( col="grey40", alpha=0.5, lwd=2)  
        }
        if (all( (st_is(poly, c("MULTIPOLYGON", "POLYGON"))) )) {
          tmout = tmout + 
            tm_shape( poly, projection=plot_crs ) +
            tm_polygons( col="grey80", alpha=0.5 )    
        }
      }
    }

    if ("rnaturalearth_countries" %in% plot_elements ) {
        if (!require("rnaturalearth")) {
          install.packages("rnaturalearth")
          if (0) {
            # or install directly
            devtools::install_github("ropensci/rnaturalearth")
            devtools::install_github("ropensci/rnaturalearthdata")
            install.packages("rnaturalearthhires",  repos = "http://packages.ropensci.org", type = "source")
          }
        }              
      require("rnaturalearth")
      world = ne_countries(type = 'countries', scale = 'large', returnclass="sf")
      tmout = tmout + 
        tm_shape( plot_elements, projection=plot_crs ) + 
        tm_borders( col="gray", alpha=0.75)  
    }

    if ("coastline" %in% plot_elements ) {
      coastline = aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs ) 
      tmout = tmout + 
        tm_shape( coastline, projection=plot_crs ) +
        tm_polygons( col="grey80", alpha=0.1 ) 
    }

    if ("isobaths" %in% plot_elements ) {
      isobaths = aegis.bathymetry::isobath_db( depths=depths, project_to=plot_crs  )
      tmout = tmout + 
        tm_shape( isobaths, projection=plot_crs ) +
        tm_lines( col="lightgray", alpha=0.75) 
    }


    if ("compass" %in% plot_elements ) {
      if (map_mode=="plot"){
        tmout = tmout + 
          tm_compass( position=compass_position) 
      }
    }

    if ("scale_bar" %in% plot_elements ) {
      tmout = tmout + 
        tm_scale_bar( position=c("right", "bottom" ), width=0.2, text.size=0.7) 
      # tm_scale_bar( position=scale_bar_position, width=0.15, text.size=0.7) +  
    }

    if ("legend" %in% plot_elements ) {
      tmout = tmout + 
        tm_legend( position=c("left", "top") ,  frame=TRUE, scale = 1 , title.size=1.5, text.size=0.80, legend.width=0.75) 
      # tm_legend( position=legend_position, frame=FALSE, scale = legend_scale , title.size=legend_title.size, text.size=legend_text.size, legend.width=legend.width) +  
    }

    tmout = tmout + 
        tm_layout( frame=FALSE ) +
        tm_view(set.view = c(tmap_zoom))


        # tm_layout( frame=FALSE, title=title ) 
    

    if ( outfilename !="" ) {

      if (outformat=="tmap") {
        # tmap_save options:
        twidth=1000
        theight=800
        tasp=0
        if (exists("twidth", ellps) )   twidth = ellps[["twidth"]]
        if (exists("theight", ellps) )   theight = ellps[["theight"]]
        if (exists("tasp", ellps) )   tasp = ellps[["tasp"]]
        tmap_mode("plot")
        tmap_save( tmout, outfilename, width=twidth, height=theight, asp=tasp )
        print(outfilename)
        return(tmout)
      }

      if (outformat=="mapview") {
        if (!require(mapview)) install.packages("mapview")
        if (!require(webshot)) {
          install.packages("webshot")
          webshot::install_phantomjs()
        }
        require(mapview)
        mapshot( tmap_leaflet(tmout), file=outfilename )
        print(outfilename)
        return(tmout)
      }

      if (outformat=="pdf") pdf( file=outfilename, width=width, height=height, bg=bg, pointsize=pointsize )
      if (outformat=="svg") svg( filename=outfilename, width=width, height=height, bg=bg, pointsize=pointsize   )
      if (outformat=="png") png( filename=outfilename, width=3072, height=2304, pointsize=pointsize, res=pres )
        print(tmout)
      dev.off()
      print(outfilename)

    }
    
    return(tmout)
  }

