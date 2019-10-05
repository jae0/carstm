  carstm_plot = function( p, res, vn, match_id="StrataID", sppoly=areal_units(p=p), intervals=NULL, ...) {
    # carstm/aegis wrapper around spplot
    require(sp)
    sppoly@data[,vn] = NA

    dimensionality = length( dim(res) )

    warning("needs to be tested")

    if (dimensionality==1) {
      ii = match( res[,match_id[1]], sppoly@data[,match_id[1]] )
      toplot = res$vn[ii,]
    }
    if (dimensionality==2) {
      ii = match( res[,match_id[1]], sppoly@data[,match_id[1]] )
      jj = match( res[,match_id[2]], sppoly@data[,match_id[2]] )
      toplot = res$vn[ii, jj]
    }
    if (dimensionality==3) {
      ii = match( res[,match_id[1]], sppoly@data[,match_id[1]] )
      jj = match( res[,match_id[2]], sppoly@data[,match_id[2]] )
      kk = match( res[,match_id[3]], sppoly@data[,match_id[3]] )
      toplot = res[ii, jj, kk]
    }
    if (dimensionality>3) {
      stop ("You need to add more methods here.. ")
    }

    if (length(ii) > 1 ) {
      sppoly@data[,vn] = toplot
      dev.new();
      p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
      sp.layout = ifelse( exists("coastLayout", p), p$coastLayout, aegis.coastline::coastline_layout(p=p) )
      mypalette = ifelse( exists("mypalette", p), p$mypalette, RColorBrewer::brewer.pal(9, "YlOrRd") )
      if (is.null(intervals)) intervals=interval_break(X=sppoly[[vn]], n=length(mypalette), style="quantile")
      spplot( sppoly, vn, main=vn,
        col.regions=mypalette,
        at=intervals,
        sp.layout=sp.layout,
        col="transparent",
        ...
      )
    }
  }

