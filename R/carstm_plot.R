
  carstm_plot = function( p, res, vn, poly_match=NULL, time_match=NULL, sppoly=areal_units(p=p), breaksat=NULL, ...) {
    # carstm/aegis wrapper around spplot
    require(sp)
    sppoly@data[,vn] = NA

    # first index is spatial strata
    data_dimensionality = length( dim(res[[vn]]) )

    if (is.null(poly_match)) poly_match = match( res$AUID, sppoly[["AUID"]] )  # should match exactly but in case a subset is sent as sppoly

    if (data_dimensionality==1) {
      sppoly@data[, vn] = res[[vn]] [ poly_match ]  # year only
    }

    if (!is.null(time_match)) {
      n_indexes = length( time_match )
      if (data_dimensionality==2) {
        if (n_indexes==1) sppoly@data[, vn] = res[[vn]] [ poly_match, time_match[[1]] ]  # year only
      }
      if (data_dimensionality==3) {
        if (n_indexes==1) sppoly@data[, vn] = res[[vn]] [ poly_match, time_match[[1]] , ]  # year only
        if (n_indexes==2) sppoly@data[, vn] = res[[vn]] [ poly_match, time_match[[1]], time_match[[2]] ] # year/subyear
      }
    }


    if (length(poly_match) > 1 ) {
      dev.new();
      p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
      if ( exists("coastLayout", p)) {
        sp.layout = p$coastLayout
      } else {
        sp.layout = aegis.coastline::coastline_layout(p=p)
      }
      if ( exists("mypalette", p)) {
        mypalette = p$mypalette
      } else {
        mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")
      }
      if (is.null(breaksat)) breaksat=interval_break(X=sppoly[[vn]], n=length(mypalette), style="quantile")
      spplot( sppoly, vn, main=vn,
        col.regions=mypalette,
        at=breaksat,
        sp.layout=sp.layout$coastLayout,
        col="transparent",
        ...
      )
    }
  }

