  carstm_plot = function( p, res, vn, match_id="StrataID", sppoly=areal_units(p=p), breaksat=NULL, ...) {
    # carstm/aegis wrapper around spplot
    require(sp)
    sppoly@data[,vn] = NA

    time_dimensionality = length( dim(res) ) - 1 # 1 = strata (space)

    if (time_dimensionality==1) { # year
      ii = match( res[,match_id[1]], sppoly@data[,match_id[1]] )
      sppoly@data[, vn] = res[ii, vn]
    }
    if (time_dimensionality==2) { # year/subyear
      ii = match( res[,match_id[1]], sppoly@data[,match_id[1]] )
      jj = match( res[,match_id[2]], sppoly@data[,match_id[2]] )
      sppoly@data[, vn] = res[ii, jj, vn]
    }
    if (time_dimensionality>2) { #year/subyear/?
      ii = match( res[,match_id[1]], sppoly@data[,match_id[1]] )
      jj = match( res[,match_id[2]], sppoly@data[,match_id[2]] )
      kk = match( res[,match_id[3]], sppoly@data[,match_id[3]] )
      sppoly@data[, vn] = res[ii, jj, kk, vn]
      warning ("You need to add more methods here.. ")
    }

    if (length(ii) > 1 ) {
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

