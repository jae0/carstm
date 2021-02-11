
  carstm_map = function( res, vn, poly_match=NULL, time_match=NULL, coastline=NULL, isobaths=NULL, plot_crs=NULL, ...) {

    # carstm/aegis wrapper around spplot, will forward any additional spplot args 
    # TODO move to an sf-based plot

    require(sp)
    sppoly = res$sppoly

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
    if (!is.null(plot_crs)) sppoly = st_transform( sppoly, crs=st_crs(plot_crs) )


    ellp = list(...)

    if (length(poly_match) > 1 ) {
      plot( sppoly[vn], reset = FALSE, border = "lightslateblue", lwd = 0.8, ...  )
      if (!is.null(coastline) ) {
        coastline = st_transform( coastline, crs=st_crs(sppoly) )
        do.call(plot, list(x=coastline, col="whitesmoke", lwd=1.0, add=TRUE) )
      }
      
      if (!is.null(isobaths) ) {
        isobaths = st_transform( isobaths, crs=st_crs(sppoly) )
        do.call(plot, list(x=isobaths, col="lightgray", lwd=0.25, add=TRUE) )
      }
    }
  }

