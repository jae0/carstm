
  carstm_map = function( res, vn, poly_match=NULL, time_match=NULL, ...) {

    # carstm/aegis wrapper around spplot, will forward any additional spplot args 
    # TODO move to an sf-based plot

    require(sp)
    sppoly = as( res$sppoly, "Spatial")
    slot(sppoly, "data")[,vn] = NA

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

    slot(sppoly, "data")[, vn] = toplot

    maintitle = list(label="my plot title", cex=2)

    if (length(poly_match) > 1 ) {
      dev.new();
      ellp = list(...)
      if ( !exists("main", ellp ) )  ellp[["main"]]=maintitle
      ellp$obj = sppoly
      ellp$zcol=vn
      ellp$col="transparent"
      do.call(spplot, ellp )
    }
  }

