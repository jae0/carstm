
  carstm_plotxy = function( X, vn=c( "fit", "summary.random", "time" ), transf=identity, subtype="xy", h=NULL, v=NULL, errorbar_labels=NULL, adj=NULL, offs=0, 
   outfilename=NULL, width_in=9, height_in=7, width_pts=1000, height_pts=800, pres=192,  bg="white", pointsize=12, ... ) {

    dtype = vn[1]
    if ( dtype %in% c("sims", "res", "fit") ) {
      vns = vn[2:length(vn)]
    } else {
      dtype="res"
      vns = vn
    }
    dta = carstm_results_unpack( X, vns ) 

    if (!exists("ID", dta)) {
      dta$ID = row.names(dta)
    }

    if (dtype=="fit") {
      xv = dta$ID
      yv = transf(dta$mean)
      yv_lb = transf(dta[,"0.025quant"])
      yv_ub = transf(dta[,"0.975quant"])
    }

    if (dtype=="res") {
      xv = dta$ID
      yv = transf(dta$mean)
      yv_lb = transf(dta[,"quant0.025"])
      yv_ub = transf(dta[,"quant0.975"])
    } 

    if (dtype=="sims") {
      # vn = "habitat"
      # dta = RES[[mf]][[vn]] # aggregate summaries 
      xv = dta$ID
      yv = transf(dta$mean)
      yv_lb = transf(dta[,"quant0.025"])
      yv_ub = transf(dta[,"quant0.975"])
    }

    if ( !is.null(outfilename) ) {
    
      if (!file.exists( dirname(outfilename) ))  dir.create( dirname(outfilename), recursive=TRUE, showWarnings=FALSE )
      if (grepl("pdf", outfilename)) pdf( file=outfilename, width=width_in, height=height_in, bg=bg, pointsize=pointsize )
      if (grepl("svg", outfilename)) svg( filename=outfilename, width=width_in, height=height_in, bg=bg, pointsize=pointsize   )
      if (grepl("png", outfilename)) png( filename=outfilename, width=width_pts, height=height_pts, pointsize=pointsize, res=pres )
    }

    if (subtype=="xy") {
      plot( yv ~ xv,  ... )
      lines( yv_lb ~ xv, col="gray", lty="dashed")
      lines( yv_ub ~ xv, col="gray", lty="dashed")
    }

    if (subtype=="errorbar" ){
      xv = 1:nrow(dta)
      loffset = xv[1] * 0.175
      plot( yv ~ xv, axes=FALSE, xlim=range(xv)+c(-0.25, 0.25), ...)
      axis(2)
      arrows(x0=xv, y0=yv_lb, x1=xv, y1=yv_ub, code=3, angle=90, length=0.1)
      if (is.null(adj)) {
        text( x=I(xv+loffset), y=yv, labels=errorbar_labels, adj=0.5, offset=4, srt=90 )
      } else {
        text( x=I(xv), y=offs, labels=errorbar_labels, adj=adj, offset=4, srt=90 )

      }
    }

    if (!is.null(h)) {
      for (i in h) abline( h=i,  col="slategray", lty="dashed" )
    }

    if (!is.null(v)) {
      for (i in v) abline( v=i,  col="slategray", lty="dashed" )
    }

    if ( !is.null(outfilename) ) {
      dev.off()
      print(outfilename)
    }

    return(dta)
  }
 