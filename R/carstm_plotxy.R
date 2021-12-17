
  carstm_plotxy = function( X, vn=c( "fit", "summary.random", "time" ), transf=identity, subtype="xy", h=NULL, v=NULL, errorbar_labels=NULL, ... ) {
    dtype = vn[1]
    if ( dtype %in% c("sims", "res", "fit") ) {
      vns = vn[2:length(vn)]
    } else {
      dtype="res"
      vns = vn
    }
    dta = carstm_results_unpack( X, vns ) 
    dev.new()

    if (subtype=="xy") {
      if (dtype=="fit") {
        plot( transf(dta$mean) ~ ID, dta, ... )
        lines( transf(dta[,4]) ~ ID, dta, col="gray", lty="dashed")
        lines( transf(dta[,6]) ~ ID, dta, col="gray", lty="dashed")
      }

      if (dtype=="res") {
        plot( transf(mean) ~ ID , dta, ... )
        lines( transf(quant0.025) ~ ID, dta, col="gray", lty="dashed")
        lines( transf(quant0.975) ~ ID, dta, col="gray", lty="dashed")
      }

      if (dtype=="sims") {
        plot( transf(mean) ~ ID , dta, ... )
        lines( transf(quant0.025) ~ ID, dta, col="gray", lty="dashed")
        lines( transf(quant0.975) ~ ID, dta, col="gray", lty="dashed")

        vn = "habitat"
        ppa = RES[[mf]][[vn]] # aggregate summaries 

        units = attr( ppa, "units")
        plot( ppa[["mean"]] ~ RES$yr, lty=1, lwd=2.5, col="slategray", type="b", main=mf, ylab="Probability", xlab="Year", ylim=c(0.2,0.8), pch=19  )
        lines( ppa[["mean"]] ~ RES$yr, lty=1, lwd=2.5, col="slategray" )
        lines( ppa[["quant0.025"]] ~ RES$yr, lty="dotted", lwd=1, col="slategray"  )
        lines( ppa[["quant0.975"]] ~ RES$yr, lty="dotted", lwd=1, col="slategray"  )
        abline( h=0.5, lty="dashed",  col="slategray" )
      

      }

    }

    if (subtype=="errorbar" ){
      if (dtype=="fit") {
        x = dta$ID
        y = transf(dta$mean)
        y0 = transf(dta[,4])
        y1 = transf(dta[,6])

        plot( y ~ x,  ...)
        arrows(x0=x, y0=y0, x1=x, y1=y1, code=3, angle=90, length=0.1)
        axis(2 )
        par(las=2, mgp=c(3,0,-4), mai=c(2,1,1,1) )
        axis(1, at=x, labels=errorbar_labels, lty="blank", lwd=0.2 )
        abline( h=0.5, lty="dashed",  col="slategray" )
      }
    }

    if (!is.null(h)) {
      for (i in h) abline( h=i,  col="slategray", lty="dashed" )
    }

    if (!is.null(v)) {
      for (i in h) abline( v=i,  col="slategray", lty="dashed" )
    }

    return(dta)
  }
 