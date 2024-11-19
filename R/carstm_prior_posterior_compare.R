
carstm_prior_posterior_compare = function( fit=NULL, vn=NULL, xrange=NULL, transf=TRUE, 
    outputdir = tempdir() ) {

    # extracted from INLA:::plot.inla()
    if (0) {
        # usage: 
        fit = carstm_model(  p=p, DS="modelled_fit" )      
		    ( vns = names(fit$marginals.hyperpar) )
        for (i in 1:length( fit$marginals.hyperpar )) {
            o = plot_prior_posterior( fit, vn=vns[i])
            dev.new() 
            print(o)
        }

        # or: 

        res = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and summary

        ( vns = c( names( res$hypers), names(res$fixed) ) )
  
        for (i in 1:length(vns) ) {
            o = carstm_prior_posterior_compare( res, vn=vns[i] )  
            dev.new(); print(o)
        }     

    }

    if ( exists("marginals.hyperpar", fit) ) {  
        # fit object from INLA    
        hypers = fit$marginals.hyperpar
        fixed = fit$marginals.fixed
    } else if ( exists("hypers", fit) ) {  
        # carstm results summary
        hypers = fit$hypers
        fixed = fit$fixed
    }

    all.hypers = INLA:::inla.all.hyper.postprocess(fit$all.hyper)
 
    iposterior = iprior = NULL

    if (!is.null(fixed)) {
        if (exists(vn, fixed)) {
            iposterior = inla.smarginal(fixed[[vn]])
            if (is.null(xrange)) xrange = range(iposterior$x) 
            iprior = INLA:::inla.get.prior.xy(
                section = "fixed", 
                hyperid = vn, 
                all.hyper = all.hypers, 
                range = xrange
            )
        }
    }

    if (!is.null(hypers)) { 
        if ( exists(vn, hypers )) {
            hyperid = attr(hypers[[vn]], "hyperid")
            if (grepl("|", hyperid)) hyperid = strsplit(hyperid, "\\|")
            id = unlist(hyperid)
            section = tolower(id[2])
            hyperid = id[1]

            iposterior = inla.smarginal(hypers[[vn]])
            
            if (is.null(xrange)) xrange = range(iposterior$x) 
            
            iprior = INLA:::inla.get.prior.xy( 
                section = section, hyperid=hyperid, all.hyper= all.hypers,  
                range = xrange, intern = FALSE, debug = FALSE
            )
            
            if (transf)   {
                if (grepl("precision", vn, ignore.case =TRUE )) {
                    iposterior = inla.tmarginal( fun=function(y) 1/sqrt( y ), iposterior, method = "linear" )
                    iprior = inla.tmarginal( fun=function(y) 1/sqrt( y ), iprior, method = "linear" )
                    vn = gsub( "Precision", "SD", vn)
                }
            }
            
        } 
    }

    ipo_xysum = 1
    ipo = NULL

    if (!is.null(iposterior)) {
        ipo = as.data.frame(iposterior)
        ipo$tag = "posterior"

        dx = diff(ipo$x)
        ipo$xy = c(median(dx), dx ) * ipo$y

        ipo_xysum = sum(ipo$xy, na.rm=TRUE)
    }
 
    ipr = as.data.frame(iprior)
    ipr$tag ="prior"

    dx = diff(ipr$x)
    ipr$xy = c(median(dx), dx ) * ipr$y
    ipr_xysum = sum(ipr$xy, na.rm=TRUE)
    ipr$y = ipr$y * ( ipo_xysum / ipr_xysum)

    o = rbind( ipr, ipo )

    plt = ggplot(o, aes(x=x, y=y, color=tag)) + 
    geom_line( alpha=0.9, linewidth=1.2  ) + 
    ylab('Density') + xlab(vn) +
    theme_light( base_size = 22 ) +
    theme( legend.position="inside", legend.position.inside=c(0.8, 0.9), 
        legend.title=element_blank()) 
    
    fn = file.path(outputdir, paste("prior_posterior_compare_", vn, ".png", sep="" ) )
    ggsave(filename=fn, plot=plt, width=12, height = 8)
    print(fn)

    return(plt)
  
}
