
carstm_prior_posterior_compare = function( hypers, all.hypers, i=1, vn=NULL, xrange=NULL, transf=TRUE ) {
    
    # extracted from INLA:::plot.inla()
    # all.hypers = INLA:::inla.all.hyper.postprocess(fit$all.hyper)
    # hypers = fit$marginals.hyperpar
    # names(hypers)
    if (is.null(vn)) {
        vn = INLA:::inla.nameunfix(names(hypers)[i])
    }  
    label = vn
    hh = hypers[[vn]]

    id = unlist(strsplit(attr(hh, "hyperid"), "\\|"))
    section = tolower(id[2])
    hyperid = id[1]
    iposterior = inla.smarginal(hh)

    if (is.null(xrange)) {
        xrange = range(iposterior$x) 
    }
    
    iprior = INLA:::inla.get.prior.xy( 
        section = section, hyperid=hyperid, all.hyper= all.hypers,  
        range = range(xrange), intern = FALSE, debug = FALSE
    )

    if (transf)   {
        if (grepl("precision", vn, ignore.case =TRUE )) {
            iposterior = inla.tmarginal( fun=function(y) 1/sqrt( y ), iposterior )
            iprior = inla.tmarginal( fun=function(y) 1/sqrt( y ), iprior )
            label = gsub( "Precision", "SD", label)
        }
    }

    ipo = as.data.frame(iposterior)
    ipo$tag = "posterior"

    ipr = as.data.frame(iprior)
    ipr$tag ="prior"

    o = rbind( ipr, ipo )
    
    plt = ggplot(o, aes(x=x, y=y, color=tag)) + 
        geom_line( alpha=0.9, linewidth=1.2  ) + 
        ylab('Density') + xlab(label) +
        theme_light( base_size = 22 ) +
        theme( legend.position="inside", legend.position.inside=c(0.8, 0.9), 
            legend.title=element_blank()) 
    return(plt)
}
