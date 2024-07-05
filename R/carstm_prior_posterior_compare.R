
carstm_prior_posterior_compare = function( hypers, all.hypers, i=1 ) {
    
    # extracted from INLA:::plot.inla()
    # all.hypers = INLA:::inla.all.hyper.postprocess(fit$all.hyper)
    # hypers = fit$marginals.hyperpar
    # names(hypers)

    label = INLA:::inla.nameunfix(names(hypers)[i])
    hh = hypers[[i]]
    id = unlist(strsplit(attr(hh, "hyperid"), "\\|"))
    section = tolower(id[2])
    hyperid = id[1]
    iposterior = inla.smarginal(hh)
    
    xrange = range(iposterior$x)
    xrange = xrange + xrange * 0.1 * c(-1, 1)
    iprior = INLA:::inla.get.prior.xy( 
        section = section, hyperid=hyperid, all.hyper= all.hypers,  
        range = range(iposterior$x), intern = FALSE, debug = FALSE
    )
 
    ipo = as.data.frame(iposterior)
    ipo$tag = "posterior"

    ipr = as.data.frame(iprior)
    ipr$tag ="prior"

    o = rbind( ipr, ipo )
    if (grepl("precision", label,ignore.case =TRUE )) o$x = 1/sqrt(o$x)
    label = gsub( "Precision", "SD", label)
    
    plt = ggplot(o, aes(x=x, y=y, color=tag)) + 
        geom_line( alpha=0.9, linewidth=1.2  ) + 
        ylab('Density') + xlab(label) +
        theme_light( base_size = 22 ) +
        theme( legend.position="inside", legend.position.inside=c(0.8, 0.9), 
            legend.title=element_blank()) 
    return(plt)
}
