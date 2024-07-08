carstm_posterior_predictive_check = function( p, M, logtransf=FALSE )  {

    vn = p$variabletomodel
    iobs = which(M$tag == "observations")
 
    fit = carstm_model( p=p, DS="modelled_fit") #,  sppoly = sppoly )

    pld = data.table(
      observed = M[iobs , ..vn] , 
      fitted = fit$summary.fitted.values[["mean"]] [iobs]
    )
    
    names(pld) = c("observed", "fitted")
     
    if (logtransf) pld$fitted = exp(pld$fitted)

    anno1 = paste( "Pearson correlation: ", round( cor( pld$fitted, pld$observed, use="pairwise.complete.obs" ), 3))
    # cor( fitted, observed, use="pairwise.complete.obs", "spearman" )

    out = ggplot(pld, aes(x =  observed, y = fitted )) +
      geom_abline(slope=1, intercept=0, color="darkgray", lwd=1.4 ) +
      geom_point(color="slategray") +
      labs(caption=anno1, color="slateblue") +
      theme( plot.caption = element_text(hjust = 0, size=12 ) )# move caption to the left 
   
    outputdir = file.path( p$modeldir, p$carstm_model_label )
    fn = file.path(outputdir, "posterior_predictive_check.png" )
    ggsave(filename=fn, plot=out, device="png", width=12, height = 8)
    print(fn)  
    print(out)  
    return(fn)
}
