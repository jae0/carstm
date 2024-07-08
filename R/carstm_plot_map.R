carstm_plot_map = function( p, outputdir, fn_root, additional_features, 
  toplot="random_spatial", 
  probs=c(0.025,0.975), legend.position=c( 0.08, 0.865 ), brks =NULL, palette="-RdYlBu",
  transf=NULL ) {
  
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
 
  if ( grepl("probability", fn_root) | grepl("presence_absence", fn_root)  ) {
    prefix = "habitat"
    tf = identity
  } else if ( grepl("numerical", fn_root) ) {
    prefix = "number"
    tf = log10
  } else if ( grepl("biomass", fn_root) ) {
    prefix = "biomass"
    tf = log10
  } else if ( grepl("weight", fn_root) |  grepl("meansize", fn_root)  ) {
    prefix = "weight"
    tf = log10
  } else {
    prefix = ""
    tf = identity
  }
  if (!is.null(transf)) tf=transf
 
  res = carstm_model( p=p,  DS="carstm_randomeffects" ) 

  if (toplot=="random_spatial") {
    vn=c(  "space", "re_total" ) 
    fn = file.path( outputdir, paste(fn_root, paste0(vn, collapse="_"), "png", sep=".") )
    toplot = carstm_results_unpack( res, vn )
    
    if (is.null(brks)) {
      qn = quantile( tf(toplot[,"mean"]), probs=probs, na.rm=TRUE ) 
      brks = pretty( qn )
    }

    carstm_map(  
        res=res, 
        vn=vn,   
        breaks = brks,
        palette=palette,
        additional_features=additional_features,
        annotation=paste( p$carstm_model_label, "persistent spatial effect" ), 
        legend.position=legend.position,
        outfilename=fn, 
        transformation= tf 
    )
    return("done")
  }


  if (toplot=="predictions") {
    res = carstm_model( p=p,  DS="carstm_predictions" ) 
    vn="predictions"
    toplot = carstm_results_unpack( res, vn )
    qn = quantile( tf(toplot[,,"mean"]), probs=probs, na.rm=TRUE )
    brks = pretty(  qn  )

    if (p$dimensionality == "space") {

      for (y in res$time_name){
        tmatch = as.character(y) 
        u = res$cyclic_id[7]
        fn = file.path( outputdir,  paste(fn_root, paste0(vn, collapse="_"), tmatch, "png", sep=".") )
        carstm_map(  
          res=res, 
          vn=vn, tmatch=tmatch, umatch=as.character(u),
          breaks=brks,
          palette=palette,
          additional_features=additional_features,
          annotation=paste( p$carstm_model_label, "  ", paste0(tmatch, collapse="-") ), 
          legend.position=legend.position,
          outfilename=fn, 
          transformation = tf 
        )
      }

    } else if (p$dimensionality == "space_time") {

      for (y in res$time_name){
        tmatch = as.character(y) 
        u = res$cyclic_id[7]
        fn = file.path( outputdir,  paste(fn_root, paste0(vn, collapse="_"), tmatch, "png", sep=".") )
        carstm_map(  
          res=res, 
          vn=vn, tmatch=tmatch, umatch=as.character(u),
          breaks=brks,
          palette=palette,
          additional_features=additional_features,
          annotation=paste( p$carstm_model_label, "  ", paste0(tmatch, collapse="-") ), 
          legend.position=legend.position,
          outfilename=fn, 
          transformation = tf 
        )
      }
 
    } else if (p$dimensionality == "space_time_season") {
 
    }


    return("done")
  }


    oeffdir = file.path(p$modeldir, p$carstm_model_label, "figures")
    fn_root_prefix = "Predicted_numerical_abundance"
    carstm_plot_marginaleffects( p, oeffdir, fn_root_prefix ) 


# maps
    if (0) {
      ## >>>>>>>>>>>> reconsider file save locations / names:  <<<<<<<<<<<<<<
      outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )
      carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="random_spatial", probs=c(0.025, 0.975) ) 
      carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="predictions", probs=c(0.1, 0.9)) 
    }

}
