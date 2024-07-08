carstm_plot_map = function( p=NULL, 
  toplot="random_spatial", 
  outputdir=tempdir(), fn=NULL, fn_root_prefix=NULL, 
  additional_features=NULL, 
  annotation = NULL,
  probs=c(0.025,0.975), legend.position=c( 0.08, 0.865 ), brks =NULL, 
  colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")), 
  # colors=RColorBrewer::brewer.pal(5, "YlOrRd"),
  transf=NULL ) {

  require(ggplot2)

  if (is.null(p)) stop("Require parameter list 'p'")

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
 
  prefix = ""
  tf = identity
  if (!is.null(fn_root_prefix)) {

    if ( grepl("probability", fn_root_prefix) | grepl("presence_absence", fn_root_prefix)  ) {
      prefix = "habitat"
      tf = identity
    } else if ( grepl("numerical", fn_root_prefix) ) {
      prefix = "number"
      tf = log10
    } else if ( grepl("biomass", fn_root_prefix) ) {
      prefix = "biomass"
      tf = log10
    } else if ( grepl("weight", fn_root_prefix) |  grepl("meansize", fn_root_prefix)  ) {
      prefix = "weight"
      tf = log10
    } 
  }

  if (!is.null(transf)) tf=transf
 
  res = carstm_model( p=p,  DS="carstm_randomeffects" ) 

  if (toplot=="random_spatial") {
    vn=c(  "space", "re_total" ) 
    
    if (is.null(fn)) {
      
      fnx = paste(paste0(vn, collapse="_"), "png", sep=".")
      if (!is.null(fn_root_prefix)) fnx = paste(fn_root_prefix, fnx, sep="." )

      fn = file.path( outputdir, fnx  )
    }

    if (is.null(annotation)) {
      annotation=paste( p$carstm_model_label, "persistent spatial effect" )
    }

    toplot = carstm_results_unpack( res, vn )
    
    if (is.null(brks)) {
      qn = quantile( tf(toplot[,"mean"]), probs=probs, na.rm=TRUE ) 
      brks = pretty( qn )
    }

    carstm_map( res=res, vn=vn, 
      breaks = brks, colors=colors, additional_features=additional_features,
      legend.position.inside=legend.position, transformation=tf,
      annotation=annotation, 
      outfilename=fn
    )
    return("done")
  }


  if (toplot=="predictions") {
    res = carstm_model( p=p,  DS="carstm_predictions" ) 
    vn="predictions"
    toplot = carstm_results_unpack( res, vn )
  
    if (p$dimensionality == "space") {

      if (is.null(brks)) {
        qn = quantile( tf(toplot[,"mean"]), probs=probs, na.rm=TRUE ) 
        brks = pretty( qn )
      }
      
      if (is.null(fn)) {
        fnx = paste(paste0(vn, collapse="_"), "png", sep=".")
        if (!is.null(fn_root_prefix)) fnx = paste( fn_root_prefix, fnx, sep="." )
        fn = file.path( outputdir, fnx) 
      }

      if (is.null(annotation)) {
        annotation=paste( p$carstm_model_label, vn )
      }

      carstm_map( res=res, vn=vn,  tmatch=tmatch, umatch=as.character(u),
        breaks = brks, colors=colors, additional_features=additional_features,
        legend.position.inside=legend.position, transformation=tf,
        annotation=annotation, 
        outfilename=fn
      )

    } else if (p$dimensionality == "space_time") {

      for (y in res$time_name){

        if (is.null(brks)) {
          qn = quantile( tf(toplot[,,"mean"]), probs=probs, na.rm=TRUE )
          brks = pretty(  qn  )
        }

        tmatch = as.character(y) 
        u = res$cyclic_id[7]
        fn = file.path( outputdir,  paste(fn_root_prefix, paste0(vn, collapse="_"), tmatch, "png", sep=".") )
        annotation = paste( p$carstm_model_label, "  ", paste0(tmatch, collapse="-") )

        carstm_map( res=res, vn=vn,  tmatch=tmatch, umatch=as.character(u),
          breaks = brks, colors=colors, additional_features=additional_features,
          legend.position.inside=legend.position, transformation=tf,
          annotation=annotation, 
          outfilename=fn
        )
      }
 
    } else if (p$dimensionality == "space_time_season") {

      for (y in res$time_name){
      for ( u in res$cyclic_name  ){
        tmatch = as.character(y) 
        umatch = as.character(u)
        u = res$cyclic_id[7]
        fn = file.path( outputdir,  paste(fn_root_prefix, paste0(vn, collapse="_"), tmatch, umatch, "png", sep=".") )
        annotation = paste( p$carstm_model_label, "  ", paste0(tmatch, collapse="-") )
        carstm_map( res=res, vn=vn,  tmatch=tmatch, umatch=umatch,
          breaks = brks, colors=colors, additional_features=additional_features,
          legend.position.inside=legend.position, transformation=tf,
          annotation=annotation, 
          outfilename=fn
        )
      }}
 
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
