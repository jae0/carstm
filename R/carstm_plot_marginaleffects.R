carstm_plot_marginaleffects = function( p, outputdir=tempdir(), fn_root_prefix=NULL, prefix = "" ) {

  require(ggplot2)

  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
 
  
  if (prefix=="") {
    if (!is.null(fn_root_prefix)) {

      if ( grepl("probability", fn_root_prefix,  ignore.case =TRUE) | grepl("presence_absence", fn_root_prefix,  ignore.case =TRUE)  ) {
        prefix = "habitat"
      } else if ( grepl("numerical", fn_root_prefix,  ignore.case =TRUE ) ) {
        prefix = "number"
      } else if ( grepl("biomass", fn_root_prefix,  ignore.case =TRUE ) ) {
        prefix = "biomass"
      } else if ( grepl("weight", fn_root_prefix,  ignore.case =TRUE) |  grepl("meansize", fn_root_prefix,  ignore.case =TRUE)  ) {
        prefix = "weight"
      } else {
        prefix = fn_root_prefix
      }
    }
  }
  
  res = carstm_model( p=p,  DS="carstm_summary" )  # parameters in p and direct summary
  names_res = names(res)

  # add more as required: 

  # annual component 
  i = grep("\\<time\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    res$time$yr = as.numeric(p$time_name[res$time$ID])
    plt = ggplot( res$time, aes(x=yr, y=mean)) +  geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=0.5) +
      labs(title="Annual effect", x="Year", y =prefix) +
      theme_light( base_size=22) 
    print(plt)    
    fn = "time.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)    
  }
  
  # seasonal component 
  i = grep("\\<cyclic\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    res$cyclic$seas = as.numeric( p$cyclic_name[res$cyclic$ID] )
    plt = ggplot( res$cyclic, aes(x=seas, y=mean) ) + geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=0.05) +
      labs(title="Seasonal effect", x="Year (fraction)", y =prefix) +
      theme_light( base_size=22) 
    fn = "cyclic.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)    
  }

  # relationship with depth
  i = grep("\\<z\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    # vn = paste( "inla.group(z, method = \"quantile\", n = ", n_discret, ")", sep="")
    plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0 ) +
      labs(title="Depth effect", x="Depth (log; m)", y =prefix) +
      theme_light( base_size=22) 
    fn = "depth.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)    
  }


  i = grep("\\<t\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    # vn = paste( "inla.group(t, method = \"quantile\", n = ", n_discret, ")", sep="")
    plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
      labs(title="Temperature effect", x="Bottom temperature (deg C)", y =prefix) +
      theme_light( base_size=22) 
    fn = "temperature.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)  
  }


  i = grep("\\<pca1\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    # vn= paste( "inla.group(pca1, method = \"quantile\", n = ", n_discret, ")", sep="")
    plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
      labs(title="PCA1 effect", x="PCA1", y =prefix) +
      theme_light( base_size=22) 
    fn = "pca1.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt) 
  }

  i = grep("\\<pca2\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    # vn = paste( "inla.group(pca2, method = \"quantile\", n = ", n_discret, ")", sep="")
    plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
      labs(title="PCA2 effect", x="PCA2", y =prefix) +
      theme_light( base_size=22) 
    fn = "pca2.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)    
  }


  i = grep("\\<substrate.grainsize\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i]
    # vn = paste( "inla.group(substrate.grainsize, method = \"quantile\", n = ", n_discret, ")", sep="") 
    plt = ggplot( res[[vn]], aes(x=ID, y=mean ))+ geom_line(color="gray", linewidth=1.75) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=5) +
      labs(title="Depth effect", x="Depth (m)", y =prefix) +
      theme_light( base_size=22) 
    fn = "substrate.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)    
  }


  i = grep("\\<gear\\>", names_res )
  if (length(i)==1) {
    vn = names_res[i] 
    gears = c("Western IIA", "Yankee #36", "US 4seam 3beam",  "Engle", "Campelen 1800", "Nephrops" )
    res$gear$grs = gears[res$gear$ID]
    plt = ggplot( res$gear, aes(x=grs, y=mean)) + geom_point( size=3, color="slategray") +
      geom_errorbar(aes(ymin=quant0.025, ymax=quant0.975), color="slategray",  linewidth=1.0, width=0.5) +
      labs(title="Gear effect", x="Year", y =prefix) +
      theme_light( base_size=22) 
    fn = "gear.png"
    if (prefix != "") fn = paste(prefix, fn, sep="_")
    fn_plt = file.path( outputdir, fn) 
    ggsave(filename=fn_plt, plot=plt, dpi=144, width=12, height = 8)
    print(fn_plt)
    print(plt)    
  }   

 
}
