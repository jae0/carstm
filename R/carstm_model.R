
carstm_model = function( p, M=NULL, DS="redo", ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  sppoly = areal_units( p=p )  # required by car fit
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

  fn_fit = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
  outputdir = dirname(fn_fit)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  if (DS=="carstm_modelled_fit") {
    if (file.exists(fn_fit)) {
      load( fn_fit )
      return( fit )
    }
  }

  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )

  # INLA does not like duplicates ... causes optimizer to crash frequently
  eps = exp( log( .Machine$double.eps ) / 2)  # ~ 1.5e-8
  M[, p$variabletomodel]  = M[, p$variabletomodel]  + runif( nrow(M), -eps, eps )

  if (exists("data_transformation", p)) M[, p$variabletomodel]  = p$data_transformation$forward( M[, p$variabletomodel] ) # make all positive


  mrange = NULL
  # get hyper param scalings
  if ( grepl("inla", p$carstm_modelengine) ) {

    # hyperparms
    j = which( is.finite(M[,p$variabletomodel]) )
    mrange = range( M[ j, p$variabletomodel ]  )  # on data scale not internal
    
    if (!exists("p$carstm_model_family", p )) p$p$carstm_model_family = "normal"

    if ( grepl( ".*lognormal", p$carstm_model_family)) {
      m = log( M[ j, p$variabletomodel ])
    } else if ( grepl( ".*poisson", p$carstm_model_family)) {
      m = log( M[ j, p$variabletomodel ] / M[ j, "data_offset" ]  )
      mrange = range( M[ j, p$variabletomodel ]/ M[ j, "data_offset" ]  )  # on data scale not internal
      mrange = mrange * median(M[ M$tag=="predictions", "data_offset" ] )
    } else if ( grepl( "*.binomial", p$carstm_model_family)) {
      m = M[ j, p$variabletomodel ]
    } else {
      m = M[,p$variabletomodel]
    }

    H = carstm_hyperparameters( sd(m), alpha=0.5, median(m) )
    m = NULL

  
    p = parameters_add_without_overwriting( p,
      options.control.family = inla.set.control.family.default()
    )



    gc()

    fit  = NULL

      if (!exists("options.control.inla", p )) p$options.control.inla = list(
        list( optimise.strategy="smart", stupid.search=TRUE, strategy="adaptive", h=0.0001, cmin=0), # default h=0.02
        list( optimise.strategy="smart", stupid.search=TRUE, strategy="adaptive", h=0.001, cmin=0), # default h=0.02 ?or 0.01
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.05, cmin=0, tolerance=1e-9),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.1, cmin=0, tolerance=1e-9),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive"), # default h=0.02
        list( optimise.strategy="smart", h=0.1 ),
        list( optimise.strategy="smart", h=0.2 ),
        list( optimise.strategy="smart", h=0.4 ),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="laplace", fast=FALSE, step.factor=0.1),
        list( stupid.search=TRUE, h=0.001, cmin=0)
      )


      for ( civ in 1:length(p$options.control.inla)) {
        res = try( inla( p$carstm_model_formula , data=M, family=p$carstm_model_family,
          control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE, config=TRUE),
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list(compute=FALSE, link=1 ),
          control.fixed= list(mean.intercept=0, prec.intercept=0.001, mean=0, prec=0.001),
          control.family = p$options.control.family,
          control.inla   = p$options.control.inla[[civ]],
          verbose=TRUE
        ))
        if (!inherits(res, "try-error" )) break()
      }
      if (inherits(res, "try-error" )) {
        message("If you are using MSWindows and you get a popup complaining about inla stopped working,")
        message("you can prevent this by setting the flag in the following link to 1, using regedit. Be careful.")
        message("e.g., see: https://monitormyweb.com/guides/how-to-disable-stopped-working-message-in-windows")
        stop( "solution did not converge")
      }
  
  }
  if ( grepl("glm", p$carstm_modelengine) ) {
    res = try( glm( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )
  }

  if ( grepl("gam", p$carstm_modelengine) ) {
    res = try( gam( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )
  }

  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")

  save( fit, file=fn_fit, compress=TRUE )

  res = carstm_summary( p=p, operation="compute", fit=fit, M=M, sppoly=sppoly, mrange=mrange )

  if (DS!="carstm_modelled_fit") fit = fn_fit

  return(fit)
}
