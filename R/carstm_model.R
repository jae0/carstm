
carstm_model = function( p, M=NULL, DS="redo", improve.hyperparam.estimates=FALSE, file_compress_method=FALSE, ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  sppoly = areal_units( p=p )  # required by car fit
  areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

  fn_fit = carstm_filenames( p=p, returntype="carstm_modelled_fit", areal_units_fn=areal_units_fn )
  fn_res = carstm_filenames( p=p, returntype="carstm_modelled_summary", areal_units_fn=areal_units_fn )
  outputdir = dirname(fn_fit)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
  
  fit = NULL
  if (DS=="carstm_modelled_fit") {
    if (file.exists(fn_fit)) load( fn_fit )
    if (is.null(fit)) message("carstm modelled fit not found.")
    return( fit )
  }

  O = NULL
  if (DS=="carstm_modelled_summary") {  # carstm_model.*carstm_modelled
    if (file.exists(fn_res)) load( fn_res)
    if (is.null(O)) message("carstm summary not found.")
    return( O )
  }
 
  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )


  vn = p$variabletomodel

  if (exists("data_transformation", p)) M[, vn]  = p$data_transformation$forward( M[, vn] ) # make all positive

  mrange = NULL
  if ( grepl("inla", p$carstm_modelengine) ) {

    if (!exists("carstm_model_family", p )) p$carstm_model_family = "gaussian"
    
    if ( p$carstm_model_family == "lognormal" ) {
      lnk_function = inla.link.log
    } else if ( grepl( ".*poisson", p$carstm_model_family)) {
      lnk_function = inla.link.log
    } else if ( p$carstm_model_family =="binomial" )  {
      lnk_function = inla.link.logit
    } else {
      lnk_function = inla.link.identity
    }


    # get hyper param scalings
    j = which( is.finite(M[,vn]) )
    if ( grepl( "data_offset", p$carstm_model_formula)) {
      m = lnk_function( M[ j, vn ] / M[ j, "data_offset" ] )
    } else {
      m = lnk_function( M[ j, vn ] )
    }


    ii = which(is.finite(m))
    mrange = range( m[ii]  )  # on data scale not internal
    H = carstm_hyperparameters( sd(m[ii] ), alpha=0.5, median(m[ii], na.rm=TRUE) )
    m = NULL

    p = parameters_add_without_overwriting( p,
      options.control.family = inla.set.control.family.default()
    )


    gc()

    fit  = NULL

      if (!exists("options.control.inla", p )) p$options.control.inla = list(
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive"), # default h=0.02
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.05, cmin=0, tolerance=1e-9),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.1, cmin=0),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.001, cmin=0), # default h=0.02 ?or 0.01
        list( optimise.strategy="smart", stupid.search=TRUE, strategy="adaptive", h=0.0001, cmin=0), # default h=0.02
        list( optimise.strategy="smart", h=0.2 ),
        list( optimise.strategy="smart", h=0.4 ),
        list( stupid.search=TRUE, fast=FALSE, step.factor=0.1),
        list( stupid.search=TRUE, cmin=0)
      )


      for ( civ in 1:length(p$options.control.inla)) {
        fit = try( inla( 
          p$carstm_model_formula , 
          data=M, 
          family = p$carstm_model_family,
          control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE),
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list(compute=TRUE, link=1 ),
          # control.fixed= list(mean.intercept=0, prec.intercept=0.001, mean=0, prec=0.001),
          control.family = p$options.control.family,
          control.inla   = p$options.control.inla[[civ]],
          verbose=TRUE
        ))
        if (!inherits(fit, "try-error" )) break()
      }

      if (inherits(fit, "try-error" )) {
        message("If you are using MSWindows and you get a popup complaining about inla stopped working,")
        message("try setting the flag in the following link to 1, using regedit. Be careful.")
        message("e.g., see: https://monitormyweb.com/guides/how-to-disable-stopped-working-message-in-windows")
        stop( "solution did not converge")
      }

      # to improve hyper param estimates..
      if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.25, diff.logdens=18 )  # get improved estimates for the hyperparameters

  }

  if ( grepl("glm", p$carstm_modelengine) ) {
    fit = try( glm( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )
  }

  if ( grepl("gam", p$carstm_modelengine) ) {
    fit = try( gam( p$carstm_model_formula , data=M, family=p$carstm_model_family ) )
  }

  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")

  message( "Saving carstm fit (this can be slow): ", fn_fit )

  save( fit, file=fn_fit, compress=file_compress_method )


  # ----------------
  # summarize
  # do the computations here as fit can be massive ... best not to copy, etc ..
  message( "Computing summaries ..." )


  # results go here
  O = list( M=M, dimensionality=p$aegis_dimensionality, summary=summary(fit), sppoly=sppoly, fn_res=fn_res )
  
  # row indices for predictions
  S = sppoly[["AUID"]]
  nAUID = length(S)

  if ( p$aegis_dimensionality == "space") {
    ipred = which( M$tag=="predictions" & M$AUID %in% S )  # filter by S and T in case additional data in other areas and times are used in the input data
    mfrom = list( S=M$AUID[ipred] )
    mto   = list( S=S )
  }

  if ( p$aegis_dimensionality == "space-year") {
    T = as.character( p$yrs )
    M$year = as.character(M$year)
    # filter by AUID and years in case additional data in other areas and times are used in the input data
    ipred = which( M$tag=="predictions" & M$AUID %in% S & M$year %in% T )  
    mfrom = list( S=M$AUID[ipred], year=M$year[ipred] )
    mto   = list( S=S, T=T )
  }

  if ( p$aegis_dimensionality == "space-year-season") {
    T = as.character( p$yrs )
    U = as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) )

    M$year = as.character(M$year)
    M$dyear = as.character( discretize_data( M$dyear, p$discretization[["dyear"]] ) )

    ipred = which( M$tag=="predictions" & M$AUID %in% S & M$year %in% T )  # ignoring U == predict at all seassonal components ..
    mfrom = list( S=M$AUID[ipred], T=M$year[ipred], U=M$dyear[ipred] )
    mto   = list( S=S, T=T, U=U )
  }


  if ( grepl("glm", p$carstm_modelengine) |  grepl("gam", p$carstm_modelengine) ) {

    if ( p$aegis_dimensionality == "space") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      MM = paste( M$AUID, M$year, sep=":")
      ipred = match(withsolutions, MM)
      mfrom = list( S=M$AUID[ipred] )
    }

    if ( p$aegis_dimensionality == "space-year") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( M$AUID, M$year, sep=":")
      ipred = ipred[ match(withsolutions, MM) ]
      mfrom = list( S=M$AUID[ipred], year=M$year[ipred]  )
    }

    if ( p$aegis_dimensionality == "space-year-dyear") {
      withsolutions = names( which(is.finite(fit$coefficients)) )
      withsolutions = withsolutions[ grep( "AUID.*:year", withsolutions ) ]
      withsolutions = gsub("AUID", "", withsolutions )
      withsolutions = gsub("dyear", "", withsolutions )
      withsolutions = gsub("year", "", withsolutions )
      MM = paste( M$AUID, M$year, M$dyear, sep=":")
      ipred = match(withsolutions, MM)
      mfrom = list( S=M$AUID[ipred], year=M$year[ipred], dyear=M$dyear[ipred] )
    }

    if ( grepl("glm", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[ipred,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
      if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - 1.96*preds$se.fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
      if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + 1.96*preds$se.fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
      if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive
    }


    if ( grepl("gam", p$carstm_modelengine) ) {

      preds = predict( fit, newdata=M[ipred,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      vn =  paste(p$variabletomodel, "predicted", sep=".")
      input = preds$fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
      if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_se", sep=".")
      input = preds$se.fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )

      vn =  paste(p$variabletomodel, "predicted_lb", sep=".")
      input = preds$fit - 1.96*preds$se.fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
      if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive

      vn =  paste(p$variabletomodel, "predicted_ub", sep=".")
      input = preds$fit + 1.96*preds$se.fit
      O[[vn]] = reformat_to_array( input =input, matchfrom=mfrom, matchto=mto )
      if ( grepl( ".*poisson", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive
    }
  }


  if ( grepl("inla", p$carstm_modelengine) ) {


    invlink = function(x) lnk_function( x,  inverse=TRUE ) 

    summary_inv = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )

    summary_inv_prec = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(pmax(y,1e-12)), x) , silent=TRUE  )
    # summary_inv_prec_1024 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=1024L) , silent=TRUE  )
    # summary_inv_prec_512 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=512L) , silent=TRUE  )

    list_simplify = function(x) as.data.frame( t( as.data.frame( x )))

    exceedance_prob = function(x, threshold)  {1 - inla.pmarginal(q = threshold, x)}

    deceedance_prob = function(x, threshold)  { inla.pmarginal(q = threshold, x)}

    tokeep =  c("mean", "sd", "quant0.025", "quant0.5", "quant0.975")

    exceedance_threshold=1 
    deceedance_threshold=1

    if (exists("deceedance_threshold", p)) deceedance_threshold=p[["deceedance_threshold"]]
    if (exists("exceedance_threshold", p)) exceedance_threshold=p[["exceedance_threshold"]]

    O[["summary"]] = list()
    O[["random"]] = list()


    toget = c("summary", "random_spatial", "random_spatiotemporal" , "predictions", "predictions_adjusted")

    if ( "summary" %in% toget) {
        
        # parameters
        # back-transform from marginals
        W = cbind ( t (sapply( fit$marginals.fixed, FUN=summary_inv ) ) )

        O[["summary"]][["fixed_effects"]] = W [, tokeep, drop =FALSE]

        # hyperpar (variance components)   
        j = grep( "^Precision.*", rownames(fit$summary.hyperpar), value=TRUE )
        if (length(j) > 0) {
          precs = try( list_simplify( sapply( fit$marginals.hyperpar[j], FUN=summary_inv_prec ) ), silent=TRUE )  # prone to integration errors ..
        #  if (inherits(precs, "try-error")) precs = try( list_simplify( sapply( fit$marginals.hyperpar[j], FUN=summary_inv_prec_1024 ) ), silent=TRUE )
           if (inherits(precs, "try-error")) message( "Try alt parameterization for prec -> sd or smaller number of n or masking negative values")

          # precs[,"mode"] =  1/sqrt( fit$summary.hyperpar[j,"mode"]  )
          toadd = setdiff( colnames(O[["summary"]]), colnames(precs) ) 
          precs[,toadd] = NA
          rownames(precs) = gsub("Precision for", "SD", rownames(precs) )
          O[["summary"]][["random_effects"]] = precs[, tokeep, drop =FALSE]
        }

        j = grep( ".*Rho.*", rownames(fit$summary.hyperpar), value=TRUE )
        if (length(j) > 0) {
          rhos = list_simplify( sapply( fit$marginals.hyperpar[j], FUN=function(x) inla.zmarginal( x, silent=TRUE  ) ) ) 
          #  rhos[,"mode"] = sapply( fit$marginals.hyperpar[j], FUN=function(x) inla.mmarginal( x )) 
          O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], rhos[, tokeep, drop =FALSE] )
        }

        # update phi's
        j = grep( "^Phi.*", rownames(fit$summary.hyperpar), value=TRUE )
        if (length(j) > 0) {
          phis = fit$summary.hyperpar[j, ,drop=FALSE]
          colnames(phis) = c( tokeep, "mode")
          O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], phis[, tokeep, drop =FALSE] )		
        }
                
    }



  if ("random_spatial" %in% toget) {
    
      if (exists("summary.random", fit)) {

        if ( exists("space", fit$summary.random) ) {
          
          O[["random"]] [["space"]] = list()  # space as a main effect

          # bym2 output index locations
          s_bym2 = expand.grid( space=S, type = c("iid", "bym"), stringsAsFactors =FALSE )

          u = which(s_bym2$type=="bym")
          v = which(s_bym2$type=="iid")
          
          g = fit$marginals.random[["space"]]
          m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  
          # m = list_simplify ( sapply( g, summary_inv ) ) 
          
          #  iid main effects  
          W = array( NA,  dim=c( length( S), length(names(m)) ), dimnames=list( space=O[["space_name"]], stat=names(m) ) ) 
          for (k in 1:length(names(m))) {
            W[,k] = reformat_to_array( input = unlist(m[v,k]), matchfrom = list( space=s_bym2[["space"]][v] ), matchto = list( space=S  ) )
          }
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random"]] [["space"]] [["iid"]] = W [, tokeep, drop =FALSE]

          #  spatial main effects  
          W = array( NA, dim=c( length( S), length(names(m)) ), dimnames=list( space=O[["space_name"]], stat=names(m) ) ) 
          for (k in 1:length(names(m))) {
            W[,k] = reformat_to_array(  input = unlist(m[u,k]), matchfrom = list( space=s_bym2[["space"]][u]  ), matchto = list( space=S  ) )
          } 
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random"]] [["space"]] [["bym"]] = W [, tokeep, drop =FALSE]

          # space == iid + bym combined:
          aa = inla.posterior.sample( nposteriors, fit, selection=list(space=0), add.names =FALSE )  # 0 means everything matching space
          g = sapply( aa, function(x) invlink(x$latent[v] + x$latent[u]) )
          mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
          mm = apply( g, 1, mean, na.rm =TRUE)
          ms = apply( g, 1, sd, na.rm =TRUE)
          W = cbind(mm, ms, mq)
          attr(W, "dimnames") = list( space=O[["space_name"]], stat=tokeep  ) 
          O[["random"]] [["space"]] [["combined"]] = W
          

          if (!is.null(exceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( 
              input = m, 
               matchfrom = list( space=s_bym2[["space"]] ), 
              matchto = list( space=S  )
            )
            dimnames( W )[["space"]] = O[["space_name"]]
            O[["random"]] [["space"]] [["exceedance"]] = W
          }

          if (!is.null(deceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( 
              input = m, 
               matchfrom = list( space=s_bym2[["space"]] ), 
              matchto = list( space=S  )
            )
            dimnames( W )[["space"]] = O[["space_name"]]
            O[["random"]] [["space"]] [["deceedance"]] = W
          }

        }
      }

  }


  if ("random_spatiotemporal"  %in% toget ) {

        #  2. iid interaction effects grouped by time
        if ( exists( "space_time", fit$summary.random) ) {

          O[["random"]] [["spacetime"]] = list()

          st_bym2 = expand.grid( space=S, type = c("iid", "bym"), time=as.character(T), stringsAsFactors =FALSE )
     
          u = which(st_bym2$type=="bym")
          v = which(st_bym2$type=="iid")
          
          g = fit$marginals.random[["space_time"]]
          m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  
          # m = list_simplify ( sapply( g, summary_inv ) ) 
     
          #  spatiotemporal interaction effects  bym
          W = array( NA, dim=c( length( S), length(T), length(names(m)) ), dimnames=list( space=O[["space_name"]], time=T, stat=names(m) ) ) 
          for (k in 1:length(names(m))) {
            W[,,k] = reformat_to_array(  input = unlist(m[u,k]), matchfrom = list( space=st_bym2[["space"]][u],  time=st_bym2[["time"]][u]  ), matchto = list( space=S, time=T  ) )
          } 
          W = W [,, tokeep, drop =FALSE]
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random"]] [["spacetime"]] [["bym"]] = W


          #  spatiotemporal interaction effects  iid
          W = array( NA, dim=c( length( S), length(T), length(names(m)) ), dimnames=list( space=O[["space_name"]], time=T, stat=names(m) ) ) 
          for (k in 1:length(names(m))) {
            W[,,k] = reformat_to_array(  input = unlist(m[v,k]), matchfrom = list( space=st_bym2[["space"]][v],  time=st_bym2[["time"]][v]  ), matchto = list( space=S, time=T  ) )
          } 
          W = W [,, tokeep, drop =FALSE]
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random"]] [["spacetime"]] [["iid"]] = W


          # space == iid + bym combined:
          aa = inla.posterior.sample( nposteriors, fit, selection=list(space_time=0), add.names =FALSE )  # 0 means everything matching space
          g = sapply( aa, function(x) invlink(x$latent[v] + x$latent[u]) )
          mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
          mm = apply( g, 1, mean, na.rm =TRUE)
          ms = apply( g, 1, sd, na.rm =TRUE)
          m = data.frame( cbind(mm, ms, mq) )
          names(m) = tokeep
          W = array( NA, dim=c( length( S), length(T), length(names(m)) ), dimnames=list( space=O[["space_name"]], time=T, stat=names(m) ) ) 
          for (k in 1:length(names(m))) {
            W[,,k] = reformat_to_array(  input = m[,k], matchfrom = list( space=st_bym2[["space"]][v],  time=st_bym2[["time"]][v]  ), matchto = list( space=S, time=T  ) )
          } 
          W = W [,, tokeep, drop =FALSE]
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random"]] [["spacetime"]] [["combined"]] = W
          

          if (!is.null(exceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( 
              input = m, matchfrom = list( space=st_bym2[["space"]][v],  time=st_bym2[["time"]][v]  ), matchto = list( space=S, time=T  )
            )
            dimnames( W )[["space"]] = O[["space_name"]]
            O[["random"]] [["spacetime"]] [["exceedance"]] = W
          }

          if (!is.null(deceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( 
              input = m, matchfrom = list( space=st_bym2[["space"]][v],  time=st_bym2[["time"]][v]  ), matchto = list( space=S, time=T  )
            )
            dimnames( W )[["space"]] = O[["space_name"]]
            O[["random"]] [["spacetime"]] [["deceedance"]] = W
           }
        }
  }



  if ("space_iid"  %in% toget ) {

      # IID random effects .. not used .. might need tweaking
    
        if (exists("space", fit$summary.random)) {
          # bym2 output index locations
          s_bym2 = expand.grid( space=S, type = c("iid", "bym"), stringsAsFactors =FALSE )

          i = which(s_bym2$type=="iid") 
    
          W = reformat_to_array( 
            input = fit$summary.random$space[ i, "mean" ], 
            matchfrom = list( space=s_bym2[["space"]][i], time=s_bym2[["time"]][i]), 
            matchto =   list( space=S, time=T )
          )
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random_space_iid"]] = W


          W = reformat_to_array( 
            input = fit$summary.random$space[ i, "sd" ], 
            matchfrom = list( space=s_bym2[["space"]][i], time=s_bym2[["time"]][i]), 
            matchto =   list( space=S, time=T )
          )
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random_space_iid_sd"]] = W
        }

    
        if (exists("space_time", fit$summary.random)) {

          st_bym2 = expand.grid( space=S, type = c("iid", "bym"), time=as.character(T), stringsAsFactors =FALSE )

          i = which(st_bym2$type=="iid") 
    
          W = reformat_to_array( 
            input = fit$summary.random$space_time[ i, "mean" ], 
            matchfrom = list( space=st_bym2[["space"]][i], time=st_bym2[["time"]][i] ), 
            matchto =   list( space=S, time=T )
          )
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random_space_time_iid"]] = W

          W = reformat_to_array( 
            input = fit$summary.random$space_time[ i, "sd" ], 
            matchfrom = list( space=st_bym2[["space"]][i], time=st_bym2[["time"]][i] ), 
            matchto =   list( space=S, time=T )
          )
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["random_space_time_iid_sd"]] = W
       }

  
  }

  if ("predictions"  %in% toget ) {
    
      # revert these to character as they are used for matching
      SP = levels( O[["space_factor"]] ) [ O[["data"]][,"space"] ]
      TP = levels( O[["time_factor"]] ) [ O[["data"]][,"time"] ]

      i = which(  O$data$tag=="predictions" & SP %in% S & TP %in% T )  

      # predictions
      if (exists("summary.fitted.values", fit)) {

        m = list_simplify ( sapply( fit$marginals.fitted.values[i], inla.zmarginal, silent=TRUE ) )  # already backfransformed by link=1
        W = array( NA,  dim=c( length( S), length(T), length(names(m)) ),  dimnames=list( space=O[["space_name"]], time=T, stat=names(m) ) ) 
        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array( 
            input = unlist(m[,k]), 
            matchfrom = list( space=SP[ i ], time=TP[i] ), 
            matchto =   list( space=S, time=T )
          )
        } 
        O[["predictions"]] = W[,, tokeep, drop =FALSE]

      }

  }

  if ("predictions_adjusted"  %in% toget ) {
      
        o = invlink( O$data[[ O$vnO ]] )  # offsets

        g = fit$marginals.fitted.values
        for (e in i) g[[e]][,1] = g[[e]][,1] / o[e] 

        m = list_simplify ( sapply( g[i], inla.zmarginal, silent=TRUE ) )  # already backfransformed by link=1
        W = array( NA, dim=c( length(S), length(T), length(names(m)) ),  dimnames=list( space=O[["space_name"]], time=T, stat=names(m) ) ) 
        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array( 
            input = unlist(m[,k]), 
            matchfrom = list( space=SP[ i ], time=TP[i] ), 
            matchto =   list( space=S, time=T )
          )
        } 
        dimnames( W )[["space"]] = O[["space_name"]]
        O[["predictions_adjusted"]] = W[,, tokeep, drop =FALSE]

        # offsets on user scale
        W = reformat_to_array( 
          input = o[i], 
          matchfrom = list( space=SP[ i ], time=TP[i] ), 
          matchto =   list( space=S, time=T )
        )
        dimnames( W )[["space"]] = O[["space_name"]]
        O[["prediction_offsets"]] = W

        O[["predictions_adjusted_direct"]] = O[["predictions"]][,,"mean"] /  O[["prediction_offsets"]]   ## relative rate as a check

        if (!is.null(exceedance_threshold)) {
          m = list_simplify ( sapply( g[i], FUN=exceedance_prob, threshold=exceedance_threshold ) )
          W = reformat_to_array( 
            input = unlist(m ), 
            matchfrom = list( space=SP[ i ], time=TP[i] ), 
            matchto =   list( space=S, time=T )
          )
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["exceedance_probability"]] = W
        }

        if (!is.null(deceedance_threshold)) {
          m = list_simplify ( sapply( g[i], FUN=deceedance_prob, threshold=deceedance_threshold ) )
          W = reformat_to_array( 
            input = unlist(m ), 
            matchfrom = list( space=SP[ i ], time=TP[i] ), 
            matchto =   list( space=S, time=T )
          )
          dimnames( W )[["space"]] = O[["space_name"]]
          O[["deceedance_probability"]] = W
        }

        if ( length(T) > 2 ) {
          if ( exists("predictions_adjusted", O ) ) {
            ti = as.numeric( T )
            lmslope = function( x ) summary( lm(x~ti) )$coefficients["ti",1:2] 
            W = t ( apply( O[["predictions_adjusted"]][,,"mean"], 1, lmslope ) )  # relative rate per year
            dimnames( W )[["space"]] = O[["space_name"]]
            O[["time_slope"]] = W
          }
        }
  
    }

    if (exists("summary.random", fit)) {

      if (exists("auid", fit$summary.random)) {

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single nonspatial effect (no grouping across time)
          Z = expand.grid( S=S, type = c("nonspatial", "spatial") )
          O$i_nonspatial = which(Z$type=="nonspatial")
          O$ns_matchfrom = list( S=Z$S[O$i_nonspatial]  )
          O$ns_matchto   = list( S=S  )

        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          #  nonspatial effects grouped by year
          Z = expand.grid( S=S, type = c("nonspatial", "spatial"), T=p$yrs )
          O$i_nonspatial = which(Z$type=="nonspatial")
          O$ns_matchfrom = list( S=Z$S[O$i_nonspatial], year=Z$T[O$i_nonspatial] )
          O$ns_matchto   = list( S=S, T=T  )
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # nonspatial at all time slices
          Z = expand.grid( S=S, type = c("nonspatial", "spatial"), T=p$yrs, U=p$dyears )
          O$i_nonspatial = which(Z$type=="nonspatial")
          O$ns_matchfrom = list( S=Z$S[O$i_nonspatial], year=as.character(Z$T[O$i_nonspatial]), dyear=as.character( discretize_data( Z$dyear[O$i_nonspatial], p$discretization[["dyear"]] ) ) )
          O$ns_matchto   = list( S=S, T=T, dyear=O$dyear )
        }

        if (nrow(fit$summary.random$auid) == nAUID*2) {
          # a single spatial effect (no grouping across time)
          Z = expand.grid( S=S, type = c("nonspatial", "spatial") )
          O$i_spatial = which(Z$type=="spatial")
          O$sp_matchfrom = list( S=Z$S[O$i_spatial]  )
          O$sp_matchto   = list( S=S  )

        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$ny ) {
          # spatial effects grouped by year
          Z = expand.grid( S=S, type = c("nonspatial", "spatial"), T=p$yrs )
          O$i_spatial = which(Z$type=="spatial")
          O$sp_matchfrom = list( S=Z$S[O$i_spatial], year=Z$T[O$i_spatial] )
          O$sp_matchto   = list( S=S, T=T  )
        } else if (nrow(fit$summary.random$auid) == nAUID*2 * p$nt ) {
          # at every time slice
          Z = expand.grid( S=S, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          O$i_spatial = which(Z$type=="spatial")
          O$sp_matchfrom = list( S=Z$S[O$i_spatial], year=as.character(Z$T[O$i_spatial]), dyear=as.character( discretize_data( Z$dyear[O$i_spatial], p$discretization[["dyear"]] ) ) )
          O$sp_matchto   = list( S=S, T=T, dyear=O$dyear )
        }

      }
    }


    # predictions

    if (exists("summary.fitted.values", fit)) {

      vn = paste( p$variabletomodel, "predicted", sep=".")
      input = fit$summary.fitted.values[ ipred, "mean" ]
      O[[vn]] = reformat_to_array( input=input, matchfrom=mfrom, matchto=mto )

      NA_mask = NULL
      if ( exists("carstm_predict_force_range", p)) {
        if (p$carstm_predict_force_range) {
          if (!is.null(mrange)) {
            # only an issue for factorial models ... missing locations get predicted to absurd values as there is no information .. filter out
            NA_mask = which(O[[vn]]  > max(mrange, na.rm=TRUE) )
          }
        }
      }
      if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive


      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      input = fit$summary.fitted.values[ ipred, "0.025quant" ]
      O[[vn]] = reformat_to_array( input=input, matchfrom=mfrom, matchto=mto )
      if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      input = fit$summary.fitted.values[ ipred, "0.975quant" ]
      O[[vn]] = reformat_to_array( input=input, matchfrom=mfrom, matchto=mto )
      if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA
      if ( grepl( ".*lognormal", p$carstm_model_family)) O[[vn]] = exp(O[[vn]])
      if (exists("data_transformation", p) ) O[[vn]] = p$data_transformation$backward( O[[vn]] ) # make all positive

      vn = paste( p$variabletomodel, "predicted_se", sep=".")
      input = fit$summary.fitted.values[ ipred, "sd" ]
      O[[vn]] = reformat_to_array( input=input, matchfrom=mfrom, matchto=mto )
      if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA

    }

    ## --------- predictions complete ------



    ## --------- start random effects -------

    # match conditions for random effects .. i_preds are locations of predictions in "fit"
    # random effects results ..
    if (exists("summary.random", fit)) {

      if (exists("iid_error", fit$summary.random)) {
        # IID random effects
        vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[ipred, "mean" ]
        O[[vn]] = reformat_to_array( input=input, matchfrom=mfrom, matchto=mto )
        if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA
      }

      if (exists("auid", fit$summary.random)) {

        input = fit$summary.random$auid[ O$i_nonspatial, "mean" ]
        vn = paste( p$variabletomodel, "random_auid_nonspatial", sep=".")
        O[[vn]] = reformat_to_array( input=input, matchfrom=O$ns_matchfrom, matchto=O$ns_matchto )
        if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA
        # carstm_map( O=O, vn=vn, time_match=list(year="2000", dyear="0.85" ) )

        input = fit$summary.random$auid[ O$i_spatial, "mean" ]  # offset structure due to bym2
        vn = paste( p$variabletomodel, "random_auid_spatial", sep=".")
        O[[vn]] = reformat_to_array( input=input, matchfrom=O$sp_matchfrom, matchto=O$sp_matchto )
        if (!is.null(NA_mask)) O[[vn]][NA_mask] = NA
        # carstm_map( O=O, vn=vn, time_match=list(year="2000", dyear="0.85" ) )

      }
    }
  }

  save( O, file=fn_res, compress=file_compress_method )

  message( "carstm summary saved as: ", fn_res )

  print( O$summary )

  # if ( grepl("inla", p$carstm_modelengine) ) {

  #   # a few plots
  #   # copied from INLA::plot.inla and streamlined to reduce access to the fit object which can be huge
  #   print( plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=TRUE ) )

  # } else {
  #   print( plot(fit) )
  # }

  return( fit )
}
