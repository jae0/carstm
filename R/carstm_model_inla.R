
carstm_model_inla = function(
  O, 
  DS = "",
  sppoly =NULL,
  space.id = NULL,
  time.id = NULL,
  cyclic.id=NULL,
  fn_fit=tempfile(pattern="fit_", fileext=".rdata"), 
  fn_res=NULL, 
  compress=TRUE,
  scale_offsets = FALSE,
  redo_fit = TRUE,
  update_results = FALSE,
  toget = c("summary", "fixed_effects", "random_other", "random_spatial", "random_spatiotemporal" , "predictions"), 
  nposteriors=NULL, 
  exceedance_threshold=NULL, 
  deceedance_threshold=NULL, 
  exceedance_threshold_predictions=NULL,
  deceedance_threshold_predictions=NULL,
  improve.hyperparam.estimates=NULL,  ... ) {
  
  if (0) {
    # for debugging
    P = list()
    DS = ""
    sppoly =NULL
    space.id = NULL
    time.id =NULL
    cyclic.id = NULL
    fn_fit=tempfile(pattern="fit_", fileext=".rdata")
    compress="gzip"
    redo_fit = TRUE
    update_results = FALSE
    toget = c("summary", "fixed_effects", "random_other", "random_spatial", "random_spatiotemporal" , "predictions")
    exceedance_threshold=NULL
    deceedance_threshold=NULL
    exceedance_threshold_predictions=NULL
    deceedance_threshold_predictions=NULL
    nposteriors=NULL
    improve.hyperparam.estimates=NULL

    if (0) {
        # usual variable names used in aegis .. char / num
        O$vn$S = "space"  # "space"
        O$vn$T = "time"
        O$vn$U = "cyclic" # "dyri"  # sub annual time

        # alt charracter descrptions of vars (if any, taken from main above if not)
        # O$vn$S0 = "AUID"  # as character
        # O$vn$T0 = "yr"  # as character
        # O$vn$U0 = "dyear"

        # copies of main effects for inla in formulae
        O$vn$ST = "space_time"  # vnST = "space_time" (copy of vnS)
        O$vn$TS = "time_space"  # vnTS = "time_space" (copy of vnT)
    }
  }

  P = list(...)  # INLA options to be passed directly to it  

  if (DS=="modelled_fit") {
    if (!is.null(fn_fit)) {
      fit  = NULL
      message("Loading fit: ", fn_fit )
      if (file.exists(fn_fit)) load( fn_fit )
      if (is.null(fit)) message("modelled fit not found.")
    }
    return( fit )
  }

  if (DS=="modelled_summary") {  # model.*modelled
    if (!is.null(fn_res)) {
      message("Loading  data summary:  ", fn_res )
      O = NULL
      if (file.exists(fn_res)) load( fn_res)
      if (is.null(O)) message(" summary not found.")
      return( O )
    }
  }

  num.threads = "1:1"
  if (exists("num.threads", O)) num.threads = O[["num.threads"]]
  if (exists("num.threads", P)) num.threads = P[["num.threads"]]
  inla.setOption(num.threads=num.threads)

  num.cores =  as.numeric( unlist(strsplit(num.threads, ":") )[1] )
  
  # local functions
  apply_generic = function(...)  mclapply(...,   mc.cores=num.cores ) # drop-in for lapply
  apply_simplify = function(...) simplify2array(mclapply(...,  mc.cores=num.cores ))  # drop in for sapply
  #apply_generic = function(...)  lapply(...  ) # drop-in for lapply
  #apply_simplify = function(...) simplify2array(lappy(...  ))  # drop in for sapply



  outputdir = dirname(fn_fit)
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
  

  # outputs
  fit  = NULL

  if ( is.null(O)){
    if ( update_results) {
      if (!is.null(fn_res)) {
        load (fn_res) 
      } else {
        fit = NULL
        load(fn_fit )
        if (!is.null(fit) ) {
          if (exists("results", fit)) O = fit$results
        } 
      }
    } 
  } 


  if ( inherits(O, "try-error")) O = NULL
  if ( is.null(O)) O = list()  # options
  
  if (is.null(nposteriors))  nposteriors = ifelse( exists("nposteriors", O), O$nposteriors, 5000 )

  if ( !exists("dimensionality", O ) ) O[["dimensionality"]] = "space-time" 

  if ( !exists("verbose", P ) ) P[["verbose"]] = TRUE

  # DATA
  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (!exists("data", P)) {
    if (exists("data", O)) {
      if ( P[["verbose"]]) message( "Data not passed as an argument, using data found in options, O")
      P[["data"]] = O[["data"]]  # priority to P
      O[["data"]] = NULL  # reduce mem usage
    }
  }

  if (class(P[["data"]])=="character") {
    if ( P[["verbose"]]) message( "Data is a function. Running ...")
    P[["data"]] = try( eval(parse(text=P[["data"]]) ) )
  }

  if (inherits(P[["data"]], "try-error"))  P[["data"]] = NULL

  if (is.null(P[["data"]])) stop("Data not found")

  setDF(P[["data"]])


  # FAMILY
  if ( !exists("family", P ) ) {
    if ( exists("family", O) ) {
      if ( P[["verbose"]] ) message( "Family found in options, O" )
      P[["family"]] = O[["family"]]   
    } else {
      P[["family"]] = "gaussian"
    }
  }

  if ( P[["family"]] == "gaussian" ) {
    lnk_function = inla.link.identity
  } else if ( P[["family"]] == "lognormal" ) {
    lnk_function = inla.link.log
  } else if ( grepl( ".*poisson", P[["family"]])) {
    lnk_function = inla.link.log
  } else if ( grepl( ".*binomial", P[["family"]])) {
    lnk_function = inla.link.logit
  } 

  invlink = function(x) lnk_function( x,  inverse=TRUE )

  if ( !exists("formula", P ) ) {
    if ( P[["verbose"]] ) message( "Formula found in options, O" )
    P[["formula"]] = O[["formula"]]
  }

  fm = O[["formula_parsed"]] = parse_formula( P[["formula"]] )
  
  if ( P[["verbose"]] ) {
    message( "Formula parsed as follows. Check in case of errors:" )
    print(fm)
  }


  if (!exists("vn", O)) O[["vn"]] = list() 

  # labels .. not all used but defining the here makes things simpler below
  vnO = NULL  # input data is expected to be on user scale

  if ( exists("O", O[["vn"]]) ) {

    if (exists( O[["vn"]]$O, P[["data"]] )) {
      vnO = O[["vn"]]$O
      if ( P[["verbose"]] ) message( 'setting offset variable from params to: ', O[["vn"]]$O )
    } else {
      if ( P[["verbose"]] ) message( 'offset variable: ', O[["vn"]]$O, ' not found in data')
    }
  }

  if (is.null(vnO)) {
    if ( !is.null(fm$offset_variable) )  {
      if (exists(fm$offset_variable, P[["data"]] )) {
        vnO = fm$offset_variable
        if ( P[["verbose"]] ) message( 'setting offset variable from formula to, ', fm$offset_variable )
      }
    }
  }
  if (is.null(vnO))  vnO="data_offset"  # use generic default in case parsing issues

  vnY = NULL
  if (is.null(vnY)) if ( exists("Y", O[["vn"]])) vnY = O[["vn"]]$Y
  if (is.null(vnY)) if ( exists("variabletomodel", O) ) vnY = O$variabletomodel  
  if (is.null(vnY)) if ( exists("dependent_variable", fm) ) vnY = fm$dependent_variable 

  vnS = ifelse( exists("S", O[["vn"]]), O[["vn"]]$S, "space" )
  vnS0 = ifelse( exists("S0", O[["vn"]]), O[["vn"]]$S0, "space0" ) # local storage/copy as a character
  
  vnT = ifelse( exists("T", O[["vn"]]), O[["vn"]]$T, "time" )
  vnT0 = ifelse( exists("T0", O[["vn"]]), O[["vn"]]$T0, "time0" )

  vnU = ifelse( exists("U", O[["vn"]]), O[["vn"]]$U, "cyclic" )  
  vnUI = ifelse( exists("UI", O[["vn"]]), O[["vn"]]$UI, "cyclic_iid" )  
  vnU0 = ifelse( exists("U0", O[["vn"]]), O[["vn"]]$U0, "cyclic0" )
  
  vnST = ifelse( exists("ST", O[["vn"]]), O[["vn"]]$ST, "space_time" )
  vnTS = ifelse( exists("TS", O[["vn"]]), O[["vn"]]$TS, "time_space" )
  vnSI = ifelse( exists("SI", O[["vn"]]), O[["vn"]]$SI, "space_iid" )
  vnTI = ifelse( exists("TI", O[["vn"]]), O[["vn"]]$TI, "time_iid" )
  vnSTI = ifelse( exists("STI", O[["vn"]]), O[["vn"]]$STI, "space_time_iid" )
  vnTSI = ifelse( exists("TSI", O[["vn"]]), O[["vn"]]$TSI, "time_space_iid" )
  

  if (is.null(sppoly)) if (exists("sppoly", O)) sppoly = O$sppoly
  if (is.null(sppoly)) if (exists("areal_units")) sppoly = areal_units( O ) 
 
  if (is.null(space.id)) {
    if (exists("space.id", O)) {
      if ( P[["verbose"]] ) message( "space.id was not passed, using space.id found in options, O" )
      space.id = O[["space.id"]]
    } 
  }  


  if (is.null(space.id)) {
    if (!is.null(sppoly))  {
      aa = attributes(sppoly)
      if (!is.null(aa)) {
        if ( exists("space.id", aa) ) {
          space.id =  as.character( aa$space.id )
        }
      }
    }
  }


  if (is.null(space.id)) {
    if (!is.null(sppoly))  {
      aa = attributes(sppoly)
      if (!is.null(aa)) {
        if ( exists("region.id", aa) ) {
          space.id =  as.character( aa$region.id )
        }
      }
    }
  }
 
  if (is.null(space.id)) {
    if (!is.null(sppoly))  {
      aa = attributes(sppoly)
      if (!is.null(aa)) {
        if ( exists("nb", aa) ) {
          if ( exists("space.id", attributes(aa$nb) )) {
            space.id =  as.character( slot( aa$nb, "space.id" ) )
          }
        }
      }
    }
  }

  if (is.null(space.id)) {
    if (!is.null(sppoly))  {
      aa = attributes(sppoly)
      if (!is.null(aa)) {
        if ( exists("nb", aa) ) {
          if ( exists("region.id", attributes(aa$nb) )) {
            space.id =  as.character( slot( aa$nb, "region.id" ) )
          }
        }
      }
    }
  }

  if (is.null(space.id)) {
    if (!is.null(sppoly))  {
      aa = attributes(sppoly)
      if (!is.null(aa)) {
        if ( exists("W.nb", aa) ) {
          if ( exists("region.id", attributes(aa$W.nb) )) {
            space.id =  as.character( slot( aa$W.nb, "space.id" ) )
          }
        }
      }
    }
  }

  if (is.null(space.id)) {
    if (!is.null(sppoly))  {
      aa = attributes(sppoly)
      if (!is.null(aa)) {
        if ( exists("W.nb", aa) ) {
          if ( exists("region.id", attributes(aa$W.nb) )) {
            space.id =  as.character( slot( aa$W.nb, "region.id" ) )
          }
        }
      }
    }
  }


  if (is.null(space.id)) {
    if (!is.null(sppoly)) {
      if (exists("AUID", sppoly)) {
        if ( P[["verbose"]] ) message( "space.id was not passed, using space.id found in sppoly" )
        space.id = as.character( sppoly[["AUID"]] )  # the master / correct sequence of the AU's and neighbourhood matrix index values
      }
    }
  }

  if (is.null(space.id)) stop("Not found: space.id is a required variable")
 
  # fiddling of AU and TU inputs: for bym2, etc, they need to be numeric, matching numerics of polygon id ("space.id")
  # convert space and time to numeric codes for INLA
  if (grepl("space", O[["dimensionality"]])) {
    if (!exists(vnS, O)) O[[vnS]] = as.character( space.id )  # this sequence is a master key as index matches nb matrix values
    P[["data"]][,vnS0] = as.character( P[["data"]][,vnS] ) # local copy
    P[["data"]][,vnS] = match( P[["data"]][,vnS0], O[[vnS]] )  # overwrite with numeric values that must match index of neighbourhood matrix
  }
  
browser()

  if ( grepl("time", O[["dimensionality"]]) | grepl("cyclic", O[["dimensionality"]]) ) {
 
    if (any( grepl( vnT, fm$vars )))  {
      
      if (!exists(vnT, O)) if (!is.null(time.id)) O[[vnT]] = as.character( time.id )  # this sequence is a master key
      if (!exists(vnT, O)) if (exists("yrs", O)) O[[vnT]] = as.character( O$yrs ) 
      if (!exists(vnT, O)) if (exists(vnT, P[["data"]])) O[[vnT]] = as.character( sort(unique( P[["data"]][, vnT ]) ) ) 

      P[["data"]][,vnT0] = as.character( P[["data"]][,vnT] ) # a copy for internal matching 
      if (vnT %in% fm$fixed_effects$vn ) {
        P[["data"]][,vnT] = match( P[["data"]][,vnT0], O[[vnT]] ) # convert to data numeric (ie. a numeric factor)
      }
      if (vnT %in% fm$random_effects$vn ) {
        P[["data"]][,vnT] = as.numeric( P[["data"]][,vnT] )  # in case it is sent as a character 
        # O[[vnT]] = as.numeric( O[[vnT]] )
        # nothing to do .. leave alone as numeric
      } 

    }
    # internal vars, for inla
    if (any( grepl( vnST, fm$vars )))  P[["data"]][,vnST] = P[["data"]][,vnS]
    if (any( grepl( vnST, fm$vars )))  P[["data"]][,vnTS] = P[["data"]][,vnT]
    # sub-annual time
    if (any( grepl( vnU, fm$vars )))  {

      if (!exists(vnU, O)) if (!is.null(cyclic.id)) O[[vnU]] = as.character( cyclic.id )  # this sequence is a master key
      if (!exists(vnU, O)) if (exists("dyears", O)) O[[vnU]] = as.character( O$dyears + diff(O$dyears)[1]/2)
      if (!exists(vnU, O)) if (exists(vnU, P[["data"]])) O[[vnU]] = as.character( sort(unique( P[["data"]][, vnU ]) ) ) 

      P[["data"]][,vnU0] = as.character( P[["data"]][,vnU] )  # a copy for internal matching 

      if (vnU %in% fm$fixed_effects$vn ) {
        P[["data"]][,vnU] = match( P[["data"]][,vnU0], O[[vnU]] ) # convert to data numeric (ie. a numeric factor)
      }
      if (vnU %in% fm$random_effects$vn ) {
        P[["data"]][,vnU] = as.numeric( P[["data"]][,vnU] )  # in case it is sent as a character 
        # O[[vnU]] = as.numeric( O[[vnU]] )        # nothing to do .. leave alone as numeric
      } 
    }
  }


  nre = nrow(fm$random_effects)
  if (nre > 0) {
    totestSpace = c( vnST, vnSI, vnSTI ) 
    totestTime = c( vnTS, vnTI, vnTSI  )
    for (re in 1:nre) {
      vntest = fm$random_effects$vn[re]
      if (exists(vntest, P[["data"]] )) next()
      # these are copies so can duplicate columns on the fly
      if (vntest %in% totestSpace) {
        P[["data"]][,vntest] = P[["data"]][,vnS] 
        next()
      }
      if (vntest %in% totestTime)  {
        if ( fm$random_effects$submodel[ grep( vntest, fm$random_effects$vn ) ] %in% c("group", "replicate") ) {
          P[["data"]][,vntest] = match( P[["data"]][,vnT0], O[[vnT]] ) 
        } else {
          P[["data"]][,vntest] = P[["data"]][,vnT] 
        }
        next()
      }
      message( "Variable name was not found in data:" )
      print(vntest)
      stop()
    }
  }



  # on user scale
  ii = which(is.finite(P[["data"]][ , vnY ]))
  
  mq = quantile( P[["data"]][ ii, vnY ], probs=c(0.025, 0.5, 0.975) )

  O[["data_range"]] = c( 
    mean=mean(P[["data"]][ ii, vnY ]), 
    sd=sd(P[["data"]][ ii, vnY ]), 
    min=min(P[["data"]][ ii, vnY ]), 
    max=max(P[["data"]][ ii, vnY ]),  
    lb=mq[1], 
    median=mq[2], 
    ub=mq[3]  
  )  # on data /user scale not internal link

  # prefilter/transformation (e.g. translation to make all positive)
  if (exists("data_transformation", O)) P[["data"]][, vnY]  = O$data_transformation$forward( P[["data"]][, vnY] ) 

  # get hyper param scalings

  # temp Y var on link scale:
  yl = lnk_function( P[["data"]][, vnY ] )   # necessary in case of log(0)

  # offsets need to be close to 1 in user scale ( that is log(1)==0 in internal scale ) in experimental mode .. rescale  
  
  if ( !is.null(fm$offset_variable) )  {
    # link function applied to offsets here .. do not need to send log() 
    if (grepl("log[[:space:]]*[(]", vnO)) message("Probably do not want to transform the offset .. it is done internally in , unlike glm, inla, etc")

    obs = 1:nrow(P[["data"]])
    if (exists("tag", P[["data"]])) {
      obso = which(P[["data"]]$tag=="observations")
      if (length(obso) > 3) obs = obso
      obso = NULL
    }   
    
    P[["data"]][, vnO]  = lnk_function( P[["data"]][, vnO ])

    if (scale_offsets) {
      if (exists("offset_scale", O)) {
        O$offset_scale = median( P[["data"]][obs, vnO] , na.rm=TRUE )  # required to stabilize optimization
        P[["data"]][, vnO] = P[["data"]][, vnO] - O$offset_scale  # apply to all and overwrite, centering upon 0 (in user space 1)
        O$offset_scale_revert = function(x) { x - O$offset_scale }  # used for Y value and so opposite sign of "-" but applied to numerator is returns it to "-"
      }
    }    
    obs = NULL
    yl = yl - P[["data"]][, vnO]
  } 

  ll = which(is.finite(yl))
  O[["reference_sd"]] = sd(yl[ll] )

  H = hyperparameters(  O[["reference_sd"]], alpha=0.5, median(yl[ll], na.rm=TRUE) )  # sd slightly biased due to 0's being dropped .. but use of pc.priors that shrink to 0
  
  O$priors = H

  m = yl = ii = ll = fy = ol = NULL
  gc()


  if (redo_fit) {
    
    if ( !exists("inla.mode", P ) ) P[["inla.mode"]] = "experimental"
    if ( !exists("control.inla", P ) ) P[["control.inla"]] = list( strategy='adaptive' )
    if ( !exists("control.predictor", P ) ) P[["control.predictor"]] = list(compute=TRUE, link=1  ) #everything on link scale
    if ( !exists("control.compute", P ) ) P[["control.compute"]] = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE )


    if ( P[["inla.mode"]] != "experimental") {
      # location of this option has moved ... might move again
      # if ( !exists("control.results", P ) ) P[["control.results"]] = list(return.marginals.random=TRUE, return.marginals.predictor=TRUE )
      P[["control.compute"]][["return.marginals.predictor"]] = NULL
    }


    # control.fixed= list(mean.intercept=0, prec.intercept=0.001, mean=0, prec=0.001),
    # control.inla = list( strategy='adaptive', int.strategy='eb' )

    if (P[["verbose"]])  {
      message("Running model fit using the following data and options: \n")
      str(P)
    }

    fit  = NULL
    fit = try( do.call( inla, P ) )      

    if (inherits(fit, "try-error" )) {
      message("If you are using MSWindows and you get a popup complaining about 'inla stopped working',")
      message("try setting the flag in the following link to 1, using regedit. Be careful.")
      message("e.g., see: https://monitormyweb.com/guides/how-to-disable-stopped-working-message-in-windows")
      stop( "solution did not converge")
    }

    # to improve hyper param estimates..
    if (!is.null(improve.hyperparam.estimates)) if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.6, diff.logdens=9  )  # get improved estimates for the hyperparameters

    if (is.null(fit)) warning("model fit error")
    if ("try-error" %in% class(fit) ) warning("model fit error")

    if (!is.null(fn_res)) {
      # then save as separate files (fit, results)
      if (P[["verbose"]])  message( "Saving fit (this can be slow): ", fn_fit )
      save( fit, file=fn_fit, compress=compress )
    }

  }

  if (is.null(fit)) load( fn_fit )
  
  if (is.null(fit)) {
    message( "fit file not found: ", fn_fit )
    stop()
  }

  if (P[["verbose"]]) {
    print( summary(fit) )
    message( "   --- NOTE: parameter estimates are on link scale and not user scale")
  }
  
  # do the computations here as fit can be massive ... best not to copy, etc ..
  if (P[["verbose"]])  message( "\nComputing summaries and computing from posterior simulations (can be longer than model fitting depending upon no of posterior sims: 'nposteriors' ) ..." )

  truncate_upperbound = function( b, upper_limit, eps=1e-12 ) {
    # not used
    k = which( b[,1] > upper_limit )
    if (length(k) > 0) b[k,2] = 0
    return( b )
  }

  if (exists("data_transformation", O))  {
    backtransform = function( b ) {
      b[,1] =  O$data_transformation$backward( b[,1]   )
      return( b )
    }
  } 


  list_simplify = function(x) as.data.frame( t( as.data.frame( x )))
  exceedance_prob = function(x, threshold)  {1 - inla.pmarginal(q = threshold, x)}
  deceedance_prob = function(x, threshold)  { inla.pmarginal(q = threshold, x)}


  tokeep =  c("mean", "sd", "quant0.025", "quant0.5", "quant0.975")
 
  if (is.null(deceedance_threshold)) if (exists("deceedance_threshold", O)) deceedance_threshold = O[["deceedance_threshold"]]
  if (is.null(exceedance_threshold)) if (exists("exceedance_threshold", O)) exceedance_threshold = O[["exceedance_threshold"]]

  if (is.null(deceedance_threshold_predictions)) if (exists("deceedance_threshold_predictions", O)) deceedance_threshold_predictions = O[["deceedance_threshold_predictions"]]
  if (is.null(exceedance_threshold_predictions)) if (exists("exceedance_threshold_predictions", O)) exceedance_threshold_predictions = O[["exceedance_threshold_predictions"]]


  if ( "summary" %in% toget) {

    if (P[["verbose"]])  message("Extracting parameter summaries from marginals"  )
    if (!exists("summary", O)) O[["summary"]] = list()

    O[["summary"]][["direct"]] = summary(fit)

    # parameters
    # back-transform from marginals

    if (exists( "marginals.fixed", fit)) {
      V = fit$marginals.fixed  # make a copy to do transformations upon
      
      fi = which( grepl("Intercept", names(V) ))
      if (length(fi) > 0) {
        if (scale_offsets) {
          if ( exists("offset_scale_revert", O))  V[[fi]] = inla.tmarginal( O$offset_scale_revert, V[[fi]])  # on link scale
        }
      }

      V = lapply( V, function(x) inla.tmarginal( invlink, x)  )

      if (length(fi) > 0) {
        if ( exists("data_transformation", O))  V[[fi]] = inla.tmarginal( O$data_transformation$backward, V[[fi]]  ) # on user scale
      }

 
      W = NULL
      W = cbind ( t (apply_simplify( V, FUN=inla.zmarginal, silent=TRUE ) ) )  # 
      O[["summary"]][["fixed_effects"]] = W [, tokeep, drop =FALSE]
      W = NULL
      V = NULL

      if (P[["verbose"]]) {
        message( "")
        message( "Fixed effects")
        print(  O[["summary"]][["fixed_effects"]] )    
        message( "")
      } 

    }

    if (exists( "marginals.hyperpar", fit)) {
      
      # hyperpar (variance components)
      hyps =  rownames(fit$summary.hyperpar)

      prcs = grep( "^Precision.*", hyps, value=TRUE )
      if (length(prcs) > 0) {

        summary_inv_prec = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(pmax(y, 1e-12)), x) , silent=TRUE  )
        summary_inv_prec_1024 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=1024L) , silent=TRUE  )
        # summary_inv_prec_512 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=512L) , silent=TRUE  )

        precs = try( list_simplify( apply_simplify( fit$marginals.hyperpar[prcs], FUN=summary_inv_prec ) ), silent=TRUE )  # prone to integration errors ..
        if (any( inherits(precs, "try-error"))) precs = try( list_simplify( apply_simplify( fit$marginals.hyperpar[prcs], FUN=summary_inv_prec_1024 ) ), silent=TRUE )
        if (any( inherits(precs, "try-error")))  {
          if (P[["verbose"]])  {
            message( "NAN and/or Inf values encountered in marginals of paramater estimates.") 
            message( "Try an alternate parameterization as model may be over parameterized. ")
            message( "Copying fit summaries directly... ")
          }
          precs = fit$summary.hyperpar[prcs,1:5]
          precs[,c(1,3:5)] = 1/sqrt( precs[,c(1,3:5)] )
          rownames(precs) = gsub("Precision for", "SD", rownames(precs) )
          colnames(precs) = tokeep
          O[["summary"]][["random_effects"]] = precs

        } else {
          # precs[,"mode"] =  1/sqrt( fit$summary.hyperpar[prcs,"mode"]  )
          toadd = setdiff( colnames(O[["summary"]]), colnames(precs) )
          precs[,toadd] = NA
          rownames(precs) = gsub("Precision for", "SD", rownames(precs) )
          O[["summary"]][["random_effects"]] = precs[, tokeep, drop =FALSE]
        }
      }
      precs = NULL

      # update phi's, lambda's (used in besagproper2 -- Leroux model) .. etc
      params_to_match =  "^Rho.*|^GroupRho.*|^Phi.*|^Lambda.*|^Diagonal.*"
      known = grep( params_to_match, hyps, value=TRUE )
      unknown = setdiff( hyps, c(prcs, known) )
      if (length(unknown) > 0 ) {
        message( "Additional hyperparameters encountered. They will be treated as a normal hyperparameter (untransformed): ", unknown )
        k = unique( c(known, unknown) )
      } else {
        k = known
      }
      if (length(k) > 0) {
        ptm = try( list_simplify( apply_simplify( fit$marginals.hyperpar[k], FUN=function(x) inla.zmarginal( x, silent=TRUE  ) ) ), silent=TRUE )
        if (any( inherits(ptm, "try-error"))) {
          if (P[["verbose"]])  message( "Model may be over parameterized. NAN and Inf values encountered. Try alt parameterizations or smaller number of n or masking negative values")
        } else {
          #  alternatively: ptm[,"mode"] = apply_simplify( fit$marginals.hyperpar[k], FUN=function(x) inla.mmarginal( x ))
          O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], ptm[, tokeep, drop =FALSE] )
        }
      }
      k = params_to_match = known = unknown = ptm = NULL
      gc()

      if (P[["verbose"]])  {
        message( "")
        message( "Random effects:")
        print(  O[["summary"]][["random_effects"]] )   
        message( "   --- NOTE: 'SD *' are on link scale and not user scale")
        message( "")
      }

    }
    
  }


  if (any( grepl("random", toget) )) {
    if (!exists("random", O)) O[["random"]] = list()
    if (exists("marginals.random", fit)) {

      if ("random_other" %in% toget) {
        summary_inv_random = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )
        raneff = setdiff( names( fit$marginals.random ), c(vnS, vnST, vnSI, vnSTI ) )
        for (re in raneff) {
          if (P[["verbose"]])  message("Extracting random covariates from marginals:  ", re  )
          g = fit$marginals.random[[re]]
          O[["random"]] [[re]] = list_simplify ( apply_simplify( g, summary_inv_random ) )  [, tokeep, drop =FALSE]
          O[["random"]] [[re]]$ID = fit$summary.random[[re]]$ID
        }
        g = raneff = NULL
      }
      gc()

      if ("random_spatial" %in% toget) {
        # space only
        summary_inv_random = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )


        bym = iid = NULL
        matchto = list( space=O[[vnS]] )
        matchfrom0 = NULL  # used for matching spatial effects in exceed/deceed .. required as bym2 uses 2x nsp for by and iid

        if ( exists(vnSI, fit$marginals.random)  | exists(vnS, fit$marginals.random) ) {
    
          if (P[["verbose"]])  message("Extracting random spatial errors from marginals"  )

          if ( exists(vnSI, fit$marginals.random) ) {
            O[["random"]] [[vnSI]] = list()  # space as a main effect
            model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnSI) ]  # should be iid
            m = list_simplify ( apply_simplify( fit$marginals.random[[vnSI]], summary_inv_random ) )

            # single spatial effect (eg in conjucyion with besag) .. indexing not needed but here in case more complex models ..
            Z = expand.grid( space=O[[vnS]], type=model_name, stringsAsFactors =FALSE )

            iid =  which(Z$type==model_name)
            matchfrom0 = matchfrom = list( space=Z[["space"]][iid] )
            W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..

            for (k in 1:length(names(m))) {
              W[,k] = reformat_to_array( input = unlist(m[iid,k]), matchfrom=matchfrom, matchto=matchto )
            }
            O[["random"]] [[vnSI]] [[model_name]] = W [, tokeep, drop =FALSE]
            W = NULL
          }

          if ( exists(vnS, fit$marginals.random) ) {
            O[["random"]] [[vnS]] = list()  # space as a main effect
            model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnS) ]
            m = list_simplify ( apply_simplify( fit$marginals.random[[vnS]], summary_inv_random ) )

            if ( model_name %in% c("bym", "bym2") ) {
              # bym2 effect: bym and iid simultaneously
              Z = expand.grid( space=O[[vnS]], type = c("iid", model_name), stringsAsFactors =FALSE )

              #  extract iid main effects
              iid = which(Z$type=="iid")
              matchfrom0 = matchfrom = list( space=Z[["space"]][iid] )
              W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) )
              names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
              for (k in 1:length(names(m))) {
                W[,k] = reformat_to_array( input = unlist(m[iid,k]), matchfrom=matchfrom, matchto=matchto )
              }
              O[["random"]] [[vnS]] [["iid"]] = W [, tokeep, drop =FALSE]
              W = NULL

            } else {
              # single spatial effect (eg besag, etc)
              Z = expand.grid( space=O[[vnS]], type=model_name, stringsAsFactors =FALSE )
            }

            bym = which(Z$type==model_name)
            matchfrom = list( space=Z[["space"]][bym] )
            if (!is.null(matchfrom0)) matchfrom0 = matchfrom # ie besag etc

            #  extract spatial effect ... besag, etc main effects or bym2 
            W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            for (k in 1:length(names(m))) {
              W[,k] = reformat_to_array( input = unlist(m[bym,k]), matchfrom=matchfrom, matchto=matchto )
            }
            O[["random"]] [[vnS]] [[model_name]] = W [, tokeep, drop =FALSE]
            m = NULL
          }


          # POSTERIOR SIMS
          if ( exists(vnSI, fit$marginals.random ) &  exists(vnS, fit$marginals.random ) )  {
            # sep besag + iid  
            selection=list()
            selection[vnSI] = 0  # 0 means everything matching space
            selection[vnS] = 0  # 0 means everything matching space
            aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE , num.threads=num.threads )  # 0 means everything matching space
            # order of effects gets messed up .. must use names
            aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
            aa_iid = which(aa_rn==vnSI)
            aa_bym = which(aa_rn==vnS)
            g = apply_simplify( aa, function(x) {invlink(x$latent[aa_iid] + x$latent[aa_bym] ) } )
          } else if ( exists(vnSI, fit$marginals.random ) & ! exists(vnS, fit$marginals.random ) )  {
            # iid  only
            selection=list()
            selection[vnSI] = 0  # 0 means everything matching space
            aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE , num.threads=num.threads )  # 0 means everything matching space
            g = apply_simplify( aa, function(x) {invlink(x$latent[iid] ) } )
          } else if ( ! exists(vnSI, fit$marginals.random ) & exists(vnS, fit$marginals.random ) ) {
            # besag only or bym/bym2
            selection=list()
            selection[vnS] = 0  # 0 means everything matching space
            aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE , num.threads=num.threads )  # 0 means everything matching space
            if ( model_name %in% c("bym", "bym2") ) {
              g = apply_simplify( aa, function(x) {invlink(x$latent[iid] + x$latent[bym] ) } )
            } else {
              g = apply_simplify( aa, function(x) {invlink(x$latent[bym] ) } )
            }
          }
          aa = NULL

          mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
          mm = apply( g, 1, mean, na.rm =TRUE)
          ms = apply( g, 1, sd, na.rm =TRUE)
          W = cbind(mm, ms, mq)
          mq = mm = ms = NULL

          attr(W, "dimnames") = list( space=O[[vnS]], stat=tokeep  )
          O[["random"]] [[vnS]] [["combined"]] = W
          W = NULL

          if (!is.null(exceedance_threshold)) {
            if (P[["verbose"]])  message("Extracting random spatial errors exceedence"  )
            for ( b in 1:length(exceedance_threshold)) {
              m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold[b]) ) ) / nposteriors
              W = reformat_to_array( input = m, matchfrom=matchfrom0, matchto = matchto )
              names(dimnames(W))[1] = vnS
              dimnames( W )[[vnS]] = O[[vnS]]
              O[["random"]] [[vnS]] [["exceedance"]] [[as.character(exceedance_threshold[b])]] = W
              W = m = NULL
            }
          }

          if (!is.null(deceedance_threshold)) {
            if (P[["verbose"]])  message("Extracting random spatial errors deceedance"  )
            for ( b in 1:length(deceedance_threshold)) {
              m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold[b]) ) ) / nposteriors
              W = reformat_to_array( input = m, matchfrom=matchfrom0, matchto=matchto  )
              names(dimnames(W))[1] = vnS
              dimnames( W )[[vnS]] = O[[vnS]]
              O[["random"]] [[vnS]] [["deceedance"]] [[as.character(deceedance_threshold[b])]] = W
              W = m = NULL
            }
          }
        }
      }

      Z = g = NULL
      gc()


      if ("random_spatiotemporal"  %in% toget ) {
        # space-time
        summary_inv_random = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )

        g = NULL
        
        bym = iid = NULL
        matchto = list( space=O[[vnS]], time=O[[vnT]]  )
        matchfrom0 = NULL

        if (exists(vnST, fit$marginals.random ) | exists(vnSTI, fit$marginals.random ) ) {
          if (P[["verbose"]])  message("Extracting random spatiotemporal errors from marginals"  )

          if (exists(vnSTI, fit$marginals.random )) {
            O[["random"]] [[vnSTI]] = list()
            model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnSTI) ]  # should be iid
            m = list_simplify ( apply_simplify( fit$marginals.random[[vnSTI]], summary_inv_random ) )

            Z = expand.grid( space=O[[vnS]], type=model_name, time=O[[vnT]], stringsAsFactors =FALSE )

            #  spatiotemporal interaction effects  iid
            iid =  which(Z$type==model_name)
            matchfrom0 = matchfrom = list( space=Z[["space"]][iid], time=Z[["time"]][iid]  )
            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = unlist(m[iid,k]), matchfrom=matchfrom, matchto=matchto )
            }
            O[["random"]] [[vnSTI]] [[model_name]] = W [,, tokeep, drop =FALSE]
            W = NULL
          }

          if (exists(vnST, fit$marginals.random )) {

            O[["random"]] [[vnST]] = list()
            model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnST) ]
            m = list_simplify ( apply_simplify( fit$marginals.random[[vnST]], summary_inv_random ) )

            if ( model_name %in% c("bym", "bym2") ) {
              # bym2 effect: bym and iid with annual results
              Z = expand.grid( space=O[[vnS]], type = c("iid", model_name), time=O[[vnT]], stringsAsFactors =FALSE )
 
              #  spatiotemporal interaction effects  iid
              iid = which(Z$type=="iid")
              matchfrom0 = matchfrom = list( space=Z[["space"]][iid], time=Z[["time"]][iid]  )
              W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
              names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
              names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
              for (k in 1:length(names(m))) {
                W[,,k] = reformat_to_array(  input = unlist(m[iid,k]), matchfrom = matchfrom, matchto = matchto )
              }
              O[["random"]] [[vnST]] [["iid"]] = W [,, tokeep, drop =FALSE]
              W = NULL
            } else {
              # besag effect: with annual results
              Z = expand.grid( space=O[[vnS]], type =model_name, time=O[[vnT]], stringsAsFactors =FALSE )
            }

            #  spatiotemporal interaction effects  bym
            bym =  which(Z$type==model_name)
            matchfrom = list( space=Z[["space"]][bym], time=Z[["time"]][bym]  )
            if (!is.null(matchfrom0)) matchfrom0 = matchfrom # ie besag etc

            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
            matchfrom = list( space=Z[["space"]][bym],  time=Z[["time"]][bym]  )
            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = unlist(m[bym,k]), matchfrom = matchfrom, matchto = matchto  )
            }
            O[["random"]] [[vnST]] [[model_name]] = W [,, tokeep, drop =FALSE]
            W = NULL
            m = NULL


            # combined effects from posterior simulations
            if ( exists(vnSTI, fit$marginals.random ) &  exists(vnST, fit$marginals.random ) ) {
              # sep besag +iid  
                selection=list()
                selection[vnST] = 0  # 0 means everything matching space_time
                selection[vnSTI] = 0  # 0 means everything matching space_time iid
                aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE, num.threads=num.threads  )  # 0 means everything matching space
                # order of effects gets messed up .. must use names
                aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
                aa_iid = which(aa_rn==vnSTI)
                aa_bym = which(aa_rn==vnST)
                g = apply_simplify( aa, function(x) {invlink(x$latent[aa_iid] + x$latent[aa_bym] ) } )
 
            } else if (exists(vnSTI, fit$marginals.random ) & ! exists(vnST, fit$marginals.random )) {
              # iid  only
                selection=list()
                selection[vnSTI] = 0  # 0 means everything matching space_time
                aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE, num.threads=num.threads )  # 0 means everything matching space
                g = apply_simplify( aa, function(x) {invlink(x$latent[iid] ) } )
 
            } else if ( !exists(vnSTI, fit$marginals.random ) & exists(vnST, fit$marginals.random ) ) {
              # besag  only or bym/bym2
                selection=list()
                selection[vnST] = 0  # 0 means everything matching space_time
                aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE, num.threads=num.threads  )  # 0 means everything matching space
                if ( model_name %in% c("bym", "bym2") ) {
                  g = apply_simplify( aa, function(x) {invlink(x$latent[iid] + x$latent[bym] ) } )
                } else {
                  g = apply_simplify( aa, function(x) {invlink(x$latent[bym] ) } )
                }
            } 

            aa = NULL

            mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
            mm = apply( g, 1, mean, na.rm =TRUE)
            ms = apply( g, 1, sd, na.rm =TRUE)

            m = data.frame( cbind(mm, ms, mq) )
            names(m) = tokeep
            mm = ms = mq = NULL
            
      
            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
            
            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = m[,k], matchfrom=matchfrom0, matchto=matchto )
            }
            O[["random"]] [[vnST]] [["combined"]] = W [,, tokeep, drop =FALSE]
            W = NULL
          }
        }

        if (!is.null(g)) { 
          if (!is.null(exceedance_threshold)) {
            if (P[["verbose"]])  message("Extracting random spatiotemporal errors exceedence"  )
            for ( b in 1:length(exceedance_threshold)) {
              m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold[b] ) ) ) / nposteriors
              W = reformat_to_array( input=m, matchfrom=matchfrom0,  matchto=matchto )
              names(dimnames(W))[1] = vnS
              dimnames( W )[[vnS]] = O[[vnS]]
              m = NULL
              if (O[["dimensionality"]] == "space-time"  ) {
                names(dimnames(W))[2] = vnT
                dimnames( W )[[vnT]] = O[[vnT]]
              }
              if (O[["dimensionality"]] == "space-time-cyclic" ) {
                names(dimnames(W))[3] = vnU
                dimnames( W )[[vnU]] = O[[vnU]]
              }
              O[["random"]] [[vnST]] [["exceedance"]] [[as.character(exceedance_threshold[b])]] = W
              W = NULL
            }
          }

          if (!is.null(deceedance_threshold)) {
            if (P[["verbose"]])  message("Extracting random spatiotemporal errors deceedance"  )

            for ( b in 1:length(deceedance_threshold)) {
              m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
              W = reformat_to_array( input = m, matchfrom=matchfrom0,  matchto=matchto )
              names(dimnames(W))[1] = vnS
              dimnames( W )[[vnS]] = O[[vnS]]
              m = NULL
              if (O[["dimensionality"]] == "space-time"  ) {
                names(dimnames(W))[2] = vnT
                dimnames( W )[[vnT]] = O[[vnT]]
              }
              if (O[["dimensionality"]] == "space-time-cyclic" ) {
                names(dimnames(W))[3] = vnU
                dimnames( W )[[vnU]] = O[[vnU]]
              }
              O[["random"]] [[vnST]] [["deceedance"]] [[as.character(deceedance_threshold[b])]] = W
              W = NULL
            }
          }  
          # end space-time
        }
        Z = g = NULL
        gc()
      }
    }  # end random effects
  }


  if ("predictions"  %in% toget ) {

    if (!exists("predictions", O)) O[["predictions"]] = list()

    summary_inv_predictions = function(x) inla.zmarginal( x, silent=TRUE  )
    
    if (!exists("tag", P[["data"]])) P[["data"]]$tag="predictions" # force predictions for all data

    # adjusted by offset
    if (exists("marginals.fitted.values", fit)) {

      if (P[["verbose"]])  message("Extracting predictions from marginals"  )

      if (  O[["dimensionality"]] == "space" ) {
        ipred = which( P[["data"]]$tag=="predictions"  &  P[["data"]][,vnS0] %in% O[[vnS]] )  # filter by S and T in case additional data in other areas and times are used in the input data
        g = fit$marginals.fitted.values[ipred]
        if (scale_offsets) {
          if ( exists("offset_scale_revert", O) ) g = apply_generic( g, function(u) {inla.tmarginal( O$offset_scale_revert, u) } )
        }    
        g = apply_generic( g, function(u) {inla.tmarginal( invlink, u) } )    
        if ( exists("data_transformation", O))  g = apply_generic( g, backtransform )

        m = list_simplify ( apply_simplify( g, summary_inv_predictions ) )
        W = array( NA, dim=c( length(O[[vnS]]),  length(names(m)) ),  dimnames=list( space=O[[vnS]], stat=names(m) ) )
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        
        matchfrom = list( space=P[["data"]][ ipred, vnS0] ) 
        matchto = list( space=O[[vnS]] )

        for (k in 1:length(names(m))) {
          W[,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=matchfrom, matchto=matchto )
        }
        O[["predictions"]] = W[, tokeep, drop =FALSE]
        W = m = NULL
      }


      if (O[["dimensionality"]] == "space-time"  ) {
        ipred = which( P[["data"]]$tag=="predictions" & P[["data"]][,vnS0] %in% O[[vnS]] & P[["data"]][,vnT0] %in% O[[vnT]] )
        g = fit$marginals.fitted.values[ipred]   
        if (scale_offsets) {
          if ( exists("offset_scale_revert", O) ) g = apply_generic( g, function(u) {inla.tmarginal( O$offset_scale_revert, u) } )
        }

        g = apply_generic( g, function(u) {inla.tmarginal( invlink, u) } )    
        if (exists("data_transformation", O)) g = apply_generic( g, backtransform )

        m = list_simplify ( apply_simplify( g, summary_inv_predictions ) )
        W = array( NA, dim=c( length(O[[vnS]]), length(O[[vnT]]), length(names(m)) ),  dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
      
        matchfrom = list( space=P[["data"]][ipred,vnS0], time=P[["data"]][ipred,vnT0] )
        matchto = list( space=O[[vnS]], time=O[[vnT]] )

        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=matchfrom, matchto=matchto)
        }
        O[["predictions"]] = W[,, tokeep, drop =FALSE]
        W = m = NULL
      }


      if ( O[["dimensionality"]] == "space-time-cyclic" ) {
        ipred = which( P[["data"]]$tag=="predictions" & P[["data"]][,vnS0] %in% O[[vnS]]  &  P[["data"]][,vnT0] %in% O[[vnT]] &  P[["data"]][,vnU0] %in% O[[vnU]])  # ignoring U == predict at all seassonal components ..
        g = fit$marginals.fitted.values[ipred]   
        if (scale_offsets) {
          if ( exists("offset_scale_revert", O) ) g = apply_generic( g, function(u) {inla.tmarginal( O$offset_scale_revert, u) } )
        }    
        g = apply_generic( g, function(u) {inla.tmarginal( invlink, u) } )    
        if (exists("data_transformation", O)) g = apply_generic( g, backtransform )

        m = list_simplify ( apply_simplify( g, summary_inv_predictions ) )
        W = array( NA, dim=c( length(O[[vnS]]), length(O[[vnT]]), length(O[[vnU]]), length(names(m)) ),  dimnames=list( space=O[[vnS]], time=O[[vnT]], cyclic=O[[vnU]], stat=names(m) ) )
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
        names(dimnames(W))[3] = vnU  # need to do this in a separate step ..

        matchfrom = list( space=P[["data"]][ipred, vnS0], time=P[["data"]][ipred, vnT0], cyclic=P[["data"]][ipred, vnU0] )
        matchto = list( space=O[[vnS]], time=O[[vnT]], cyclic=O[[vnU]] )

        for (k in 1:length(names(m))) {
          W[,,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=matchfrom, matchto=matchto )
        }
        O[["predictions"]] = W[,,, tokeep, drop =FALSE]
        W = m = NULL
      }

      if (!is.null(exceedance_threshold_predictions)) {
        for (b in exceedance_threshold_predictions) {  
          m = list_simplify ( apply_simplify( g, FUN=exceedance_prob, threshold=exceedance_threshold_predictions[b] ) )
          W = reformat_to_array(  input = unlist(m ), matchfrom = matchfrom, matchto = matchto )
          names(dimnames(W))[1] = vnS
          dimnames( W )[[vnS]] = O[[vnS]]
          if (O[["dimensionality"]] == "space-time"  ) {
            names(dimnames(W))[2] = vnT
            dimnames( W )[[vnT]] = O[[vnT]]
          }
          if (O[["dimensionality"]]=="space-time-cyclic") {
            names(dimnames(W))[3] = vnU
            dimnames( W )[[vnU]] = O[[vnU]]
          }
          O[["predictions"]] [["exceedance"]] [[ as.character(exceedance_threshold[b]) ]]= W
          W = m = NULL
        }
      }

      if (!is.null(deceedance_threshold_predictions)) {
        for (b in deceedance_threshold_predictions) {  
          m = list_simplify ( apply_simplify( g, FUN=deceedance_prob, threshold=deceedance_threshold_predictions[b] ) )
          W = reformat_to_array(  input = unlist(m ), matchfrom = matchfrom, matchto = matchto )
          names(dimnames(W))[1] = vnS
          dimnames( W )[[vnS]] = O[[vnS]]
          if (O[["dimensionality"]] == "space-time"  ) {
            names(dimnames(W))[2] = vnT
            dimnames( W )[[vnT]] = O[[vnT]]
          }
          if (O[["dimensionality"]]=="space-time-cyclic") {
            names(dimnames(W))[3] = vnU
            dimnames( W )[[vnU]] = O[[vnU]]
          }
          O[["predictions"]] [["deceedance"]] [[ as.character(deceedance_threshold[b]) ]]= W
          W = m = NULL
        }
      }
    }


  }

  # copy the modified data in case needed for plotting ..
  if (!is.null(P[["data"]])) {
    O[["data"]] = P[["data"]]
    P[["data"]] = NULL
    gc()
  }

  if (!is.null(sppoly)) O[["sppoly"]] = sppoly
  if (!is.null(space.id)) O[["space.id"]] = space.id

  if (!is.null(fn_res)) {
    # then save as separate files (fit, results)
    save( O, file=fn_res, compress=compress )
    if (P[["verbose"]])  message( "Summary saved as: ", fn_res )
    fit$results = O

  } else {
    # save as a single file
    fit$results = O
    save( fit, file=fn_fit, compress=compress )
    if (P[["verbose"]])  message( "fit and summary saved as: ", fn_fit )
  }
  return(fit)

}
