
carstm_model_inla = function(
  O, # p parameter list
  DS = "",
  sppoly =NULL,
  fit = NULL,
  space_id=NULL, time_id=NULL, cyclic_id=NULL, 
  fn_fit=tempfile(pattern="fit_", fileext=".rdata"), 
  fn_res=NULL, 
  redo_fit = TRUE,
  theta=NULL,
  compress="gzip",
  compression_level=1,
  toget = c("summary", "random_spatial", "random_spatiotemporal", "predictions"), 
  nposteriors=NULL, 
  posterior_simulations_to_retain=NULL,
  exceedance_threshold=NULL, 
  deceedance_threshold=NULL, 
  exceedance_threshold_predictions=NULL,
  deceedance_threshold_predictions=NULL,
  debug=FALSE,
  eps = 1e-32,
  ... ) {
    
  if (DS=="modelled_fit") {
    if (!is.null(fn_fit)) {
      if (file.exists(fn_fit)) {
        if (grepl("\\.RDS$", fn_fit)) {
          fit = readRDS(fn_fit)
        } else {
          load( fn_fit )
        }
      }
      if (is.null(fit)) message("modelled fit not found.")
    }
    return( fit )
  }

  if (DS=="modelled_summary") {  # model.*modelled
    if (!is.null(fn_res)) {
      O = NULL
      if (file.exists(fn_res)) {
        if (grepl("\\.RDS$", fn_fit)) {
          res = readRDS(fn_res)
        } else {
          load( fn_res)
        }
      }
      if (is.null(O)) message(" summary not found.")
      return( O )
    } else {
      if (file.exists(fn_fit)){
        if (grepl("\\.RDS$", fn_fit)) {
          fit = readRDS(fn_fit)
        } else {
          load( fn_fit )
        }
      } 
      if (!is.null(fit)) {
        if (exists( "results", fit)) {
          return( fit$results )
        }
      }
      message("modelled results not found. .. try to run extraction: '' ")
    }
  }
  
  if (exists("debug")) if (is.logical(debug)) if (debug) browser()

  run_start  = Sys.time()

 
  ### 1. Prepare inla_args and vars used in both modelling and extraction
  
  inla_args = list(...)  # INLA options to be passed directly to it  
  
  # some checks
  if (!exists("verbose", inla_args)) inla_args[["verbose"]]=FALSE
  be_verbose = inla_args[["verbose"]]


  if (exists("num.threads", inla_args)) {
    num.threads = inla_args[["num.threads"]]
  } else {
    num.threads = "1:1"
  }

  inla.setOption(num.threads=num.threads)
  
  # n cores to use for posterior marginals in mcapply .. can be memory intensive so make it a bit less than "num.threads" ..
  if (exists("mc.cores", inla_args)) {
    mc.cores = inla_args[["mc.cores"]]
  } else if ( !is.null(num.threads) ) {
    mc.cores = as.numeric( unlist(strsplit(num.threads, ":") )[1] )
  } else {
    mc.cores = 1
  }

    # Formula parsing
    if ( !exists("formula", inla_args ) ) {
      if (exists("formula", O ) ) {
        inla_args[["formula"]] = O[["formula"]]
        if ( be_verbose ) message( "Formula found in options, O" )
      } else {
        stop( "Model formula was not provided")
      }
    }

    O[["fm"]] = parse_formula( inla_args[["formula"]] )
    
    if ( be_verbose ) {
      message( "Formula parsed as follows. Check in case of errors:" )
      print(O[["fm"]])
    }


    # labels .. not all used but defining the here makes things simpler below
    # shorten var names
    re = O[["fm"]]$random_effects

    vO = O[["fm"]]$offset_variable
    vY = O[["fm"]]$dependent_variable 

    vS = O[["fm"]]$vn$S
    vT = O[["fm"]]$vn$T
    vU = O[["fm"]]$vn$U
    vS2 = O[["fm"]]$vn$S2 
    vT2 = O[["fm"]]$vn$T2
    vU2 = O[["fm"]]$vn$U2
    vS3 = O[["fm"]]$vn$S3 

    # family related  
    if ( !exists("family", inla_args ) ) {
      if ( exists("family", O) ) {
        if ( be_verbose ) message( "Family found in options, O" )
        inla_args[["family"]] = O[["family"]]   
      } else {
        inla_args[["family"]] = "gaussian"
      }
    }

    if (!exists("family", inla_args))  stop( "Model family was not provided")
    O[["family"]] = inla_args[["family"]]  # copy to O in case it was provided as inla_args
 
    if ( inla_args[["family"]] == "gaussian" ) {
      lnk_function = inla.link.identity
    } else if ( inla_args[["family"]] == "lognormal" ) {
      lnk_function = inla.link.log
    } else if ( grepl( ".*poisson", inla_args[["family"]])) {
      lnk_function = inla.link.log
    } else if ( grepl( ".*nbinomial", inla_args[["family"]])) {
      lnk_function = inla.link.log
      invlink_id  =  "exp"
    } else if ( grepl( ".*binomial", inla_args[["family"]])) {
      lnk_function = inla.link.logit
    } 

    O[["invlink"]] = invlink = function(x) lnk_function( x,  inverse=TRUE )


    # changes in INLA data structures noted in 2023:
    # posterior sims is on link scale ((for both experimental and classical))
    # marginals.fitted.values are on response scale (for both experimental and classical)
    # summary.fitted.values on response scale ((for both experimental and classical))


  ### START Modelling
 
  if (redo_fit) {

    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (!exists("data", inla_args)) {
      if (exists("data", O)) {
        if ( be_verbose) message( "Data not passed as an argument, using data found in options, O")
        inla_args[["data"]] = O[["data"]]  # priority to inla_args
        O[["data"]] = NULL  # reduce mem usage
      } else {
        stop("No data")
      }
    }
    
    if (any(class(inla_args[["data"]])=="character")) {
      if ( be_verbose) message( "Data is a function. Running ...")
      inla_args[["data"]] = try( eval(parse(text=inla_args[["data"]]) ) )
    }

    if (inherits(inla_args[["data"]], "try-error"))  inla_args[["data"]] = NULL

    if (is.null(inla_args[["data"]])) stop("Data not found")

    setDT(inla_args[["data"]])
 
    if ( !exists("carstm_model_label", O) ) {
      if ( exists("model_label", O) ) {
        O[["carstm_model_label"]] = O[["model_label"]] 
      } else {
        stop("carstm_model_label should be specified")
      }
    }

    DM = unique( re$dimensionality )  
    dims = c("s", "t", "c")
    RR = setdiff( dims, setdiff( dims, DM ) )
    SS = paste0(
      c("space", "time", "cyclic")[match(RR, dims)], 
      collapse="-"
    ) 

    if ( exists("dimensionality", O) ) {
      if (O[["dimensionality"]] != SS ) {
        message("Dimensionality parameters specified (left) do not match those in formula (right): " )
        message( O[["dimensionality"]], " vs ", SS )
        message( "This is OK when predicting to one time slice, but verify that this is your intent ...")
      }
    } else  {
      O[["dimensionality"]] = SS
      message("Check dimensionality guessed from formula:  ", O[["dimensionality"]]) 
    }

    if ( !is.null(O[["fm"]][["vn"]][["S"]]) ) {
      # this is about model data inputs not dimensionality of outputs

      if (is.null(sppoly)) if (exists("sppoly", O)) sppoly = O$sppoly
      if (is.null(sppoly)) if (exists("areal_units")) {
        sppoly = areal_units( O )   
        if (!exists(sppoly)) sppoly = NULL
      }
      if (is.null(sppoly)) stop( "sppoly is required") 
      
      # test to see if sppoly has been altered:
      sp_nb_auid = NULL
      sp_nb_auid = try( attributes(attributes(sppoly)$NB_graph)$region.id, silent=TRUE)  # order of neighbourhood graph 
      if (!inherits("try-error", sp_nb_auid)) {
        o = match( sppoly$AUID, sp_nb_auid)  
        if (length(unique(diff(o))) !=1) {
          stop( "Neighbourhood indices and polygons AUID order are mismatched .. this will not end well" )
        }
        o = NULL
        O[["space_id"]] = sp_nb_auid
        sp_nb_auid = NULL
      }

      O[["sppoly"]] = sppoly  # copy in case mapping directly from O
      
      # the master / correct sequence of the AU's and neighbourhood matrix index values
      if (!is.null(space_id)) {
        O[["space_id"]] = space_id
        space_id = NULL
      } else {
        if (exists("space_id", attributes(sppoly)) ) {
          O[["space_id"]] = as.character( attributes(sppoly)$space_id )
        } else if (exists("region.id", attributes(sppoly)) ) {
          O[["space_id"]] = as.character( attributes(sppoly)$region.id )
        } else if (exists("AUID", sppoly) ) {
          O[["space_id"]] = as.character( sppoly[["AUID"]] )
        }
      }

      if (!exists("space_id", O)) stop( "space_id could not be determined from data")
      if (!exists("space_id", attributes(sppoly)) ) attributes(sppoly)$space_id = O[["space_id"]]  # copy as attribute in case
      
      O[["space_n"]] = length( O[[ "space_id" ]] )

      # better formatted areal unit names (for reporting or plotting) (outside of carstm))
      if (!exists("space_name", O)) {
        if (exists("space_name", attributes(sppoly))) {
          O[["space_name"]] = attributes(sppoly)$space_name  
        } else {
          O[["space_name"]] = O[["space_id"]]
        }
      }
      if (!exists("space_name", attributes(sppoly)) ) attributes(sppoly)$space_name = O[["space_name"]]  # copy as attribute in case

      if (exists("AUID_label", sppoly)) {
        # long name useful in some cases (for plotting)
          O[["space_label"]] = sppoly[["AUID_label"]]
          if (!exists("space_label", attributes(sppoly)) ) attributes(sppoly)$space_label = O[["space_label"]]  # copy as attribute in case
      }

      if ( is.null(O[["fm"]][["vn"]][["T"]]) ) {
        # force to carry a "time" to make dimensions of predictions simpler to manipulate 
        inla_args[["data"]][["time"]] = -1  
      }
  
      missingS = which(is.na(inla_args[["data"]][[vS]] ))
      if (length( missingS ) > 0 ) {
        warning( "Data areal units and space_id (from sppoly) do not match ... this should not happen:")
        print( head(missingS) )
        warning( "Dropping them from analysis and assuming neighbourhood relations are OK even though this is unlikely!")
        inla_args[["data"]] =  inla_args[["data"]] [ -missingS , ]
      }
      missingS = NULL
    
      missingS = unique( setdiff( O[["space_id"]], unique(  inla_args[["data"]][[vS]] ) ) )
      if (length(missingS) > 0) {
        warning( "No. of areal unique units in data do not match those in sppoly:", paste0( head(missingS), sep=", "))
      }
      missingS = NULL

    }

    if ( !is.null(O[["fm"]][["vn"]][["T"]])  | !is.null(O[["fm"]][["vn"]][["U"]])  ) {
      # carstm internal ID
      if (!is.null(time_id)) {
        O[["time_id"]] = time_id
        time_id = NULL
      } else {
        if (exists("yrs", O)) {
          if (is.factor(O[["yrs"]])) {
            O[["time_id"]] = as.numeric(O[["yrs"]])
          } else {
            O[["time_id"]] = as.numeric(O[["yrs"]])  # redundant but make explicit
          }
        }
      }
      if (!exists("time_id", O)) stop( "time_id or yrs need to be provided")
      
      O[["time_n"]] = length( O[[ "time_id" ]] )

      # for plot labels, etc .. time (year) gets swapped out for time index (outside of carstm)
      if (!exists("time_name", O)) {
        if (exists("yrs", O)) {
          O[["time_name"]] = as.character(O[[ "yrs" ]] )
        } else {
          O[["time_name"]] = as.character(O[[ "time_id" ]] ) 
        }  
      }
  
      # sub-annual time 
      if (!is.null(O[["fm"]][["vn"]][["U"]]) )  {
        # carstm internal ID
        if (!is.null(cyclic_id)) {
          O[["cyclic_id"]] = cyclic_id
          cyclic_id = NULL
        } else {
          if (exists("cyclic_levels", O)) {
            if (is.factor(O[["cyclic_levels"]])) {
              O[["cyclic_id"]] = as.numeric(O[["cyclic_levels"]])
            } else {
              O[["cyclic_id"]] = O[["cyclic_levels"]]
            }
          }
        }
        if (!exists("cyclic_id", O)) stop( "cyclic_id or cyclic_levels need to be provided")
        
        O[["cyclic_n"]] = length( O[[ "cyclic_id" ]] )

        # for plot labels, etc .. season (time of year) gets swapped out for cyclic index   (outside of carstm)
        if (!exists("cyclic_name", O)) {
          if (exists("cyclic_levels", O)) {
            O[["cyclic_name"]] = as.character(O[[ "cyclic_levels" ]] )
          } else {
            O[["cyclic_name"]] = as.character(O[[ "cyclic_id" ]] ) 
          }  
        }
      }
    } 

    ii = which(is.finite(inla_args[["data"]][[vY]]))

    if ( be_verbose ) {
      # dev.new()
      # hist( inla_args[["data"]][[vY]][ii] , main="Histogram of input variable to model" )
    }

    mq = quantile( inla_args[["data"]][[vY]][ii] , probs=c(0.025, 0.5, 0.975) )

    O[["data_range"]] = c( 
      mean=mean(inla_args[["data"]][[vY]][ii] ), 
      sd=sd(inla_args[["data"]][[vY]][ii] ), 
      min=min(inla_args[["data"]][[vY]][ii] ), 
      max=max(inla_args[["data"]][[vY]][ii] ),  
      lb=mq[1], 
      median=mq[2], 
      ub=mq[3]  
    )  # on data /user scale not internal link
    mq = NULL

    # prefilter/transformation (e.g. translation to make all positive)
    if (exists("data_transformation", O)) inla_args[["data"]][[vY]]  = O$data_transformation$forward( inla_args[["data"]][[vY]] ) 

    # get hyper param scalings

    # temp Y var on link scale:
    # on user scale
    yl = inla_args[["data"]][[vY]]

    if ( grepl( ".*binomial", inla_args[["family"]])) {
      # for binomial, prob=0,1, becomes infinite so minor fix for hyper parameter prior approximation
      tweak = 0.05 # tail truncation prob
      if (exists("habitat.threshold.quantile", O)) tweak = O[["habitat.threshold.quantile"]]
      yl [ yl==1 ] = 1 - tweak
      yl [ yl==0 ] = tweak
    }
    

    yl = lnk_function( yl )   # necessary in case of log(0)
    
    # offsets need to be close to 1 in user scale ( that is log(1)==0 in internal scale ) in experimental mode .. rescale  
    
    if ( !is.null(vO) )  {
      # link function applied to offsets here .. do not need to send log() 
      if (grepl("log[[:space:]]*[(]", vO)) message("Probably do not want to transform the offset .. it is done internally in , unlike glm, inla, etc")

      obs = 1:nrow(inla_args[["data"]])
      if (exists("tag", inla_args[["data"]])) {
        obso = which(inla_args[["data"]][["tag"]]=="observations")
        if (length(obso) > 3) obs = obso
        obso = NULL
      }   
      
      inla_args[["data"]][[vO]]  = lnk_function( inla_args[["data"]][[vO]])
  
      obs = NULL
      yl = yl - inla_args[["data"]][[vO]]

    } 

    ll = which(is.finite(yl))

    mqi = quantile( yl[ll], probs=c(0.025, 0.5, 0.975) )

    O[["data_range_internal"]] = c( 
      mean=mean(yl[ll]), 
      sd=sd(yl[ll]), 
      min=min(yl[ll]), 
      max=max(yl[ll]),  
      lb=mqi[1], 
      median=mqi[2], 
      ub=mqi[3]  
    )  # on data /user scale not internal link
    
    mqi = NULL

    O[["predictions_range"]] = range( which( inla_args[["data"]][["tag"]] == "predictions" ) )
    O[["priors"]] = H = inla_hyperparameters(  reference_sd = O[["data_range_internal"]][["sd"]], alpha=0.5, median(yl[ll], na.rm=TRUE) )  # sd slightly biased due to 0's being dropped .. but use of pc.priors that shrink to 0
    
    m = yl = ii = ll = fy = ol = NULL
    gc()
  
    if (exists("debug")) if (is.character(debug)) if ( debug =="fit") browser()
   
    # access H, cyclic_values, etc
    inla_args[[".parent.frame"]]=environment()
     

    # check INLA options
    if (!exists("inla.mode", inla_args)) inla_args[["inla.mode"]] = "compact" # changed inla options in 2023
    if (!exists("control.inla", inla_args)) inla_args[["control.inla"]] = list( strategy='adaptive', cmin=0 ) #int.strategy='eb'
    if (!exists("control.predictor", inla_args)) inla_args[["control.predictor"]] = list( compute=TRUE, link=1  ) #everything on link scale
    if (!exists("control.mode", inla_args ) ) inla_args[["control.mode"]] = list( restart=FALSE ) 
    if (!is.null(theta) ) inla_args[["control.mode"]]$theta= theta
    
    if (!exists("control.compute", inla_args)) inla_args[["control.compute"]] = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE )
      if ( inla_args[["inla.mode"]] == "classic") {
        # if ( !exists("control.results", inla_args ) ) inla_args[["control.results"]] = list(return.marginals.random=TRUE, return.marginals.predictor=TRUE )
        inla_args[["control.compute"]]$return.marginals.predictor  = TRUE  # location of this option has moved ... might move again
      }

    # if (!exists("control.fixed", inla_args)) inla_args[["control.fixed"]] = list(mean.intercept=0.0, prec.intercept=0.001, mean=0, prec=0.001)
    if (!exists("control.fixed", inla_args)) inla_args[["control.fixed"]] = H$fixed

    O[["inla.mode"]] = inla_args[["inla.mode"]]  # copy for later
    
    setDF(inla_args[["data"]]) # in case .. INLA requires this ?
 
 
    fit = try( do.call( inla, inla_args ) )      

    if (inherits(fit, "try-error" )) {
      inla_args[["control.inla"]] = list( int.strategy='eb', cmin=0 )
      fit = try( do.call( inla, inla_args ) )      
    }

    if (inherits(fit, "try-error" )) {
      inla_args[["safe"]] = TRUE
      inla_args[["inla.mode"]] = "classic"
      inla_args[["control.inla"]] = list( int.strategy='eb', cmin=0 )
      fit = try( do.call( inla, inla_args ) )      
    }

    if (inherits(fit, "try-error" )) stop( "Model fit error" )
 
    # to improve hyper param estimates..
    if (exists("improve.hyperparam.estimates", O)) {
      if (O[["improve.hyperparam.estimates"]]) {
        fit = inla.hyperpar(fit, dz=0.6, diff.logdens=9  )  # get improved estimates for the hyperparameters
      }
    }

    if (is.null(fit)) warning("model fit error")
    if ("try-error" %in% class(fit) ) warning("model fit error")


    run_fit_end  = Sys.time()

    message( "Fitted model saved as: ", fn_fit )

    outputdir = dirname(fn_fit)
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  
    # collect a few objects for ease of extraction

    setDT(inla_args[["data"]]) # revert to DT for faster / efficient operations

    # make these temporary indices here to drop inla_args and reduce RAM usage and make things easier later
 
    O[["ipred"]] = which( inla_args[["data"]][["tag"]]=="predictions" )  

    if (  O[["dimensionality"]] == "space" ) {
        # filter by S and T in case additional data in other areas and times are used in the input data
      O[["matchfrom"]] = list( 
        space=O[["space_id"]][inla_args[["data"]][[vS]][O[["ipred"]]]] 
      ) 
    }

    if (O[["dimensionality"]] == "space-time"  ) {
      O[["matchfrom"]] = list( 
        space=O[["space_id"]][inla_args[["data"]][[vS]][O[["ipred"]]]] , 
        time=O[["time_id"]][inla_args[["data"]][[vT]][O[["ipred"]]] ]
      )
    }

    if ( O[["dimensionality"]] == "space-time-cyclic" ) {
      O[["matchfrom"]] = list( 
        space=O[["space_id"]][inla_args[["data"]][[vS]][O[["ipred"]]]] , 
        time=O[["time_id"]][inla_args[["data"]][[vT]][O[["ipred"]]]],
        cyclic=O[["cyclic_id"]][inla_args[["data"]][[vU]][O[["ipred"]]]]
      )
    }

    # 2023 change in behaviour .. predictions do not require offsets in experimental and classical modes .. already incorporated
    # if (!is.null(vO)) {
    #   if ( O[["inla.mode"]] == "experimental" ) {
    #     O[["Offset"]] = inla_args[["data"]][[vO]][ O[["ipred"]] ] 
    #   }
    # }

    O[["carstm_prediction_surface_parameters"]] = NULL
    
    fit$modelinfo = O  # store in case a restart is needed

    fit$.args = NULL

    inla_args= NULL; gc()

    carstm_saveRDS( fit, file=fn_fit, compress=compress, compression_level=compression_level )

  }


  
  ### END Modelling

  if (!exists("run_fit_end")) run_fit_end  = run_start  # in case modelling is skipped

  run_extract_start  = Sys.time()

  if (is.null(fit)) fit =readRDS( fn_fit )
  
  if (is.null(fit)) {
    message( "fit file not found: ", fn_fit )
    stop()
  }
 
  # do the computations here as fit can be massive ... best not to copy, etc ..
  if (be_verbose)  message( "\nComputing summaries and extracting posterior simulations ..." )

  if (!exists("modelinfo", fit)) stop("modelinfo not in fit, fit object needs to be re-run")
  O = fit$modelinfo
   
  fit$modelinfo = NULL
  gc()


  if (exists("debug")) if (is.character(debug)) if (debug=="extract") browser()


  if (exists("data_transformation", O))  {
    backtransform = function( b ) {
      b[,1] =  O$data_transformation$backward( b[,1]   )
      return( b )
    }
  } 

  list_simplify = function(x) as.data.frame( t( as.data.frame( x )))
  exceedance_prob = function(x, threshold)  {1 - inla.pmarginal(q = threshold, x)}
  deceedance_prob = function(x, threshold)  { inla.pmarginal(q = threshold, x)}

  
  # local functions
  list_to_dataframe = function(Y ){
    if (!is.vector(Y)) {
      Y = as.data.frame(Y)
    }
    Z = data.frame(lapply(Y, function(x) Reduce(c, x)))
    row.names(Z) = row.names(Y)
    names(Z) = names(Y)
    return(Z)
  } 


  apply_generic = function(...)  mclapply(...,   mc.cores=mc.cores ) # drop-in for lapply
  apply_simplify = function(...) simplify2array(mclapply(...,  mc.cores=mc.cores ), higher = FALSE )  # drop in for sapply
  
  # not used, for reference
  apply_generic_serial = function(...)  lapply(...  ) # drop-in for lapply .. serial
  apply_simplify_serial = function(...) simplify2array(lapply(...  ))  # drop in for sapply .. serial

  sqrt_safe = function( a, eps=eps )  sqrt( pmin( pmax( a, eps ), 1/eps ) )


  marginal_clean = function( w) {
    i <- which(!is.finite(rowSums(w)) )
    if ( length(i) > 0) w = w[-i,] 
    w = w[order(w[,1]),]
    return(w)
  }

  test_for_error = function( Z ) {
    if ( "try-error" %in%  class(Z) ) return("error")
    if (is.list(Z)) {
      if (any( unlist(lapply(Z, function(o) inherits(o, "try-error"))) )) return("error")
    } else if (is.vector(Z) ){
      if (any( inherits(Z, "try-error")))  {
        return("error")
      } else if (any(grepl("Error", m))) {
        return("error")
      }  
    }
    return( "good" )    
  }

  marginal_summary = function(Z, invlink=NULL ) {
    
    if (!is.null(invlink)) {
      Z = try( apply_generic( Z, inla.tmarginal, fun=invlink, n=4096L) )
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      Z = try( apply_generic( Z, marginal_clean ) )
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      Z = try( apply_generic( Z, inla.zmarginal, silent=TRUE  ) )
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      Z = try( simplify2array( Z ), silent=TRUE) 
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      Z = try( list_simplify( Z), silent=TRUE )
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      return(Z)
    
    } else {

      Z = try( apply_generic( Z, inla.zmarginal, silent=TRUE  ), silent=TRUE)
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }
    
      Z = try( simplify2array( Z ), silent=TRUE) 
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      Z = try( list_simplify( Z ), silent=TRUE)
      if (test_for_error(Z) =="error") {
        class(m) = "try-error"
        return(m) 
      }

      return(Z)
    }
  }


  inla_tokeep = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  tokeep =      c("mean", "sd", "quant0.025", "quant0.5", "quant0.975")
 
  if (is.null(deceedance_threshold)) if (exists("deceedance_threshold", O)) deceedance_threshold = O[["deceedance_threshold"]]
  if (is.null(exceedance_threshold)) if (exists("exceedance_threshold", O)) exceedance_threshold = O[["exceedance_threshold"]]

  if (is.null(deceedance_threshold_predictions)) if (exists("deceedance_threshold_predictions", O)) deceedance_threshold_predictions = O[["deceedance_threshold_predictions"]]
  if (is.null(exceedance_threshold_predictions)) if (exists("exceedance_threshold_predictions", O)) exceedance_threshold_predictions = O[["exceedance_threshold_predictions"]]

  if (exists("debug")) if (is.character(debug)) if (debug=="summary") browser()

  if (!is.null(posterior_simulations_to_retain)) {
 
    posterior_simulations_to_retain = unique( c( "summary", posterior_simulations_to_retain))
 
    if (is.null(nposteriors)) {
      if (exists("nposteriors", O)) {
        nposteriors = O$nposteriors
      } else {
        message("nposteriors not found, defaulting to 1000")
        nposteriors = 1000 
      }
    }
    message( "Sampling from joint posteriors: n = ", nposteriors )
 
    S = try( inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores ), silent=TRUE)
    if ( "try-error" %in%  class(S) ) {
      S = try( inla.posterior.sample( nposteriors, fit, add.names=FALSE, use.improved.mean=FALSE, skew.corr =FALSE, seed=123L ), silent=TRUE)
    }
    if ( "try-error" %in%  class(S) ) stop("posterior sampling error")

    message( "Sampling complete ... now reformatting parameters and extracting required components" )

    for (z in c("tag", "start", "length") ) assign(z, attributes(S)[[".contents"]][[z]] )  # index info

      if (0) {
        # not used at the moment .. 
        nlogdens = length(S[[1]]$logdens)
        logdens = array(NA, dim=c( nlogdens, nposteriors  ) )
        for (i in 1:nposteriors) {
          logdens[,i] = unlist(S[[i]]$logdens)
        }
        logdens_names =  names(S[[1]]$logdens)
        logdens = format_results( logdens, labels=logdens_names  ) # same seq as space_id ( == attributes(space)$row.names )
      }
  }

  names.fixed = fit$names.fixed
  summary.fixed  = fit$summary.fixed
  summary.random = fit$summary.random
  summary.hyperpar = fit$summary.hyperpar 
  
  dic = fit$dic[c("dic", "p.eff", "dic.sat", "mean.deviance")]
  waic = fit$waic[c("waic", "p.eff")]
  mlik = fit$mlik[2]


  if (!exists("summary", O)) O[["summary"]] = list()
  if (!exists("random", O)) O[["random"]] = list()


  if ( "summary" %in% toget) {

    O[["summary"]][["direct"]] = summary(fit)
    # print(O[["summary"]][["direct"]])

    # remove a few unused but dense data objects
    O[["summary"]][["direct"]]$linear.predictor = NULL
    O[["summary"]][["direct"]]$waic$local.waic = NULL
    O[["summary"]][["direct"]]$waic$local.p.eff = NULL
    O[["summary"]][["direct"]]$dic$local.dic = NULL
    O[["summary"]][["direct"]]$dic$local.p.eff = NULL
    O[["summary"]][["direct"]]$dic$local.dic.sat = NULL
    O[["summary"]][["direct"]]$dic$family = NULL
    
    # extract and back transform where possible
    if (exists( "marginals.fixed", fit)) {
      m = fit$marginals.fixed  # make a copy to do transformations upon
      fi = grep("Intercept", names(m) )
 
      m = try( apply_generic( m, function(x) marginal_clean( inla.tmarginal( invlink, x, n=4096L)) ), silent=TRUE  )

      if ( exists("data_transformation", O))  {
        m[[fi]] = inla.tmarginal( O$data_transformation$backward, m[[fi]] , n=4096L ) # on user scale
      }

      m = try(apply_simplify( m, FUN=inla.zmarginal, silent=TRUE ), silent=TRUE)
      if (test_for_error(m) =="error") { 
          message("Problem with marginals for intercept (probably too variable), doing a simple inverse transform, SD is blanked out")
          W = fit$summary.fixed[ "(Intercept)", 1:5]
          colnames(W) = tokeep
          W[,c(1,3:5)] = exp( W[,c(1,3:5)] )
          W[,2] = NA

      } else {
        W = cbind ( t (m) )  # 
        W = list_to_dataframe( W [, tokeep, drop =FALSE] )
      }

      W$ID = row.names(W)
      O[["summary"]][["fixed_effects"]] = W
      m = NULL
      
      if ( "summary" %in% posterior_simulations_to_retain ) {
       # posteriors
        ftx = names.fixed 
        flabels= ftx
        nfixed = length(names.fixed)
        fixed = array(NA, dim=c( nfixed, nposteriors  ) )
        fkk = inla_get_indices(ftx, tag=tag, start=start, len=length, model="direct_match")
        fkk = unlist(fkk)
        for (i in 1:nposteriors) {
          fixed[,i] = S[[i]]$latent[fkk,]  
        }
        fixed = invlink(fixed)
        row.names(fixed) = flabels
        O[["sims"]][["fixed_effects"]] = fixed
        O[["posterior_summary"]][["fixed_effects"]] = posterior_summary(format_results( fixed, labels=flabels))
        fkk = fixed = NULL
      }      
    }
    W = NULL
    gc()
   
 
    if (exists( "marginals.hyperpar", fit)) {
      
      # hyperpar (variance components)
      hyps =  row.names(fit$summary.hyperpar)

      prcs = grep( "^Precision.*", hyps, value=TRUE )
      if (length(prcs) > 0) {
  
        m = fit$marginals.hyperpar[prcs]
        m = try( apply_generic( m, inla.tmarginal, fun=function(y) 1/sqrt_safe( y, eps ), n = 4096L ), silent=TRUE)
        m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
        m = try( simplify2array( m ), silent=TRUE)
        if (test_for_error(m) =="error") {  
          if (be_verbose)  {
            message( "NAN or Inf values encountered in marginals.") 
            message( "Try an alternate parameterization as model may be over parameterized or degenerate. ")
            message( "Copying fit summaries directly rather than from marginals ... ")
          }
          m = fit$summary.hyperpar[prcs,1:5]
          m[,c(1,3:5)] = 1/sqrt_safe( m[,c(1,3:5)], eps )
          colnames(m) = tokeep
        }
        m = try( list_simplify( m), silent=TRUE)
        rownames(m) = gsub("Precision for", "SD", rownames(m) )
        rownames(m) = gsub(" for", "", rownames(m) )
        rownames(m) = gsub(" the", "", rownames(m) )
        
        O[["summary"]][["random_effects"]] = m[, tokeep, drop =FALSE] 
      }
 
      m = NULL

      # update phi's, lambda's (used in besagproper2 -- Leroux model) .. etc
      rhos = grep( "^Rho.*|^GroupRho.*", hyps, value=TRUE )
      phis = grep( "^Phi.*", hyps, value=TRUE )
      other = grep( "^Lambda.*|^Diagonal.*|zero-probability.*", hyps, value=TRUE )

      known = c( rhos, phis, other )
      unknown = setdiff( hyps, c(prcs, known) )

      if (length(rhos) > 0) {
        m = marginal_summary( fit$marginals.hyperpar[ rhos ] )
        
        if (test_for_error(m) =="error") {  
          m = fit$summary.hyperpar[rhos, 1:5]
          colnames(m) = tokeep
        }
        
        O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]],  m[, tokeep, drop =FALSE] )

      }

      if (length(phis) > 0) {

        m = marginal_summary( fit$marginals.hyperpar[ phis ] )
        
        if (test_for_error(m) =="error") {  
          m = fit$summary.hyperpar[phis, 1:5]
          colnames(m) = tokeep
        }

        O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], m[, tokeep, drop =FALSE] )

      }

      if (length(unknown) > 0) {
        m = marginal_summary( fit$marginals.hyperpar[ unknown ] )
        
        if (is.character(m) && m=="error")  {
          #  alternatively: m[,"mode"] = apply_simplify( fit$marginals.hyperpar[ unknown ], FUN=function(x) inla.mmarginal( x ))
          m = fit$summary.hyperpar[unknown, 1:5]
          colnames(m) = tokeep
        }

        O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], m[, tokeep, drop =FALSE])

      }

      O[["summary"]][["random_effects"]] = list_to_dataframe(  O[["summary"]][["random_effects"]] )
      O[["summary"]][["random_effects"]]$ID = row.names( O[["summary"]][["random_effects"]] )
   

      if (length(fit[["marginals.random"]]) > 0) { 

        if (be_verbose)  message("Extracting random effects of covariates, if any" )
        if (exists("debug")) if (is.character(debug)) if ( debug =="random_covariates") browser()
 
        raneff = setdiff( names( fit$marginals.random ), c(vS, vS2, vS3  ) )
        for (rnef in raneff) {
          m = marginal_summary( fit$marginals.random[[rnef]], invlink=invlink )
          if (test_for_error(m) =="error") {  
             message( "failed to transform marginals .. copying directly from INLA summary instead: ", rnef)
            m = fit$summary.random[[rnef]][, inla_tokeep, drop =FALSE ]
            names(m) =  tokeep
            O[["random"]] [[rnef]] = m
            O[["random"]] [[rnef]]$ID = fit$summary.random[[rnef]]$ID
            O[["random"]] [[rnef]] = list_to_dataframe( O[["random"]] [[rnef]] )
          } else {
            O[["random"]] [[rnef]] = m[, tokeep, drop =FALSE]
            O[["random"]] [[rnef]]$ID = fit$summary.random[[rnef]]$ID
            O[["random"]] [[rnef]] = list_to_dataframe( O[["random"]] [[rnef]] )
          }
        }
        m = raneff = NULL
      }
    }


   if ( "summary" %in% posterior_simulations_to_retain ) {

      other_random = setdiff( names(summary.random), c(vS, vS2 ) )

      if (length(other_random) > 0 ) {
        nrandom = sum(sapply(summary.random[other_random], nrow))  
        random = array(NA, dim=c(nrandom, nposteriors  ) )
        rkk = inla_get_indices(other_random, tag=tag, start=start, len=length, model="direct_match"  )  # if bym2, must be decomposed  
        rkk = unlist( rkk )
        for (i in 1:nposteriors) {
          random[,i] = S[[i]]$latent[rkk,]  
        }
        rlabels = names(S[[1]]$latent[rkk,])

        random = invlink(random)
        row.names(random) = rlabels
        O[["posterior_summary"]][["random_effects"]] = posterior_summary( format_results( random, labels=rlabels ) )
      }

      nhyperpar = nrow(summary.hyperpar)
      hyperpar = array( NA,dim=c(nhyperpar, nposteriors  ) )
      for (i in 1:nposteriors) {
        hyperpar[,i] = S[[i]]$hyperpar
      }
      hyper_names = names(S[[1]]$hyperpar)
      k = grep("Precision", hyper_names)
      if (length(k) > 0) {
          hyperpar[k,] = 1 / sqrt(hyperpar[k,]) 
          hyper_names[k] = gsub("Precision for", "SD", hyper_names[k] )
      }
      hyper_names = gsub("for ", "", hyper_names )
      hyper_names = gsub("the ", "", hyper_names )
      
      row.names( hyperpar ) = hyper_names
      O[["sims"]][["hyperpars"]] =  hyperpar

      O[["posterior_summary"]][["hyperpars_summary"]] = posterior_summary( format_results( hyperpar, labels=hyper_names) )

    }


    if (be_verbose)  {
      # message( "")
      # message( "Random effects:")
      # print(  O[["summary"]][["random_effects"]]  )   
      # message( "\n--- NOTE --- 'SD *' from marginal summaries are on link scale")
      # message( "--- NOTE --- SD * from posteriors simulations are on user scale")
      # message( "")
    }

  }  # end random_effects
 
  random = hyperpar = rkk=  NULL
  k = phis = rhos = known = unknown = m  = NULL

  gc()

 
  # separate out random spatial and randomm spatiotemporal (as they can be large arrays)
  if ("random_spatial" %in% toget) {
    # space only
    Z = NULL
    iSP = which( re$dimensionality=="s" & re$level=="main")
    if (length(iSP) > 0 ) {

      if (be_verbose)  message("Extracting random spatial errors"  )
      if (exists("debug")) if (is.character(debug)) if ( debug =="random_spatial") browser()

      matchto = list( space=O[["space_id"]] )

      W = array( NA, dim=c( O[["space_n"]], length(tokeep) ), dimnames=list( space=O[["space_name"]], stat=tokeep ) )
      names(dimnames(W))[1] = vS  # need to do this in a separate step ..

      O[["random"]] [[vS]] = list()  # space as a main effect  vS==vS

      if (length(iSP) == 1) {

        # vS = re$vn[ iSP ]  # == vS as it is a single spatial effect
        
        model_name = re$model[ iSP ]  # should be iid

        m = fit$marginals.random[[vS]]
 
        m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=4096L ) , silent=TRUE )
        m = try( apply_generic( m, inla.zmarginal, silent=TRUE ), silent=TRUE )
        m = try( simplify2array( m ), silent=TRUE)
        m = try( list_simplify( m ), silent=TRUE )
        # single spatial effect (eg in conjuction with besag) .. indexing not needed but here in case more complex models ..
        if (test_for_error(m) =="error") {  
          message( "failed to transform random_spatial marginals .. copying directly from INLA summary instead")
          m = fit$summary.random[[vS]][, inla_tokeep ]
          names(m) =  tokeep
        } 
        
        if ( model_name == "bym2" ) {
          # bym2 effect is coded by INLA as a double length vector: bym and iid simultaneously
          # this first part captures the iid part while the part outside of the if * captures the bym
          Z = expand.grid( space=O[["space_id"]], type = c("iid", model_name), stringsAsFactors =FALSE )

          #  extract iid main effects
          iid = which(Z$type=="iid")
          matchfrom = list( space=Z[["space"]][iid] )

          for (k in 1:length(tokeep)) {
            W[,k] = reformat_to_array( input = unlist(m[iid, tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vS]] [["iid"]] = W [, tokeep, drop =FALSE]

        } else {
          # single spatial effect that is not bym2
          # this is redundant with iSP being a single factor, but approach is generalizable for higher dims 
          Z = expand.grid( space=O[["space_id"]], type=model_name, stringsAsFactors =FALSE )
        }

        bym2 =  which(Z$type==model_name)
        matchfrom = list( space=Z[["space"]][bym2] )
        W[] = NA
        for (k in 1:length(tokeep)) {
          W[,k] = reformat_to_array( input = unlist(m[bym2, tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
        }
        O[["random"]] [[vS]] [[model_name]] = W [, tokeep, drop =FALSE]

      }

      if (length(iSP) == 2) {
        
        for (j in 1:length(iSP)) {
          
          vnS = re$vn[ iSP[j] ]
          model_name = re$model[ iSP[j] ]  
          
          m = fit$marginals.random[[vnS]]
         
          m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=4096L), silent=TRUE  )
          m = try( apply_generic( m, inla.zmarginal, silent=TRUE) , silent=TRUE )
          m = try( simplify2array( m ), silent=TRUE) 
          m = try( list_simplify(m), silent=TRUE  )
          if (test_for_error(m) =="error") {  
            message( "failed to transform random_spatial marginals .. copying directly from INLA summary instead")
            m = fit$summary.random[[vnS]][, inla_tokeep ]
            names(m) = c("ID", tokeep)
          } 

          # single spatial effect (eg besag, etc)
          matchfrom = list( space=O[["space_id"]] )
          W[] = NA
          for (k in 1:length(tokeep)) {
            W[,k] = reformat_to_array( input = unlist(m[, tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vnS]] [[model_name]] = data.frame( W [, tokeep, drop =FALSE], ID=row.names(W) )

        }
      }

      Z = m = matchfrom = NULL
      gc()


      if ( "random_spatial" %in% posterior_simulations_to_retain ) {

        # POSTERIOR SIMS
        slabels = O[["space_name"]]

        space = array(NA, dim=c( O$space_n, nposteriors  ) )
        row.names(space) = slabels

        space1 = space2 = NULL

        if (length(iSP) == 1) {

          stx1 = paste("^", re$vn[iSP], "$", sep="")
          if (re$model[iSP] == "bym2") {
            # special case bym2 has two factors rolled together
            space1 = array(NA, dim=c( O$space_n, nposteriors  ) )
            space2 = array(NA, dim=c( O$space_n, nposteriors  ) )
            row.names(space1) = slabels
            row.names(space2) = slabels

            skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length, model="bym2" )  # if bym2, must be decomposed  
            for (i in 1:nposteriors) {
              space1[,i] = S[[i]]$latent[skk1[["iid"]],] 
              space2[,i] = S[[i]]$latent[skk1[["bym"]],]
            }      
            space = space1 + space2

          } else {
            # single spatial effect of some kind
            skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
            skk1 = unlist(skk1)
            for (i in 1:nposteriors) {
              space[,i] = S[[i]]$latent[skk1,] 
            }      
          }
        }
 
        if (length(iSP) == 2) {
          space1 = array(NA, dim=c( O$space_n, nposteriors  ) )
          space2 = array(NA, dim=c( O$space_n, nposteriors  ) )
          row.names(space1) = slabels
          row.names(space2) = slabels

          i1 = which(re$model[iSP[1]] == "iid")
          i2 = which(re$model[iSP[2]] %in% c("besag", "bym") )    # add others as required
          if (length(i1)==0 | length(i2)==0) stop( "Unexpected situation: two spatial effects found, expecting one to be iid, and a second to be besag or bym, but it was not." )
          stx1 = paste("^", re$vn[iSP[i1]], "$", sep="")
          stx2 = paste("^", re$vn[iSP[i2]], "$", sep="")
          skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
          skk2 = inla_get_indices(stx2, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
          for (i in 1:nposteriors) {
            space1[,i] = S[[i]]$latent[skk1[[1]],]  
            space2[,i] = S[[i]]$latent[skk2[[1]],]
          }      
          space = space1 + space2
        }

        space = invlink(space) 
        row.names(space) = slabels
        matchfrom = list( space=O[["space_id"]] )

        if (!is.null(exceedance_threshold)) {
          if (be_verbose)  message("Extracting random spatial errors exceedence"  )

          for ( b in 1:length(exceedance_threshold)) {
            m = apply ( space, 1, FUN=function(x) length( which(x > exceedance_threshold[b]) ) ) / nposteriors
            m = reformat_to_array( input = m, matchfrom=matchfrom, matchto = matchto )
            names(dimnames(m))[1] = vS
            dimnames( m )[[vS]] = O[["space_name"]]
            O[["random"]] [[vS]] [["exceedance"]] [[as.character(exceedance_threshold[b])]] = data.frame( m [, tokeep, drop =FALSE], ID=row.names(m) )
            m = NULL
          }
        }

        if (!is.null(deceedance_threshold)) {
          if (be_verbose)  message("Extracting random spatial errors deceedance"  )
          # redundant but generalizable to higher dims
          for ( b in 1:length(deceedance_threshold)) {
            m = apply ( space, 1, FUN=function(x) length( which(x < deceedance_threshold[b]) ) ) / nposteriors
            m = reformat_to_array( input = m, matchfrom=matchfrom, matchto=matchto  )
            names(dimnames(m))[1] = vS
            dimnames( m )[[vS]] = O[["space_name"]]
            O[["random"]] [[vS]] [["deceedance"]] [[as.character(deceedance_threshold[b])]] = data.frame( m [, tokeep, drop =FALSE], ID=row.names(m) ) 
            m = NULL
          }
        }

        m = posterior_summary( format_results( space, labels=slabels  ) )
        W[] = NA
        for (k in 1:length(tokeep)) {
          W[,k] = reformat_to_array(  input = m[, tokeep[k]], matchfrom=matchfrom, matchto=matchto )
        }
        O[["random"]] [[vS]] [["combined"]] = W[, tokeep, drop =FALSE] 

        if ( "random_spatial" %in% posterior_simulations_to_retain ) {
          O[["sims"]] [[vS]] [["combined"]] = space  # already inverse link scale
        }

        if ( "random_spatial12" %in% posterior_simulations_to_retain ) {
          if (!is.null(space1)) O[["sims"]] [[vS]] [["iid"]]  =  invlink(space1) 
          if (!is.null(space2)) O[["sims"]] [[vS]] [["bym2"]] =  invlink(space2) 
        }
      }
    }
  }  # end random spatial effects

  matchfrom = i1 = i2= NULL
  Z = W = m = space = space1 = space2 = skk1 = skk2 = iSP = NULL
  gc()
 

  if ("random_spatiotemporal" %in% toget ) {
    # space-time

    iST = which( re$dimensionality=="st" & re$level=="main")
    
    if (length(iST) > 0 ) {

      if (be_verbose)  message("Extracting random spatiotemporal errors"  )

      if (exists("debug")) if (is.character(debug)) if ( debug =="random_spatiotemporal") browser()

      matchto = list( space=O[["space_id"]], time=O[["time_id"]]  )
      
      W = array( NA, 
        dim=c( O[["space_n"]], O[["time_n"]], length(tokeep) ), 
        dimnames=list( space=O[["space_name"]], time=O[["time_name"]], stat=tokeep ) )
      names(dimnames(W))[1] = vS  # need to do this in a separate step ..
      names(dimnames(W))[2] = vT  # need to do this in a separate step ..

      if (length(iST) == 1) {

        vnST = re$vn[ iST ]
        model_name = re$model[ iST ]   

        if (exists(vnST, fit$marginals.random )) {
 
          m = fit$marginals.random[[vnST]]
          m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=4096L) , silent=TRUE )

          m = try( apply_generic( m, inla.zmarginal, silent=TRUE), silent=TRUE )
          m = try( simplify2array( m ), silent=TRUE)
          m = try( list_simplify( m ), silent=TRUE )
          if (test_for_error(m) =="error") {  
            message( "failed to transform random_spatial marginals .. copying directly from INLA summary instead")
            m = fit$summary.random[[vnST]][, inla_tokeep ]
            names(m) = tokeep
          } 
        
          if ( model_name == "bym2" ) {
            # bym2 effect: bym and iid with annual results
            Z = expand.grid( space=O[["space_id"]], type = c("iid", model_name), time=O[["time_id"]], stringsAsFactors =FALSE )

            #  spatiotemporal interaction effects  iid
            iid = which(Z$type=="iid")
            matchfrom = list( space=Z[["space"]][iid], time=Z[["time"]][iid]  )

            for (k in 1:length(tokeep)) {
              W[,,k] = reformat_to_array(  input = unlist(m[iid, tokeep[k]]), matchfrom = matchfrom, matchto = matchto )
            }
            O[["random"]] [[vnST]] [["iid"]] =  W [,, tokeep, drop =FALSE] 

          } else {
            # besag effect: with annual results
            Z = expand.grid( space=O[["space_id"]], type =model_name, time=O[["time_id"]], stringsAsFactors =FALSE )
          }

          #  spatiotemporal interaction effects  bym
          bym2 =  which(Z$type==model_name)
          matchfrom = list( space=Z[["space"]][bym2], time=Z[["time"]][bym2]  )
            
          for (k in 1:length(tokeep)) {
            W[,,k] = reformat_to_array( input = unlist(m[bym2,tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vnST]] [[model_name]] =   W [,, tokeep, drop =FALSE] 

        }
      }
      
      if (length(iST) == 2) {

        for (j in 1:length(iST)) {
          vnST = re$vn[ iST[j] ]
          model_name = re$model[ iST[j] ]  
 
          m = fit$marginals.random[[vnST]]
          m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=4096L), silent=TRUE )
          m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE )
          m = try( simplify2array( m ), silent=TRUE)
          m = try( list_simplify( m ), silent=TRUE )
          if (test_for_error(m) =="error") { 
            message( "failed to transform random_spatiotemporal marginals .. copying directly from INLA summary instead")
            m = fit$summary.random[[vnST]][, inla_tokeep ]
            names(m) = tokeep
          } 
          
          Z = expand.grid( space=O[["space_id"]], type =model_name, time=O[["time_id"]], stringsAsFactors =FALSE )
          jst =  which(Z$type==model_name)
          matchfrom = list( space=Z[["space"]][jst], time=Z[["time"]][jst]  )
          for (k in 1:length(tokeep)) {
            W[,,k] = reformat_to_array( input = unlist(m[jst, tokeep[k] ]), matchfrom = matchfrom, matchto = matchto  )
          }
          O[["random"]] [[vnST]] [[model_name]] = W [,, tokeep, drop =FALSE]
          m = NULL
        }
      }

      if ( "random_spatiotemporal" %in% posterior_simulations_to_retain ) {
        # posterior simulations

        space_time1 = space_time2  = space_time = array( NA, 
          dim=c( O[["space_n"]] * O[["time_n"]] , nposteriors  ) )
        L = CJ( time=O[["time_id"]], space=O[["space_id"]] )  # note:: CJ has reverse order vs expand.grid
        stlabels = paste(L[["space"]], L[["time"]], sep="_")

        if (length(iST) == 1) {
          stx1 = paste("^", re$vn[iST], "$", sep="")
          if (re$model[iST] == "bym2") {
            # special case bym2 has two factors rolled together
            skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length, model="bym2" )  # if bym2, must be decomposed  
            for (i in 1:nposteriors) {
              space_time1[,i] = S[[i]]$latent[skk1[["iid"]],]  
              space_time2[,i] = S[[i]]$latent[skk1[["bym"]],]
            }    
            space_time = space_time1 + space_time2  
            row.names(space_time1) = stlabels
            row.names(space_time2) = stlabels

          } else {
            # single spatial effect of some kind
            skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
            skk1 = unlist(skk1)
            for (i in 1:nposteriors) {
              space_time[,i] = S[[i]]$latent[skk1,] 
            }      
          }
        }
        
        if (length(iST) == 2) {
          i1 = which(re$model[iST] == "iid")
          i2 = which(re$model[iST] %in% c("besag", "bym") )   # add others as required
          if (length(i1)==0 | length(i2)==0) stop( "Unexpected situation: two spatial-time effects found, expecting one to be iid, and a second to be besag or bym, but it was not." )
          stx1 = paste("^", re$vn[iST[i1]], "$", sep="")
          stx2 = paste("^", re$vn[iST[i2]], "$", sep="")
          skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
          skk2 = inla_get_indices(stx2, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
          for (i in 1:nposteriors) {
            space_time1[,i] = S[[i]]$latent[skk1[[1]],] 
            space_time2[,i] = S[[i]]$latent[skk2[[1]],]
          }      
          space_time  = space_time1 + space_time2
          row.names(space_time1) = stlabels
          row.names(space_time2) = stlabels
        }
      
        space_time = invlink(space_time)
        row.names(space_time) = stlabels

        Z = expand.grid( space=O[["space_id"]], type =model_name, time=O[["time_id"]], stringsAsFactors =FALSE )
        jst =  which(Z$type==model_name)
        matchfrom = list( space=Z[["space"]][jst], time=Z[["time"]][jst]  )

        if (!is.null(exceedance_threshold)) {
          if (be_verbose)  message("Extracting random spatiotemporal errors exceedence"  )

          for ( b in 1:length(exceedance_threshold)) {
            m = apply ( space_time, 1, FUN=function(x) length( which(x > exceedance_threshold[b] ) ) ) / nposteriors
            m = reformat_to_array( input=m, matchfrom=matchfrom,  matchto=matchto )
            names(dimnames(m))[1] = vS
            names(dimnames(m))[2] = vT
            dimnames( m )[[vS]] = O[["space_name"]]
            dimnames( m )[[vT]] = O[["time_name"]]
            O[["random"]] [[vnST]] [["exceedance"]] [[as.character(exceedance_threshold[b])]] = m
          }
          m = NULL
        }

        if (!is.null(deceedance_threshold)) {
          if (be_verbose)  message("Extracting random spatiotemporal errors deceedance"  )

          for ( b in 1:length(deceedance_threshold)) {
            m = apply ( space_time, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
            m = reformat_to_array( input = m, matchfrom=matchfrom,  matchto=matchto )
            names(dimnames(m))[1] = vS
            names(dimnames(m))[2] = O[["time_name"]]
            dimnames( m )[[vS]] = O[["space_name"]]
            dimnames( m )[[vT]] = vT
            O[["random"]] [[vnST]] [["deceedance"]] [[as.character(deceedance_threshold[b])]] = m
          }
          m = NULL
        }
  
        W[] = NA
        m = posterior_summary(format_results( space_time, labels=stlabels ))
        for (k in 1:length(tokeep)) {
          W[,,k] = reformat_to_array(  input = m[, tokeep[k]], matchfrom=matchfrom, matchto=matchto )
        }
        O[["random"]] [[vnST]] [["combined"]] = W[,, tokeep, drop =FALSE] 

        if ( "random_spatiotemporal" %in% posterior_simulations_to_retain ) {
          O[["sims"]] [[vnST]] [["combined"]] = space_time
        }
        if ( "random_spatiotemporal12" %in% posterior_simulations_to_retain ) {
          if (!is.null(space_time1)) O[["sims"]] [[vnST]] [["iid"]] =  invlink(space_time1)
          if (!is.null(space_time2)) O[["sims"]] [[vnST]] [["bym2"]] = invlink(space_time2)
        }
      } 
      Z = W = m = space_time = space_time1 = space_time2 = skk1 = skk2 = NULL
      gc()

    }
  }  # end random spatio-temporal effects

  run_predict_start  = Sys.time()


  if ("predictions" %in% toget ) {

    # marginals.fitted.values are on response scale (for both experimental and classical) and  already incorporates offsets
    # summary.fitted.values on response scale ((for both experimental and classical)) 
    # posteriors are on link scale and already incorporates offsets

    if (be_verbose) message("Extracting posterior predictions" )

    if (exists("debug")) if (is.character(debug)) if ( debug =="predictions") browser()

    # prepapre prediction simulations (from joint posteriors)
    ptx = "^Predictor"
 
    preds = O[["predictions_range"]][1]:O[["predictions_range"]][2]  
    npredictions = length(preds)
     
    pkk = inla_get_indices(ptx, tag=tag, start=start, len=length)
    pkk = unlist(pkk)[preds]
    predictions = array(NA, dim=c( npredictions, nposteriors  ) )
    for (i in 1:nposteriors)  predictions[,i] = S[[i]]$latent[pkk, ]
    predictions = invlink(predictions)  #  required for both classical and experimental

    pkk = preds = S = NULL

    gc()

    if ( exists("data_transformation", O))  predictions = O$data_transformation$backward( predictions  )
 
    if (!exists("predictions", O)) O[["predictions"]] = list()
  
    if (exists("marginals.fitted.values", fit)) {

    if (length(fit[["marginals.fitted.values"]]) > 0 ) {

      if ( O[["dimensionality"]] == "space" ) {

        m = fit$marginals.fitted.values[O[["ipred"]]]  # already incorporates offsets

        # 2023 change in behaviour .. predictions do not require offsets in experimental and classical modes .. already incorporated
        # if (!is.null(vO)) {
        #   if ( O[["inla.mode"]] == "experimental" ) {
        #     for ( i in 1:length(O[["ipred"]]) ) m[[i]][,1] = m[[i]][,1] + O[["Offset"]][i]
        #   } 
        # }
        
        # 2023 change in behaviour .. no longer needs invlink as marginals.fitted.values is on response scale
        # m = apply_generic( m, function(u) {inla.tmarginal( invlink, u) } )

        if ( exists("data_transformation", O))  m = apply_generic( m, backtransform )

        m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
        m = try( list_simplify( simplify2array( m ) ), silent=TRUE)
        if (test_for_error(m) =="error") {  
            message( "failed to summarize marginals.fitted.values .. copying directly from INLA summary instead")
            m = fit$summary.fitted.values[, inla_tokeep ]
            names(m) = tokeep
        } 

        W = array( NA, 
          dim=c( O[["space_n"]],  length(names(m)) ),  
          dimnames=list( space=O[["space_name"]], stat=names(m) ) )
        names(dimnames(W))[1] = vS  # need to do this in a separate step ..
        
        # matchfrom already created higher up
        matchto = list( space=O[["space_id"]] )

        for (k in 1:length(names(m))) {
          W[,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
        }
        O[["predictions"]] = W[, tokeep, drop =FALSE]
        m = W = NULL

        if ( "predictions" %in% posterior_simulations_to_retain ) {
          W = array( NA, 
            dim=c( O[["space_n"]], nposteriors ),  
            dimnames=list( space=O[["space_name"]], sim=1:nposteriors ) )
  
          names(dimnames(W))[1] = vS  # need to do this in a separate step ..
          for (k in 1:nposteriors ) {
            W[,k] = reformat_to_array( input=unlist(predictions[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
          }
          O[["sims"]][["predictions"]] = W[, drop =FALSE]
          W = NULL 
        }
      
      }


      if (O[["dimensionality"]] == "space-time"  ) {

        m = fit$marginals.fitted.values[O[["ipred"]]]   

        # 2023 change in behaviour .. predictions do not require offsets in experimental and classical modes .. already incorporated
        # if (!is.null(vO)) {
        #   if ( O[["inla.mode"]] == "experimental" ) {
        #     for ( i in 1:length(O[["ipred"]]) ) m[[i]][,1] = m[[i]][,1] + O[["Offset"]][i]
        #   } 
        # }

        # 2023 change in behaviour .. no longer needs invlink as marginals.fitted.values is on response scale
        # m = apply_generic( m, function(u) {inla.tmarginal( invlink, u) } )    
 
        if (exists("data_transformation", O)) m = apply_generic( m, backtransform )
        m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
        m = try( list_simplify( simplify2array( m ) ), silent=TRUE)
        if (test_for_error(m) =="error") {  
            message( "failed to summarize marginals.fitted.values .. copying directly from INLA summary instead")
            m = fit$summary.fitted.values[, inla_tokeep ]
            names(m) = tokeep
          } 

        W = array( NA, 
          dim=c( O[["space_n"]], O[["time_n"]], length(names(m)) ),  
          dimnames=list( space=O[["space_name"]], time=O[["time_name"]], stat=names(m) ) 
        )
 
        names(dimnames(W))[1] = vS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vT  # need to do this in a separate step ..

        # matchfrom already created higher up
        matchto = list( space=O[["space_id"]], time=O[["time_id"]] )
        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=O[["matchfrom"]], matchto=matchto)
        }
        O[["predictions"]] = W[,, tokeep, drop =FALSE]
    
        m = W = NULL
 
        if ( "predictions" %in% posterior_simulations_to_retain ) {

          W = array( NA, 
            dim=c( O[["space_n"]], O[["time_n"]], nposteriors ),  
            dimnames=list( space=O[["space_name"]], time=O[["time_name"]], sim=1:nposteriors ) )

          names(dimnames(W))[1] = vS  # need to do this in a separate step ..
          names(dimnames(W))[2] = vT  # need to do this in a separate step ..

          for (k in 1:nposteriors ) {
            W[,,k] = reformat_to_array( input=unlist(predictions[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
          }
          
          O[["sims"]][["predictions"]] = W[,,, drop =FALSE]
          W = NULL
        }

      }


      if ( O[["dimensionality"]] == "space-time-cyclic" ) {

        m = fit$marginals.fitted.values[O[["ipred"]]]   

        # 2023 change in behaviour .. predictions do not require offsets in experimental and classical modes .. already incorporated
        # if (!is.null(vO)) {
        # if ( O[["inla.mode"]] == "experimental" ) {
        #   # offsets required in experimental mode ... do not ask me why
        #   for ( i in 1:length(O[["ipred"]]) ) m[[i]][,1] = m[[i]][,1] + O[["Offset"]][i] 
        # } 
        #}
 
        # 2023 change in behaviour .. no longer needs invlink as marginals.fitted.values is on response scale
        # m = apply_generic( m, function(u) {inla.tmarginal( invlink, u) } )    

        if (exists("data_transformation", O)) m = apply_generic( m, backtransform )

        m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
        m = try( list_simplify( simplify2array( m ) ), silent=TRUE)
        if (test_for_error(m) =="error") {  
            message( "failed to summarize marginals.fitted.values .. copying directly from INLA summary instead")
            m = fit$summary.fitted.values[, inla_tokeep ]
            names(m) = tokeep
          } 

        W = array( NA, 
          dim=c( O[["space_n"]], O[["time_n"]], O[["cyclic_n"]], length(names(m)) ),  
          dimnames=list( space=O[["space_name"]], time=O[["time_name"]], cyclic=O[["cyclic_name"]], stat=names(m) ) )
        names(dimnames(W))[1] = vS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vT  # need to do this in a separate step ..
        names(dimnames(W))[3] = vU  # need to do this in a separate step ..

        matchto = list( space=O[["space_id"]], time=O[["time_id"]], cyclic=O[["cyclic_id"]] )

        for (k in 1:length(names(m))) {
          W[,,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
        }
        O[["predictions"]] = W[,,, tokeep, drop =FALSE]
        m = W = NULL
      
        if ( "predictions" %in% posterior_simulations_to_retain ) {
          
          W = array( NA, 
            dim=c( O[["space_n"]], O[["time_n"]], O[["cyclic_n"]], nposteriors ),  
            dimnames=list( space=O[["space_name"]], time=O[["time_name"]], cyclic=O[["cyclic_name"]], sim=1:nposteriors ) )

          names(dimnames(W))[1] = vS  # need to do this in a separate step ..
          names(dimnames(W))[2] = vT  # need to do this in a separate step ..
          names(dimnames(W))[3] = vU  # need to do this in a separate step ..
          for (k in 1:nposteriors ) {
            W[,,,k] = reformat_to_array( input=unlist(predictions[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
          }
          O[["sims"]][["predictions"]] = W[,,, ,drop =FALSE]
          W = NULL
        
        }

      }

 
      if ( "predictions" %in% posterior_simulations_to_retain ) {
    
        if (!is.null(exceedance_threshold_predictions)) {
          if (be_verbose)  message("Extracting de/exceedance of predictions"  )

          for ( b in 1:length(exceedance_threshold_predictions)) {
              m = apply ( predictions, 1, FUN=function(x) length( which(x > exceedance_threshold_predictions[b] ) ) ) / nposteriors
              W = reformat_to_array( input=m, matchfrom=O[["matchfrom"]],  matchto=matchto )
              names(dimnames(W))[1] = vS
              dimnames( W )[[vS]] = O[["space_name"]]
              m = NULL
              if ( length(O[["matchfrom"]]) > 1 ) {
                names(dimnames(W))[2] = vT
                dimnames( W )[[vT]] = O[["time_name"]]
              }
              if (length(O[["matchfrom"]]) > 2  ) {
                names(dimnames(W))[3] = vU
                dimnames( W )[[vU]] = O[["cyclic_name"]]
              }
              O[["predictions_exceedance"]] [["exceedance"]] [[ as.character(exceedance_threshold_predictions[b]) ]]= W
              W = NULL
            }
        }

        if (!is.null(deceedance_threshold_predictions)) {
      
            for ( b in 1:length(deceedance_threshold_predictions)) {
              m = apply ( predictions, 1, FUN=function(x) length( which(x < deceedance_threshold_predictions[b] ) ) ) / nposteriors
              W = reformat_to_array( input=m, matchfrom=O[["matchfrom"]],  matchto=matchto )
              names(dimnames(W))[1] = vS
              dimnames( W )[[vS]] = O[["space_name"]]
              m = NULL
              if ( length(O[["matchfrom"]]) > 1 ) {
                names(dimnames(W))[2] = vT
                dimnames( W )[[vT]] = O[["time_name"]]
              }
              if ( length(O[["matchfrom"]]) > 2 ) {
                names(dimnames(W))[3] = vU
                dimnames( W )[[vU]] = O[["cyclic_name"]]
              }
              O[["predictions_deceedance"]] [["deceedance"]] [[ as.character(deceedance_threshold_predictions[b]) ]]= W
              W = NULL
            }
        }
        
        predictions = NULL
      }

    }
    }
  }

  run_predict_end  = Sys.time()

  O[["ipred"]] = NULL
  O[["matchfrom"]] = NULL

  
  if (is.null(fn_res)) {
    message( "Saving results summary as a sublist in fit: fit$results : \n", fn_fit)
    message( "Return object is 'fit$results' (and not 'fit')")
    fit$results = O
    carstm_saveRDS( fit, file=fn_res, compress=compress, compression_level=compression_level  )
  } else {
    message( "Saving results summary as: \n", fn_res )
    carstm_saveRDS( O, file=fn_res, compress=compress , compression_level=compression_level  )
  }
 
  run_end  = Sys.time()
  
  if (be_verbose) {
    message("---------------------------------------")
    message("Run times:")
    message("Total: ", format(difftime(run_end, run_start)))
    message("Fit: ", format(difftime(run_end, run_fit_end)))
    message("Extract total: ", format(difftime(run_predict_end, run_extract_start)))
    message("Extract effects: ", format(difftime(run_predict_start, run_extract_start)))
    message("Extract predictions: ", format(difftime(run_predict_end, run_predict_start)))
    message("Remainder of time on file compression, save, etc: ", format(difftime(run_end, run_predict_end ) + difftime(run_fit_end, run_start) ))
    message("---------------------------------------")
  } 

  return(O)
  
}
