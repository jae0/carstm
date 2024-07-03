
carstm_model_inla = function(
  O, # p parameter list
  DS = "",
  sppoly =NULL,
  fit = NULL,
  space_id=NULL, time_id=NULL, cyclic_id=NULL, 
  fn_fit=tempfile(pattern="fit_", fileext=".RDS"), 
  fn_res=NULL,  
  redo_fit = TRUE,
  redo_posterior_simulations = TRUE,
  summarize_simulations=FALSE,
  theta=NULL,
  compress="qs-preset", 
  compression_level=1,
  qs_preset="high", 
  toget = c("modelinfo", "marginals", "random_spatial", "predictions"), 
  nposteriors=NULL, 
  posterior_simulations_to_retain=c( "predictions" ),
  exceedance_threshold=NULL, 
  deceedance_threshold=NULL, 
  exceedance_threshold_predictions=NULL,
  deceedance_threshold_predictions=NULL,
  debug=FALSE,
  ndiscretization = 1024L,
  eps = 1e-12,
  ... ) {
    
  if (DS=="modelled_fit") {
    if (!is.null(fn_fit)) {
      if (file.exists(fn_fit)) {
        if (grepl("\\.RDS$", fn_fit)) {
          fit = aegis::read_write_fast(fn_fit)
        } else {
          load( fn_fit )
        }
      }
      if (is.null(fit)) message("Modelled fit was not found.")
    }
    return( fit )
  }

  if (DS=="modelled_summary") {  # model.*modelled
    if (!is.null(fn_res)) {
      O = NULL
      if (file.exists(fn_res)) {
        if (grepl("\\.RDS$", fn_fit)) {
          res = aegis::read_write_fast(fn_res)
        } else {
          load( fn_res)
        }
      }
      if (is.null(O)) message("Summary of fit was not found.")
      return( O )
    } else {
      if (file.exists(fn_fit)){
        if (grepl("\\.RDS$", fn_fit)) {
          fit = aegis::read_write_fast(fn_fit)
        } else {
          load( fn_fit )
        }
      } 
      if (!is.null(fit)) {
        if (exists( "results", fit)) {
          return( fit$results )
        }
      }
      message("Modelled results not found.")
    }
  }
  

  if (DS=="carstm_modelinfo") {   
    fn_modelinfo = gsub( "_fit~", "_modelinfo~", fn_fit, fixed=TRUE )
    O = read_write_fast( file=fn_modelinfo )
    return( O )
  }

  if (DS=="carstm_marginals") {   
    fn_marginals = gsub( "_fit~", "_marginals~", fn_fit, fixed=TRUE )
    Omarginals = read_write_fast( file=fn_marginals )
    return( Omarginals )
  }

  if (DS=="carstm_randomeffects") {   
    fn_randomeffects = gsub( "_fit~", "_randomeffects~", fn_fit, fixed=TRUE )
    Orandom = read_write_fast( file=fn_randomeffects ) 
    return( Orandom )  
  }

  if (DS=="carstm_predictions") {   
    fn_preds = gsub( "_fit~", "_predictions~",  fn_fit, fixed=TRUE )
    Opredictions = read_write_fast( file=fn_preds ) 
    return( Opredictions)  
  }

  if (DS=="carstm_samples") {   
      fn_samples = gsub( "_fit~", "_samples~", fn_fit, fixed=TRUE )
      Osamples = read_write_fast(  file=fn_samples ) 
      return( Osamples )
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
  fmre = O[["fm"]]$random_effects

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
  # posterior samples is on link scale ((for both experimental and classical))
  # marginals.fitted.values are on response scale (for both experimental and classical)
  # summary.fitted.values on response scale ((for both experimental and classical))


  # --------------------------------------------
  # --------------------------------------------
  # --------------------------------------------

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
      if ( be_verbose) message( "Data is a function. Getting data ...")
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

    DM = unique( fmre$dimensionality )  
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

    if ( be_verbose ) message("Prepping data ")
      # dev.new()
      # hist( inla_args[["data"]][[vY]][ii] , main="Histogram of input variable to model" )
    

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

    if (!exists("control.inla", inla_args)) inla_args[["control.inla"]] = list( strategy='adaptive', cmin=0 ) #int.strategy='eb'
    if (!exists("control.predictor", inla_args)) inla_args[["control.predictor"]] = list( compute=TRUE, link=1  ) #everything on link scale
    if (!exists("control.mode", inla_args ) ) inla_args[["control.mode"]] = list( restart=TRUE ) 
    if (!is.null(theta) ) inla_args[["control.mode"]]$theta= theta
    
    if (!exists("control.compute", inla_args)) inla_args[["control.compute"]] = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE )


    # if (!exists("control.fixed", inla_args)) inla_args[["control.fixed"]] = list(mean.intercept=0.0, prec.intercept=0.001, mean=0, prec=0.001)
    if (!exists("control.fixed", inla_args)) inla_args[["control.fixed"]] = H$fixed

    setDF(inla_args[["data"]]) # in case .. INLA requires this ?
 
 
    fit = try( do.call( inla, inla_args ) )      

    if (inherits(fit, "try-error" )) {
      inla_args[["control.inla"]] = list( int.strategy='eb', cmin=0 )
      fit = try( do.call( inla, inla_args ) )      
    }

    if (inherits(fit, "try-error" )) {
      inla_args[["safe"]] = TRUE
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
 

    O[["carstm_prediction_surface_parameters"]] = NULL
    
    # fit$modelinfo = O  # store in case a restart is needed


    fn_modelinfo = gsub( "_fit~", "_modelinfo~", fn_fit, fixed=TRUE )
    read_write_fast( data=O, file=fn_modelinfo, compress=compress, compression_level=compression_level, qs_preset=qs_preset )
      
    message( "Model info saved as: \n", fn_modelinfo )

    gc()



    fit$.args = NULL

    inla_args= NULL; gc()
 
    read_write_fast( data=fit, file=fn_fit, compress=compress, compression_level=compression_level, qs_preset=qs_preset )
    if (be_verbose)  message( "\nModel fit saved as: \n", fn_fit )

  } # END redo fit

  run_fit_end  = Sys.time()
    
    
  ### END Modelling

  # --------------------------------------------
  # --------------------------------------------
  # --------------------------------------------


  if (is.null(fit)) {
    if (be_verbose)  message( "\nLoading fit as: redo_fit=FALSE ...\n", fn_fit )
    fit =aegis::read_write_fast( fn_fit )
  }
  
  if (is.null(fit)) {
    message( "Fit file not found, set option: redo_fit=TRUE" )
    stop()
  }
 
  # do the computations here as fit can be massive ... best not to copy, etc ..

  # --------------------------------------------
  # --------------------------------------------
  # --------------------------------------------
  fn_modelinfo = gsub( "_fit~", "_modelinfo~", fn_fit, fixed=TRUE )
  if (exists("modelinfo", fit)) {
    O = fit$modelinfo  # priority to the fit$modelinfo
  } else {
    if (file.exists(fn_modelinfo)) O = aegis::read_write_fast( fn_modelinfo )
  } 



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
      Z = try( apply_generic( Z, inla.tmarginal, fun=invlink, n=ndiscretization) )
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

  # some info required for posterior simulations
  O[["names.fixed"]] = fit$names.fixed
  O[["nposteriors"]] = nposteriors

  
  if (exists("debug")) if (is.character(debug)) if (debug=="summary") browser()

 

  if ( "summary" %in% toget) {
    
    if (be_verbose)  message("Extracting parameter summary" )
 
    dic = fit$dic[c("dic", "p.eff", "dic.sat", "mean.deviance")]
    waic = fit$waic[c("waic", "p.eff")]
    mlik = fit$mlik[2]

    O[["direct"]] = summary(fit)
    # print(O[["direct"]])

    # remove a few unused but dense data objects
    O[["direct"]]$linear.predictor = NULL
    O[["direct"]]$waic$local.waic = NULL
    O[["direct"]]$waic$local.p.eff = NULL
    O[["direct"]]$dic$local.dic = NULL
    O[["direct"]]$dic$local.p.eff = NULL
    O[["direct"]]$dic$local.dic.sat = NULL
    O[["direct"]]$dic$family = NULL
     
  }

  Omarginals = list()

  if ( "marginals" %in% toget) {

     # extract and back transform where possible
    if (exists( "marginals.fixed", fit)) {
      m = fit$marginals.fixed  # make a copy to do transformations upon
      fi = grep("Intercept", names(m) )
      m = try( apply_generic( m, function(x) marginal_clean( inla.tmarginal( invlink, x, n=ndiscretization)) ), silent=TRUE  )
      if ( exists("data_transformation", O))  {
        m[[fi]] = inla.tmarginal( O$data_transformation$backward, m[[fi]] , n=ndiscretization ) # on user scale
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
      Omarginals[["fixed_effects"]] = W
      m = NULL
      W = NULL
      gc()
    }


    if (exists( "marginals.hyperpar", fit)) {
      
      # hyperpar (variance components)
      hyps =  row.names(fit$summary.hyperpar)

      prcs = grep( "^Precision.*", hyps, value=TRUE )
      if (length(prcs) > 0) {
  
        m = fit$marginals.hyperpar[prcs]
        m = try( apply_generic( m, inla.tmarginal, fun=function(y) 1/sqrt_safe( y, eps ), n = ndiscretization ), silent=TRUE)
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
        
        Omarginals[["random_effects"]] = m[, tokeep, drop =FALSE] 
        m = NULL
        gc()
        
      }
 

      # update phi's, lambda's (used in besagproper2 -- Leroux model) .. etc
      rhos = grep( "^Rho.*|^GroupRho.*", hyps, value=TRUE )
      phis = grep( "^Phi.*", hyps, value=TRUE )
      other = grep( "^Lambda.*|^Diagonal.*|zero-probability.*", hyps, value=TRUE )

      known = c( rhos, phis, other )
      unknown = setdiff( hyps, c(prcs, known) )

      hyps = prcs = other = known = NULL

      if (length(rhos) > 0) {
        m = marginal_summary( fit$marginals.hyperpar[ rhos ] )
        if (test_for_error(m) =="error") {  
          m = fit$summary.hyperpar[rhos, 1:5]
          colnames(m) = tokeep
        }
        Omarginals[["random_effects"]] = rbind( Omarginals[["random_effects"]],  m[, tokeep, drop =FALSE] )
        m = NULL
      }
      rhos = NULL

      if (length(phis) > 0) {
        m = marginal_summary( fit$marginals.hyperpar[ phis ] )
        if (test_for_error(m) =="error") {  
          m = fit$summary.hyperpar[phis, 1:5]
          colnames(m) = tokeep
        }
        Omarginals[["random_effects"]] = rbind( Omarginals[["random_effects"]], m[, tokeep, drop =FALSE] )
        m = NULL
      }
      phis = NULL

      if (length(unknown) > 0) {
        m = marginal_summary( fit$marginals.hyperpar[ unknown ] )
        if (is.character(m) && m=="error")  {
          #  alternatively: m[,"mode"] = apply_simplify( fit$marginals.hyperpar[ unknown ], FUN=function(x) inla.mmarginal( x ))
          m = fit$summary.hyperpar[unknown, 1:5]
          colnames(m) = tokeep
        }
        Omarginals[["random_effects"]] = rbind( Omarginals[["random_effects"]], m[, tokeep, drop =FALSE])
        m = NULL
      }
      unknown = NULL


      Omarginals[["random_effects"]] = list_to_dataframe(  Omarginals[["random_effects"]] )
      Omarginals[["random_effects"]]$ID = row.names( Omarginals[["random_effects"]] )
      
      gc()

      if (length(fit[["marginals.random"]]) > 0) { 

        if (exists("debug")) if (is.character(debug)) if ( debug =="random_covariates") browser()
 
        raneff = setdiff( names( fit$marginals.random ), c(vS, vS2, vS3  ) )
        for (rnef in raneff) {
          m = marginal_summary( fit$marginals.random[[rnef]], invlink=invlink )
          if (test_for_error(m) =="error") {  
             message( "failed to transform marginals .. copying directly from INLA summary instead: ", rnef)
            m = fit$summary.random[[rnef]][, inla_tokeep, drop =FALSE ]
            names(m) =  tokeep
            Omarginals[[rnef]] = m
            Omarginals[[rnef]]$ID = fit$summary.random[[rnef]]$ID
            Omarginals[[rnef]] = list_to_dataframe( Omarginals[[rnef]] )
          } else {
            Omarginals[[rnef]] = m[, tokeep, drop =FALSE]
            Omarginals[[rnef]]$ID = fit$summary.random[[rnef]]$ID
            Omarginals[[rnef]] = list_to_dataframe( Omarginals[[rnef]] )
          }
        }
        m = raneff = NULL
        gc()
      }
    }
 

    if (be_verbose)  {
      # message( "")
      # message( "Random effects:")
      # print(  Omarginals[["random_effects"]]  )   
      # message( "\n--- NOTE --- 'SD *' from marginal summaries are on link scale")
      # message( "--- NOTE --- SD * from posteriors simulations are on user scale")
      # message( "")
    }

    # save summary
    fn_marginals = gsub( "_fit~", "_marginals~", fn_fit, fixed=TRUE )
    read_write_fast( data=Omarginals, file=fn_marginals, compress=compress, compression_level=compression_level, qs_preset=qs_preset )
    Omarginals = NULL 
    gc()
  }  # end parameters
   
  
  
  Orandom = list()  # random effects (random spatial and random spatiotemporal) from INLA's marginals

  # separate out random spatial and randomm spatiotemporal (as they can be large arrays)
  if ("random_spatial" %in% toget) {
    # space only
    Z = NULL
    iSP = which( fmre$dimensionality=="s" & fmre$level=="main")
    
    if (length(iSP) > 0 ) {

      if (be_verbose)  message("Extracting random spatial effects"  )
      if (exists("debug")) if (is.character(debug)) if ( debug =="random_spatial") browser()

      matchto = list( space=O[["space_id"]] )

      W = array( NA, dim=c( O[["space_n"]], length(tokeep) ), dimnames=list( space=O[["space_name"]], stat=tokeep ) )
      names(dimnames(W))[1] = vS  # need to do this in a separate step ..

      Orandom[[vS]] = list()  # space as a main effect  vS==vS

      if (length(iSP) == 1) {

        # vS = fmre$vn[ iSP ]  # == vS as it is a single spatial effect
        
        model_name = fmre$model[ iSP ]  # iid (re_unstructured)

        m = fit$marginals.random[[vS]]
 
        m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=ndiscretization ) , silent=TRUE )
        m = try( apply_generic( m, inla.zmarginal, silent=TRUE ), silent=TRUE )
        m = try( simplify2array( m ), silent=TRUE)
        m = try( list_simplify( m ), silent=TRUE )
        # single spatial effect (eg in conjuction with besag) .. indexing not needed but here in case more complex models ..
        if (test_for_error(m) =="error") {  
          message( "failed to transform random_spatial marginals .. copying directly from INLA summary instead")
          m = fit$summary.random[[vS]][, inla_tokeep ]
          names(m) =  tokeep
        } 
        
        if ( model_name %in% c("bym", "bym2") ) {
          # bym2 effect is coded by INLA as a double length vector: re_total and re_neighbourhood  
          Z = expand.grid( space=O[["space_id"]], type = c("re_total", "re_neighbourhood"), stringsAsFactors =FALSE )

          #  extract re_total main effects
          ire = which(Z$type=="re_total")
          matchfrom = list( space=Z[["space"]][ire] )

          for (k in 1:length(tokeep)) {
            W[,k] = reformat_to_array( input = unlist(m[ire, tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
          }
 
          Orandom[[vS]] [["re_total"]] = W [, tokeep, drop =FALSE]

        } else {
          # single spatial effect that is not bym2
          # this is redundant with iSP being a single factor, but approach is generalizable for higher dims 
          Z = expand.grid( space=O[["space_id"]], type="re_neighbourhood", stringsAsFactors =FALSE )
        }

        ine =  which(Z$type=="re_neighbourhood")
        matchfrom = list( space=Z[["space"]][ine] )
        W[] = NA
        for (k in 1:length(tokeep)) {
          W[,k] = reformat_to_array( input = unlist(m[ine, tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
        }
        Orandom[[vS]] [["re_neighbourhood"]] = W [, tokeep, drop =FALSE]

      }

      if (length(iSP) == 2) {
        
        for (j in 1:length(iSP)) {
          
          vnS = fmre$vn[ iSP[j] ]
          model_name = fmre$model[ iSP[j] ]  
          
          m = fit$marginals.random[[vnS]]
         
          m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=ndiscretization), silent=TRUE  )
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
          Orandom[[vnS]] [[model_name]] = data.frame( W [, tokeep, drop =FALSE], ID=row.names(W) )

        }
      }

      Z = m = matchfrom = NULL
      gc()

    }
  }  # end random spatial effects

  matchfrom = i1 = i2= NULL
  Z = W = m = space = space1 = space2 = skk1 = skk2 = iSP = NULL
  gc()

  if ("random_spatiotemporal" %in% toget ) {
    # space-time

    iST = which( fmre$dimensionality=="st" & fmre$level=="main")
    
    if (length(iST) > 0 ) {

      if (be_verbose)  message("Extracting random spatiotemporal effects"  )

      if (exists("debug")) if (is.character(debug)) if ( debug =="random_spatiotemporal") browser()

      matchto = list( space=O[["space_id"]], time=O[["time_id"]]  )
      
      W = array( NA, 
        dim=c( O[["space_n"]], O[["time_n"]], length(tokeep) ), 
        dimnames=list( space=O[["space_name"]], time=O[["time_name"]], stat=tokeep ) )
      names(dimnames(W))[1] = vS  # need to do this in a separate step ..
      names(dimnames(W))[2] = vT  # need to do this in a separate step ..

      if (length(iST) == 1) {

        vnST = fmre$vn[ iST ]
        model_name = fmre$model[ iST ]   

        if (exists(vnST, fit$marginals.random )) {
 
          m = fit$marginals.random[[vnST]]
          m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=ndiscretization) , silent=TRUE )

          m = try( apply_generic( m, inla.zmarginal, silent=TRUE), silent=TRUE )
          m = try( simplify2array( m ), silent=TRUE)
          m = try( list_simplify( m ), silent=TRUE )
          if (test_for_error(m) =="error") {  
            message( "!!! Failed to transform random_spatial marginals .. copying directly from INLA summary instead")
            m = fit$summary.random[[vnST]][, inla_tokeep ]
            names(m) = tokeep
          } 
        
          if ( model_name %in% c("bym", "bym2")  ) {
            # bym2 effect: re_total and re_neighbourhood with annual results
            Z = expand.grid( space=O[["space_id"]], type = c("re_total", "re_neighbourhood"), time=O[["time_id"]], stringsAsFactors =FALSE )

            #  spatiotemporal interaction effects 
            ire = which(Z$type=="re_total")
            matchfrom = list( space=Z[["space"]][ire], time=Z[["time"]][ire]  )

            for (k in 1:length(tokeep)) {
              W[,,k] = reformat_to_array(  input = unlist(m[ire, tokeep[k]]), matchfrom = matchfrom, matchto = matchto )
            }
            Orandom[[vnST]] [["re_total"]] =  W [,, tokeep, drop =FALSE] 

          } else {
            # besag effect: with annual results
            Z = expand.grid( space=O[["space_id"]], type ="re_neighbourhood", time=O[["time_id"]], stringsAsFactors =FALSE )
          }

          #  spatiotemporal interaction effects 
          ine =  which(Z$type=="re_neighbourhood")
          matchfrom = list( space=Z[["space"]][ine], time=Z[["time"]][ine]  )
            
          for (k in 1:length(tokeep)) {
            W[,,k] = reformat_to_array( input = unlist(m[ine,tokeep[k]]), matchfrom=matchfrom, matchto=matchto )
          }
          Orandom[[vnST]] [["re_neighbourhood"]] =   W [,, tokeep, drop =FALSE] 

        }
      }
      
      if (length(iST) == 2) {

        for (j in 1:length(iST)) {
          vnST = fmre$vn[ iST[j] ]
          model_name = fmre$model[ iST[j] ]  

          m = fit$marginals.random[[vnST]]
          m = try( apply_generic( m, inla.tmarginal, fun=invlink, n=ndiscretization), silent=TRUE )
          m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE )
          m = try( simplify2array( m ), silent=TRUE)
          m = try( list_simplify( m ), silent=TRUE )
          if (test_for_error(m) =="error") { 
            message( "!!! Failed to transform random_spatiotemporal marginals .. copying directly from INLA summary instead")
            m = fit$summary.random[[vnST]][, inla_tokeep ]
            names(m) = tokeep
          } 
          
          Z = expand.grid( space=O[["space_id"]], type =model_name, time=O[["time_id"]], stringsAsFactors =FALSE )
          jst =  which(Z$type==model_name)
          matchfrom = list( space=Z[["space"]][jst], time=Z[["time"]][jst]  )
          for (k in 1:length(tokeep)) {
            W[,,k] = reformat_to_array( input = unlist(m[jst, tokeep[k] ]), matchfrom = matchfrom, matchto = matchto  )
          }
          Orandom[[vnST]] [[model_name]] = W [,, tokeep, drop =FALSE]
          m = NULL
        }
      }
      Z = W = m = space_time = space_time1 = space_time2 = skk1 = skk2 = NULL
  
    }


  }  # end random spatio-temporal effects
   
  # save random effects to file
  if (length(Orandom) > 0) {
    fn_randomeffects = gsub( "_fit~", "_randomeffects~", fn_fit, fixed=TRUE )
    read_write_fast( data=Orandom, file=fn_randomeffects, compress=compress, compression_level=compression_level, qs_preset=qs_preset ) 
    Orandom = NULL 
    gc()
  }


  Opredictions = list()
  
  if ("predictions" %in% toget ) {

    # marginals.fitted.values are on response scale (for both experimental and classical) and  already incorporates offsets
    # summary.fitted.values on response scale ((for both experimental and classical)) 
    # posteriors are on link scale and already incorporates offsets

    if (be_verbose) message("Extracting posterior predictions" )

    if (exists("debug")) if (is.character(debug)) if ( debug =="predictions") browser()
 
  
    if (exists("marginals.fitted.values", fit)) {

      if (length(fit[["marginals.fitted.values"]]) > 0 ) {

        if ( O[["dimensionality"]] == "space" ) {

          m = fit$marginals.fitted.values[O[["ipred"]]]  # already incorporates offsets
  
          if ( exists("data_transformation", O))  m = apply_generic( m, backtransform )

          m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
          m = try( list_simplify( simplify2array( m ) ), silent=TRUE)
          if (test_for_error(m) =="error") {  
              message( "!!! Failed to summarize marginals.fitted.values .. copying directly from INLA summary instead")
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
          Opredictions = W[, tokeep, drop =FALSE]
          m = W = NULL
        }


        if (O[["dimensionality"]] == "space-time"  ) {

          m = fit$marginals.fitted.values[O[["ipred"]]]   
  
          if (exists("data_transformation", O)) m = apply_generic( m, backtransform )
          m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
          m = try( list_simplify( simplify2array( m ) ), silent=TRUE)
          if (test_for_error(m) =="error") {  
              message( "!!! Failed to summarize marginals.fitted.values .. copying directly from INLA summary instead")
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
          Opredictions = W[,, tokeep, drop =FALSE]
      
          m = W = NULL
        }


        if ( O[["dimensionality"]] == "space-time-cyclic" ) {

          m = fit$marginals.fitted.values[O[["ipred"]]]   

          if (exists("data_transformation", O)) m = apply_generic( m, backtransform )

          m = try( apply_generic( m, inla.zmarginal, silent=TRUE  ), silent=TRUE)
          m = try( list_simplify( simplify2array( m ) ), silent=TRUE)
          if (test_for_error(m) =="error") {  
              message( "!!! Failed to summarize marginals.fitted.values .. copying directly from INLA summary instead")
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
          Opredictions = W[,,, tokeep, drop =FALSE]
          m = W = NULL
        }
      }
    }

    # save summary
    fn_preds = gsub( "_fit~", "_predictions~",  fn_fit, fixed=TRUE )
    read_write_fast( data=Opredictions, file=fn_preds, compress=compress, compression_level=compression_level, qs_preset=qs_preset ) 
    Opredictions = NULL 
    gc()
     
  }
  
  O[["ipred"]] = NULL
 

  # --------------------------------------------
  # --------------------------------------------
  # --------------------------------------------


  run_post_samples  = Sys.time()

  if (redo_posterior_simulations) {

    if (exists("debug")) if (is.character(debug)) if ( debug =="posterior_samples") browser()

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
      S = try( inla.posterior.sample( nposteriors, fit, add.names=FALSE, use.improved.mean=FALSE ), silent=TRUE)
    }
    if ( "try-error" %in%  class(S) ) {
      S = try( inla.posterior.sample( nposteriors, fit, add.names=FALSE, skew.corr =FALSE ), silent=TRUE)
    }
    if ( "try-error" %in%  class(S) ) {
      S = try( inla.posterior.sample( nposteriors, fit, add.names=FALSE, use.improved.mean=FALSE, skew.corr =FALSE ), silent=TRUE)
    }
    
    if ( "try-error" %in%  class(S) ) stop("posterior sampling error")

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
    
    if (be_verbose)  message( "\nPosterior simulations complete. Extracting required components ... "  )


    fit = NULL; gc()  # no longer needed .. remove from memory

    Osamples = list()
    
    if (summarize_simulations) Osamples[["summary"]] = list()

    # extraction here as no need for fit object anymore (reduces RAM requirements)
    if (is.null(exceedance_threshold)) if (exists("exceedance_threshold", O)) exceedance_threshold = O[["exceedance_threshold"]]
    if (is.null(deceedance_threshold)) if (exists("deceedance_threshold", O)) deceedance_threshold = O[["deceedance_threshold"]]

    if (is.null(deceedance_threshold_predictions)) if (exists("deceedance_threshold_predictions", O)) deceedance_threshold_predictions = O[["deceedance_threshold_predictions"]]
    if (is.null(exceedance_threshold_predictions)) if (exists("exceedance_threshold_predictions", O)) exceedance_threshold_predictions = O[["exceedance_threshold_predictions"]]

    for (z in c("tag", "start", "length") ) assign(z, attributes(S)[[".contents"]][[z]] )  # index info

    if ( "summary" %in% posterior_simulations_to_retain ) {
      
      if (exists("debug")) if (is.character(debug)) if ( debug =="posterior_samples_summary") browser()
      # -- check variable names here
      # posteriors
      flabels= O[["names.fixed"]]
      nfixed = length(O[["names.fixed"]])
      Osamples[["fixed_effects"]] = array(NA, dim=c( nfixed, O[["nposteriors"]]  ) )
      fkk = inla_get_indices(O[["names.fixed"]], tag=tag, start=start, len=length, model="direct_match")
      fkk = unlist(fkk)
      for (i in 1:O[["nposteriors"]]) Osamples[["fixed_effects"]][,i] = S[[i]]$latent[fkk,]  
      Osamples[["fixed_effects"]] = invlink(Osamples[["fixed_effects"]])
      row.names(Osamples[["fixed_effects"]]) = flabels

      if (summarize_simulations) {
        Osamples[["summary"]][["fixed_effects"]] = posterior_summary(format_results( fixed, labels=flabels))
      }
      fkk = fixed = NULL
 

      other_random = setdiff( names(O$random), c(vS, vS2 ) )

      if (length(other_random) > 0 ) {
        nrandom = sum(sapply(O$random[other_random], nrow))  
        Osamples[["random_effects"]] = array(NA, dim=c(nrandom, O[["nposteriors"]]  ) )
        rkk = inla_get_indices(other_random, tag=tag, start=start, len=length, model="direct_match"  )  # if bym2, must be decomposed  
        rkk = unlist( rkk )
        for (i in 1:O[["nposteriors"]]) random[,i] = S[[i]]$latent[rkk,]  
        rlabels = names(S[[1]]$latent[rkk,])
        Osamples[["random_effects"]] = invlink(Osamples[["random_effects"]])
        row.names(Osamples[["random_effects"]]) = rlabels

        if (summarize_simulations) {
          Osamples[["summary"]][["random_effects"]] = posterior_summary( format_results( Osamples[["random_effects"]], labels=rlabels ) )
        }
        rkk = nrandom= rlabels= NULL
      }

      hyper_names = names(S[[1]]$hyperpar)
      nhyperpar = length(hyper_names)
      if (nhyperpar == 0) stop("No hyperparameters? ... control.mode(restart=TRUE) might help" )
      Osamples[["hyperpars"]] = array( NA,dim=c(nhyperpar, O[["nposteriors"]]  ) )
      for (i in 1:O[["nposteriors"]]) Osamples[["hyperpars"]][,i] = S[[i]]$hyperpar
      k = grep("Precision", hyper_names)
      if (length(k) > 0) {
          Osamples[["hyperpars"]][k,] = 1 / sqrt(Osamples[["hyperpars"]][k,]) 
          hyper_names[k] = gsub("Precision for", "SD", hyper_names[k] )
      }
      hyper_names = gsub("for ", "", hyper_names )
      hyper_names = gsub("the ", "", hyper_names )
      row.names( Osamples[["hyperpars"]] ) = hyper_names
 
      if (summarize_simulations) {
        Osamples[["summary"]][["hyperpars"]] = posterior_summary( format_results( Osamples[["hyperpars"]], labels=hyper_names) )
      }
      hyper_names = nhyperpar = k = NULL
    }  # end summary
    
    k = known = unknown = m  = NULL

    gc()
  
    if ( "random_spatial" %in% posterior_simulations_to_retain ) {
      if (exists("debug")) if (is.character(debug)) if ( debug =="posterior_samples_spatial") browser()
      # -- check variable names here

      # POSTERIOR samples
      space = array(NA, dim=c( O$space_n, O[["nposteriors"]]  ) )
      row.names(space) = O[["space_name"]]
      space1 = space2 = NULL

      fmre =  O[["fm"]]$random_effects
      iSP = which( fmre$dimensionality=="s" & fmre$level=="main")

      matchto = list( space=O[["space_id"]] )

      if (length(iSP) == 1) {
        stx1 = paste("^", fmre$vn[iSP], "$", sep="")
        if (fmre$model[iSP] %in% c("bym2", "bym") ) {
          # special case bym2 has two factors rolled together
          space1 = array(NA, dim=c( O$space_n, O[["nposteriors"]]  ) )
          space2 = array(NA, dim=c( O$space_n, O[["nposteriors"]]  ) )
          row.names(space1) = O[["space_name"]]
          row.names(space2) = O[["space_name"]]
          skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length, model="bym2" )  # if bym2, must be decomposed  
          for (i in 1:O[["nposteriors"]]) {
            space[,i]  = S[[i]]$latent[skk1[["re_total"]],]   # ICAR + IID 
            space2[,i] = S[[i]]$latent[skk1[["re_neighbourhood"]],]   # ICAR
          }      
          space1 = space - space2  # re_unstructured (IID) 
        } else {
          # single spatial effect of some kind
          skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  
          skk1 = unlist(skk1)
          for (i in 1:O[["nposteriors"]]) {
            space[,i] = S[[i]]$latent[skk1,]   # ICAR or IID 
          }
        }
      }

      if (length(iSP) == 2) {
        # Assume additive 
        space1 = array(NA, dim=c( O$space_n, O[["nposteriors"]]  ) )
        space2 = array(NA, dim=c( O$space_n, O[["nposteriors"]]  ) )
        row.names(space1) = O[["space_name"]]
        row.names(space2) = O[["space_name"]]
        model_name1 = fmre$model[ iSP[1] ] 
        model_name2 = fmre$model[ iSP[2] ] 
        stx1 = paste("^", fmre$vn[iSP[1]], "$", sep="")
        stx2 = paste("^", fmre$vn[iSP[i]], "$", sep="")
        skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
        skk2 = inla_get_indices(stx2, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
        for (i in 1:O[["nposteriors"]]) {
          space1[,i] = S[[i]]$latent[skk1[[1]],]  
          space2[,i] = S[[i]]$latent[skk2[[1]],]
        }      
        space = space1 + space2  
        message("Two spatial effects .. assuming they are are additive")
      }

      space = invlink(space)   # the overall spatial random effect 
      row.names(space) = O[["space_name"]]
      matchfrom = list( space=O[["space_id"]] )
      
      Osamples[[vS]] = list()

      if (!is.null(exceedance_threshold)) {
        if (be_verbose)  message("Extracting random spatial effects exceedence"  )
        for ( b in 1:length(exceedance_threshold)) {
          m = apply ( space, 1, FUN=function(x) length( which(x > exceedance_threshold[b]) ) ) / O[["nposteriors"]]
          m = reformat_to_array( input = m, matchfrom=matchfrom, matchto = matchto )
          names(dimnames(m))[1] = vS
          dimnames( m )[[vS]] = O[["space_name"]]
          Osamples[[vS]] [["exceedance"]] [[as.character(exceedance_threshold[b])]] = data.frame( m [, tokeep, drop =FALSE], ID=row.names(m) )
          m = NULL
        }
      }

      if (!is.null(deceedance_threshold)) {
        if (be_verbose)  message("Extracting random spatial effects deceedance"  )
        # redundant but generalizable to higher dims
        for ( b in 1:length(deceedance_threshold)) {
          m = apply ( space, 1, FUN=function(x) length( which(x < deceedance_threshold[b]) ) ) / O[["nposteriors"]]
          m = reformat_to_array( input = m, matchfrom=matchfrom, matchto=matchto  )
          names(dimnames(m))[1] = vS
          dimnames( m )[[vS]] = O[["space_name"]]
          Osamples[[vS]] [["deceedance"]] [[as.character(deceedance_threshold[b])]] = data.frame( m [, tokeep, drop =FALSE], ID=row.names(m) ) 
          m = NULL
        }
      }

      if (summarize_simulations) {
        Osamples[["summary"]] [[vS]] = list()
        m = posterior_summary( format_results( space, labels=O[["space_name"]]  ) )
        W = array( NA, dim=c( O[["space_n"]], length(tokeep) ), dimnames=list( space=O[["space_name"]], stat=tokeep ) )
        names(dimnames(W))[1] = vS  # need to do this in a separate step ..
        for (k in 1:length(tokeep)) {
          W[,k] = reformat_to_array(  input = m[, tokeep[k]], matchfrom=matchfrom, matchto=matchto )
        }
        Osamples[["summary"]] [[vS]] [["re_total"]] = W[, tokeep, drop =FALSE] 
      }

      Osamples[[vS]][["re_total"]] = space  # already inverse link scale
      
      if ( "random_spatial12" %in% posterior_simulations_to_retain ) {
        if (!is.null(space1)) Osamples [[vS]] [[model_name1]] =  invlink(space1) 
        if (!is.null(space2)) Osamples [[vS]] [[model_name2]] =  invlink(space2) 
      }
    } # END random spatial post samples

    matchfrom = i1 = i2= NULL
    Z = W = m = space = space1 = space2 = skk1 = skk2 = iSP = NULL
    gc()


    if ( "random_spatiotemporal" %in% posterior_simulations_to_retain ) {
      # posterior simulations
      if (exists("debug")) if (is.character(debug)) if ( debug =="posterior_samples_spatiotemporal") browser()

      fmre =  O[["fm"]]$random_effects
      iST = which( fmre$dimensionality=="st" & fmre$level=="main")

      matchto = list( space=O[["space_id"]], time=O[["time_id"]]  )

      space_time1 = space_time2  = space_time = array( NA, 
        dim=c( O[["space_n"]] * O[["time_n"]] , O[["nposteriors"]]  ) )
      L = CJ( time=O[["time_id"]], space=O[["space_id"]] )  # note:: CJ has reverse order vs expand.grid
      stlabels = paste(L[["space"]], L[["time"]], sep="_")

      if (length(iST) == 1) {
        stx1 = paste("^", fmre$vn[iST], "$", sep="")
        if (fmre$model[iST] %in% c("bym", "bym2") ) {
          # special case bym2 has two factors rolled together
          skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length, model="bym2" )  # if bym2, must be decomposed  
          for (i in 1:O[["nposteriors"]]) {
            space_time[,i]  = S[[i]]$latent[skk1[["re_total"]],]  
            space_time2[,i] = S[[i]]$latent[skk1[["re_neighbourhood"]],]
          }    
          space_time1 = space_time - space_time2   # re_unstructured
          row.names(space_time1) = stlabels
          row.names(space_time2) = stlabels
          row.names(space_time) = stlabels

        } else {
          # single spatial effect of some kind
          skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
          skk1 = unlist(skk1)
          for (i in 1:O[["nposteriors"]]) {
            space_time[,i] = S[[i]]$latent[skk1,] 
          }      
        }
      }
      
      if (length(iST) == 2) {
        model_name1 = fmre$model[ iST[1] ] 
        model_name2 = fmre$model[ iST[2] ] 
        stx1 = paste("^", fmre$vn[iST[1]], "$", sep="")
        stx2 = paste("^", fmre$vn[iST[2]], "$", sep="")
        skk1 = inla_get_indices(stx1, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
        skk2 = inla_get_indices(stx2, tag=tag, start=start, len=length )  # if bym2, must be decomposed  
        for (i in 1:O[["nposteriors"]]) {
          space_time1[,i] = S[[i]]$latent[skk1[[1]],] 
          space_time2[,i] = S[[i]]$latent[skk2[[1]],]
        }      
        space_time  = space_time1 + space_time2  # assume additive
        row.names(space_time1) = stlabels
        row.names(space_time2) = stlabels 
      }
    
      space_time = invlink(space_time)
      row.names(space_time) = stlabels

      Z = expand.grid( space=O[["space_id"]], type ="re_total", time=O[["time_id"]], stringsAsFactors =FALSE )
      jst =  which(Z$type=="re_total")
      matchfrom = list( space=Z[["space"]][jst], time=Z[["time"]][jst]  )

      Osamples[[vnST]] = list()

      if (!is.null(exceedance_threshold)) {
        if (be_verbose)  message("Extracting random spatiotemporal errors exceedence"  )

        Osamples[[vnST]] [["exceedance"]] = list()

        for ( b in 1:length(exceedance_threshold)) {
          m = apply ( space_time, 1, FUN=function(x) length( which(x > exceedance_threshold[b] ) ) ) / O[["nposteriors"]]
          m = reformat_to_array( input=m, matchfrom=matchfrom,  matchto=matchto )
          names(dimnames(m))[1] = vS
          names(dimnames(m))[2] = vT
          dimnames( m )[[vS]] = O[["space_name"]]
          dimnames( m )[[vT]] = O[["time_name"]]
          Osamples[[vnST]] [["exceedance"]] [[as.character(exceedance_threshold[b])]] = m
        }
        m = NULL
      }

      if (!is.null(deceedance_threshold)) {
        if (be_verbose)  message("Extracting random spatiotemporal errors deceedance"  )

        Osamples[[vnST]] [["deceedance"]] = list()

        for ( b in 1:length(deceedance_threshold)) {
          m = apply ( space_time, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / O[["nposteriors"]]
          m = reformat_to_array( input = m, matchfrom=matchfrom,  matchto=matchto )
          names(dimnames(m))[1] = vS
          names(dimnames(m))[2] = O[["time_name"]]
          dimnames( m )[[vS]] = O[["space_name"]]
          dimnames( m )[[vT]] = vT
          Osamples[[vnST]] [["deceedance"]] [[as.character(deceedance_threshold[b])]] = m
        }
        m = NULL
      }

      if (summarize_simulations) {
        W = array( NA, 
          dim=c( O[["space_n"]], O[["time_n"]], length(tokeep) ), 
          dimnames=list( space=O[["space_name"]], time=O[["time_name"]], stat=tokeep ) )
        names(dimnames(W))[1] = vS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vT  # need to do this in a separate step ..

        m = posterior_summary(format_results( space_time, labels=stlabels ))
        for (k in 1:length(tokeep)) {
          W[,,k] = reformat_to_array(  input = m[, tokeep[k]], matchfrom=matchfrom, matchto=matchto )
        }
        Osamples[["summary"]][[vnST]] = list()
        Osamples[["summary"]][[vnST]] [["re_total"]] = W[,, tokeep, drop =FALSE] 
      }

      if ( "random_spatiotemporal" %in% posterior_simulations_to_retain ) {
        Osamples [[vnST]] [["re_total"]] = space_time
      }
      if ( "random_spatiotemporal12" %in% posterior_simulations_to_retain ) {
        if (!is.null(space_time1)) Osamples [[vnST]] [[model_name1]] = invlink(space_time1)
        if (!is.null(space_time2)) Osamples [[vnST]] [[model_name2]] = invlink(space_time2)
      }

    }  # END sp-temp post samples
    Z = W = m = space_time = space_time1 = space_time2 = skk1 = skk2 = NULL
    gc()


    if ( "predictions" %in% posterior_simulations_to_retain ) {
  
      # marginals.fitted.values are on response scale (for both experimental and classical) and  already incorporates offsets
      # summary.fitted.values on response scale ((for both experimental and classical)) 
      # posteriors are on link scale and already incorporates offsets
      if (exists("debug")) if (is.character(debug)) if ( debug =="posterior_samples_predictions") browser()
      if (be_verbose) message("Extracting posterior predictions" )

      # prepapre prediction simulations (from joint posteriors)
      ptx = "^Predictor"
      preds = O[["predictions_range"]][1]:O[["predictions_range"]][2]  
      npredictions = length(preds)
      
      pkk = inla_get_indices(ptx, tag=tag, start=start, len=length)
      pkk = unlist(pkk)[preds]
      predictions = array(NA, dim=c( npredictions, O[["nposteriors"]]  ) )
      for (i in 1:O[["nposteriors"]])  predictions[,i] = S[[i]]$latent[pkk, ]
      predictions = invlink(predictions)  #  required for both classical and experimental

      pkk = preds = S = NULL
      gc()

      if ( exists("data_transformation", O))  predictions = O$data_transformation$backward( predictions  )
  
      if ( O[["dimensionality"]] == "space" ) {

          matchto = list( space=O[["space_id"]] )

          W = array( NA, 
            dim=c( O[["space_n"]], O[["nposteriors"]] ),  
            dimnames=list( space=O[["space_name"]], sim=1:O[["nposteriors"]] ) )
          names(dimnames(W))[1] = vS  # need to do this in a separate step ..
          for (k in 1:O[["nposteriors"]] ) {
            W[,k] = reformat_to_array( input=unlist(predictions[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
          }
          Osamples[["predictions"]] = W[, drop =FALSE]
          W = NULL 
      }


      if (O[["dimensionality"]] == "space-time"  ) {

          matchto = list( space=O[["space_id"]], time=O[["time_id"]] )

          W = array( NA, 
            dim=c( O[["space_n"]], O[["time_n"]], O[["nposteriors"]] ),  
            dimnames=list( space=O[["space_name"]], time=O[["time_name"]], sim=1:O[["nposteriors"]] ) )
          names(dimnames(W))[1] = vS  # need to do this in a separate step ..
          names(dimnames(W))[2] = vT  # need to do this in a separate step ..

          for (k in 1:O[["nposteriors"]] ) {
            W[,,k] = reformat_to_array( input=unlist(predictions[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
          }
          Osamples[["predictions"]] = W[,,, drop =FALSE]
          W = NULL
      }


      if ( O[["dimensionality"]] == "space-time-cyclic" ) {

          matchto = list( space=O[["space_id"]], time=O[["time_id"]], cyclic=O[["cyclic_id"]] )

          W = array( NA, 
            dim=c( O[["space_n"]], O[["time_n"]], O[["cyclic_n"]], O[["nposteriors"]] ),  
            dimnames=list( space=O[["space_name"]], time=O[["time_name"]], cyclic=O[["cyclic_name"]], sim=1:O[["nposteriors"]] ) )
          names(dimnames(W))[1] = vS  # need to do this in a separate step ..
          names(dimnames(W))[2] = vT  # need to do this in a separate step ..
          names(dimnames(W))[3] = vU  # need to do this in a separate step ..
          for (k in 1:O[["nposteriors"]] ) {
            W[,,,k] = reformat_to_array( input=unlist(predictions[,k]), matchfrom=O[["matchfrom"]], matchto=matchto )
          }
          Osamples[["predictions"]] = W[,,, ,drop =FALSE]
          W = NULL
      }
  
      
      if (!is.null(exceedance_threshold_predictions)) {
        if (be_verbose)  message("Extracting de/exceedance of predictions"  )
        
        Osamples[["predictions_exceedance"]] = list()
        Osamples[["predictions_exceedance"]] [["exceedance"]]  = list()

        for ( b in 1:length(exceedance_threshold_predictions)) {
            m = apply ( predictions, 1, FUN=function(x) length( which(x > exceedance_threshold_predictions[b] ) ) ) / O[["nposteriors"]]
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
            Osamples[["predictions_exceedance"]] [["exceedance"]] [[ as.character(exceedance_threshold_predictions[b]) ]]= W
            W = NULL
          }
      }

      if (!is.null(deceedance_threshold_predictions)) {
          
          Osamples[["predictions_deceedance"]]  = list()
          Osamples[["predictions_deceedance"]] [["deceedance"]] = list()

          for ( b in 1:length(deceedance_threshold_predictions)) {
            m = apply ( predictions, 1, FUN=function(x) length( which(x < deceedance_threshold_predictions[b] ) ) ) / O[["nposteriors"]]
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
            Osamples[["predictions_deceedance"]] [["deceedance"]] [[ as.character(deceedance_threshold_predictions[b]) ]]= W
            W = NULL
          }
      }
          
      predictions = NULL

      # save summary
      fn_samples = gsub( "_fit~", "_samples~", fn_fit, fixed=TRUE )
      read_write_fast( data=Osamples, file=fn_samples, compress=compress, compression_level=compression_level, qs_preset=qs_preset ) 
      Osamples = NULL 
      gc()
      


    } # END predictions posteriors loop
  } # END redo posterior samples


  end_post_samples  = Sys.time()


  # --------------------------------------------
  # --------------------------------------------
  # --------------------------------------------
   
 
  
  if (is.null(fn_res)) {
      message( "Return object is 'fit$results' (and not 'fit')")
      fit$results = O
      read_write_fast( data=fit, file=fn_fit, compress=compress, compression_level=compression_level, qs_preset=qs_preset )
  } else {
      message( "Saving results summary as: \n", fn_res )
      read_write_fast( data=O, file=fn_res, compress=compress, compression_level=compression_level, qs_preset=qs_preset )
  }
  

  run_end  = Sys.time()
 
  # --------------------------------------------
  # --------------------------------------------
  # --------------------------------------------

  if (be_verbose) {
    message("---------------------------------------")
    message("Run times:")
    message("Fit: ", format(difftime(run_fit_end, run_start)))
    message("Extract results from marginals: ", format(difftime(run_post_samples, run_fit_end)))
    message("Posterior simulations: ", format(difftime(end_post_samples, run_post_samples)))
    message("Total: ", format(difftime(run_end, run_start)))
    message("---------------------------------------")
  } 

  return(O)
  
}
