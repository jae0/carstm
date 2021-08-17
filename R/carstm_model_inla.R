
carstm_model_inla = function(p, M=NULL, E=NULL, sppoly=NULL, region.id=NULL,
  fn_fit=tempfile(pattern="fit_", fileext=".Rdata"), 
  fn_res=tempfile(pattern="res_", fileext=".Rdata"), 
  num.threads="1:1",
  compression_level=1,
  redo_fit = TRUE,
  update_results = FALSE,
  toget = c("summary", "fixed_effects", "random_other", "random_spatial", "random_spatiotemporal" , "predictions"), 
  exceedance_threshold=NULL, deceedance_threshold=NULL, nposteriors=NULL, 
  improve.hyperparam.estimates=NULL,  ... ) {

  ellp = list(...)

  inla.setOption(num.threads=num.threads)

  num.cores =  as.numeric( unlist(strsplit(num.threads, ":") )[1] )
  
  # local functions
  apply_generic = function(...)  mclapply(...,   mc.cores=num.cores ) # drop-in for lapply
  apply_simplify = function(...) simplify2array(mclapply(...,  mc.cores=num.cores ))  # drop in for sapply

  if (0) {
    # for debugging
    ellp = list()
    E=NULL 
    sppoly=NULL
    region.id=NULL
    fn_fit=tempfile(pattern="fit_", fileext=".Rdata")
    fn_res=tempfile(pattern="res_", fileext=".Rdata")
    compression_level=1
    redo_fit = TRUE
    update_results = FALSE
    toget = c("summary", "fixed_effects", "random_other", "random_spatial", "random_spatiotemporal" , "predictions")
    exceedance_threshold=NULL
    deceedance_threshold=NULL
    nposteriors=NULL
    improve.hyperparam.estimates=NULL

    if (0) {
        # usual variable names used in aegis .. char / num
        p$vnS = "space"  # "space"
        p$vnT = "time"
        p$vnU = "season" # "dyri"  # sub annual time

        # alt charracter descrptions of vars (if any, taken from main above if not)
        # p$vnS0 = "AUID"  # as character
        # p$vnT0 = "yr"  # as character
        # p$vnU0 = "dyear"

        # copies of main effects for inla in formulae
        p$vnST = "space_time"  # vnST = "space_time" (copy of vnS)
        p$vnTS = "time_space"  # vnTS = "time_space" (copy of vnT)
    }
  }

  # labels .. not all used but defining the here makes things simpler below
  vnS = ifelse( exists("vnS", p), p$vnS, "space" )
  vnS0 = ifelse( exists("vnS0", p), p$vnS0, "space0" ) # local storage/copy as a character
  vnSI = ifelse( exists("vnSI", p), p$vnSI, "space_iid" )
  vnSTI = ifelse( exists("vnSTI", p), p$vnSTI, "space_time_iid" )
  vnT = ifelse( exists("vnT", p), p$vnT, "time" )
  vnT0 = ifelse( exists("vnT0", p), p$vnT0, "time0" )
  vnST = ifelse( exists("vnST", p), p$vnST, "space_time" )
  vnTS = ifelse( exists("vnTS", p), p$vnTS, "time_space" )
  vnU = ifelse( exists("vnU", p), p$vnU, "season" )  
  vnU0 = ifelse( exists("vnU0", p), p$vnU0, "season0" )
  vnO = ifelse( exists("vnO", p), p$vnO, "data_offset" )  # input data is expected to be on user scale
  vnY = p$variabletomodel

  be_verbose = FALSE
  if ( exists("verbose", ellp ) ) if ( ellp[["verbose"]] )   be_verbose = TRUE

  # outputs
  O = NULL
  if ( update_results) try( load(fn_res) )
  if ( inherits(O, "try-error")) O = list()
  if ( is.null(O)) O = list()
  
  if (is.null(nposteriors))  nposteriors = ifelse( exists("nposteriors", p), p$nposteriors, 5000 )

  O[["formula"]] = p$carstm_model_formula
  O[["formula_parsed"]] = fm = carstm_parse_formula(p$carstm_model_formula )

  if ( fm$dependent_variable != vnY){
    message("Possible parsing issue. Dependent variable in params is specified as: ",  vnY )
    message("But formula is using: ", fm$dependent_variable )
    stop()
  } 

  if (!exists("carstm_model_family", p )) p$carstm_model_family = "gaussian"

  # gaussian by default
  lnk_function = inla.link.identity
  invlink_fixed = inla.link.identity
  invlink_random = inla.link.identity
  invlink_predictions = inla.link.identity ## preds are on user scale

  if ( p$carstm_model_family == "lognormal" ) {
    lnk_function = inla.link.log
    invlink_fixed = function(x) lnk_function( x,  inverse=TRUE )
    invlink_random = function(x) lnk_function( x,  inverse=TRUE )
    invlink_predictions = function(x) lnk_function( x,  inverse=TRUE )  
 
  } else if ( grepl( ".*poisson", p$carstm_model_family)) {
    lnk_function = inla.link.log
    invlink_fixed = function(x) lnk_function( x,  inverse=TRUE )
    invlink_random = function(x) lnk_function( x,  inverse=TRUE )
    invlink_predictions = function(x) lnk_function( x,  inverse=TRUE )

  } else if ( grepl( ".*binomial", p$carstm_model_family)) {
    lnk_function = inla.link.logit
    invlink_fixed = function(x) lnk_function( x,  inverse=TRUE )
    invlink_random = function(x) lnk_function( x,  inverse=TRUE )
    invlink_predictions = function(x) lnk_function( x,  inverse=TRUE )
  } 


  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )
  if (is.null(M)) if (exists("M", O)) M = O$M
  setDF(M)

  if (is.null(region.id)) {
    if (is.null(sppoly)) sppoly = areal_units( p=p )  # required by car fit
    if (is.null(sppoly)) if (exists("sppoly", O)) sppoly = O$sppoly
    if (!is.null(sppoly)) {
      if (is.null(region.id)) {
        if (exists("AUID", sppoly)) {
          region.id = as.character( sppoly[["AUID"]] )  # the master / correct sequence of the AU's and neighbourhood matrix index values
        }
      }
      if (is.null(region.id)) {
        region.id = try( as.character( slot( slot(sppoly, "nb"), "region.id" ) ) )
        if (inherits(region.id, "try-error")) region.id = NULL
      }
      if (is.null(region.id)) {
        region.id = try( as.character( slot( slot(sppoly, "W.nb"), "region.id" ) ) )
        if (inherits(region.id, "try-error")) region.id = NULL
      }
      if (is.null(region.id)) {
        region.id = try( as.character( slot( sppoly,  "region.id" ) ) )
        if (inherits(region.id, "try-error")) region.id = NULL
      }
    }
  }


  if (is.null(region.id)) stop("Not found: region.id or sppoly$AUID is a required variable ")

  nAUID = length(region.id)
  O[[vnS]] = as.character( region.id )  # this sequence is a master key as index matches nb matrix values

  fit  = NULL

  if (redo_fit) {


    # fiddling of AU and TU inputs: for bym2, etc, they need to be numeric, matching numerics of polygon id ("region.id")
    # convert space and time to numeric codes for INLA
    if (grepl("space", p$aegis_dimensionality)) {
      M[,vnS0] = as.character( M[,vnS] ) # local copy
      M[,vnS] = match( M[,vnS0], O[[vnS]] )  # overwrite with numeric values that must match index of neighbourhood matrix
    }

    if ( grepl("year", p$aegis_dimensionality) | grepl("season", p$aegis_dimensionality) ) {
      if (any( grepl( vnT, p$carstm_model_formula )))  {
        O[[vnT]] = as.character( p$yrs )
        M[,vnT0] = as.character( M[,vnT] ) # a copy for internal matching 
        M[,vnT] = match( M[,vnT0], O[[vnT]] ) # convert to numeric (ie. a numeric factor)
      }
      # internal vars, for inla
      if (any( grepl( vnST, p$carstm_model_formula )))  M[,vnST] = M[,vnS]
      if (any( grepl( vnST, p$carstm_model_formula )))  M[,vnTS] = M[,vnT]
      # sub-annual time
      if (any( grepl( vnU, p$carstm_model_formula )))  {
        O[[vnU]] = as.character( p$dyears + diff(p$dyears)[1]/2)
        M[,vnU0] = as.character( M[,vnU] )  # a copy for internal matching 
        M[,vnU] = match( M[,vnU0], O[[vnU]] )  # convert to numeric (ie. a numeric factor)
      }
    }

    # on user scale
    ii = which(is.finite(M[ , vnY ]))
    
    mq = quantile( M[ ii, vnY ], probs=c(0.025, 0.5, 0.975) )

    O[["data_range"]] = c( 
      mean=mean(M[ ii, vnY ]), 
      sd=sd(M[ ii, vnY ]), 
      min=min(M[ ii, vnY ]), 
      max=max(M[ ii, vnY ]),  
      lb=mq[1], 
      median=mq[2], 
      ub=mq[3]  
    )  # on data /user scale not internal link

    # prefilter/transformation (e.g. translation to make all positive)
    if (exists("data_transformation", p)) M[, vnY]  = p$data_transformation$forward( M[, vnY] ) 

    # get hyper param scalings

    # temp Y var on link scale:
    yl = lnk_function( M[, vnY ] )   # necessary in case of log(0)

    # offsets need to be close to 1 in user scale ( that is log(1)==0 in internal scale ) in experimental mode .. rescale  
    if ( !is.null(fm$offset_variable) )  {
      if ( fm$offset_variable != vnO ){
        message("Possible parsing issue. Offset variable in params is specified as: ",  vnO )
        message("But formula is using: ", fm$offset_variable )
        stop()
      } 

      # link function applied to offsets here .. do not need to send log() 
      if (grepl("log", vnO)) message("Probably do not want to transform the offset .. it is done internally in carstm, unlike glm, inla, etc")

      obs = 1:nrow(M)
      if (exists("tag", M)) {
        obso = which(M$tag=="observations")
        if (length(obso) > 3) obs = obso
        obso = NULL
      }   
      
      ol = lnk_function( M[, vnO ])
      O$offset_scale = median( ol[obs] , na.rm=TRUE )  # required to stabilize optimization
      obs = NULL
      ol = ol - O$offset_scale  # apply to all and overwrite, certing upon 0 (in user space 1)
      M[ , vnO ] = ol 
      O$offset_scale_revert = function(x) { x - O$offset_scale }  # used for Y value and so opposite sign of "-" but applied to numerator is returns it to "-"
      yl = yl - ol
    } 

    ll = which(is.finite(yl))
    H = carstm_hyperparameters( sd(yl[ll] ), alpha=0.5, median(yl[ll], na.rm=TRUE) )  # sd slightly biased due to 0's being dropped

    m = yl = ii = ll = fy = ol = NULL
    gc()

    

    ellp[["formula"]] = p$carstm_model_formula
    ellp[["family"]] = p$carstm_model_family
    
    ellp[["data"]] = M
    M = NULL
    
    if (!is.null(E)) ellp[["E"]] = E
    E = NULL
    
    gc()

    if ( !exists("inla.mode", ellp ) ) ellp[["inla.mode"]] = "experimental"
    if ( !exists("control.inla", ellp ) ) ellp[["control.inla"]] = list( strategy='adaptive' )
    if ( !exists("control.predictor", ellp ) ) ellp[["control.predictor"]] = list(compute=TRUE, link=1  ) #everything on link scale
    if ( !exists("control.compute", ellp ) ) ellp[["control.compute"]] = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE )


    if ( ellp[["inla.mode"]] != "experimental") {
      # location of this option has moved ... might move again
      if ( !exists("control.results", ellp ) ) ellp[["control.results"]] = list(return.marginals.random=TRUE, return.marginals.predictor=TRUE )
      ellp[["control.compute"]][["return.marginals.predictor"]] = NULL
    }


    # control.fixed= list(mean.intercept=0, prec.intercept=0.001, mean=0, prec=0.001),
    # control.inla = list( strategy='adaptive', int.strategy='eb' )
 
    fit = try( do.call( inla, ellp ) )      

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
    
    M = ellp[["data"]]  # rename for followup
    ellp = NULL; gc()
      
    if (be_verbose)  message( "Saving carstm fit (this can be slow): ", fn_fit )

    save( fit, file=fn_fit, compression_level=compression_level )

  }

  if (is.null(fit)) load( fn_fit )
  
  if (is.null(fit)) {
    message( "fit file not found: ", fn_fit )
    stop()
  }

  # do the computations here as fit can be massive ... best not to copy, etc ..
  if (be_verbose)  message( "\nComputing summaries and computing from posterior simulations (can be longer than model fitting depending upon no of posterior sims: 'nposteriors' ) ..." )

  
  list_simplify = function(x) as.data.frame( t( as.data.frame( x )))
  exceedance_prob = function(x, threshold)  {1 - inla.pmarginal(q = threshold, x)}
  deceedance_prob = function(x, threshold)  { inla.pmarginal(q = threshold, x)}


  tokeep =  c("mean", "sd", "quant0.025", "quant0.5", "quant0.975")

  # exceedance_threshold=1
  # deceedance_threshold=1

  if (exists("deceedance_threshold", p)) deceedance_threshold = p[["deceedance_threshold"]]
  if (exists("exceedance_threshold", p)) exceedance_threshold = p[["exceedance_threshold"]]


  if ( "summary" %in% toget) {

    if (be_verbose)  message("Extracting from marginals: parameter summaries"  )
    if (!exists("summary", O)) O[["summary"]] = list()

    O[["summary"]][["direct"]] = summary(fit)

    # parameters
    # back-transform from marginals

    if (exists( "marginals.fixed", fit)) {
      V = fit$marginals.fixed

      if (exists("offset_scale_revert", O)) {
        fi = which( grepl("Intercept", names(V) ))
        if (length(fi)>0) V[[fi]] = inla.tmarginal( O$offset_scale_revert, V[[fi]]) 
      }

      summary_inv_fixed = function(x) inla.zmarginal( inla.tmarginal( invlink_fixed, x) , silent=TRUE  )
      W = NULL
      W = cbind ( t (apply_simplify( V, FUN=summary_inv_fixed ) ) )  # 
      O[["summary"]][["fixed_effects"]] = W [, tokeep, drop =FALSE]
      W = NULL
    }

    # hyperpar (variance components)
    j = grep( "^Precision.*", rownames(fit$summary.hyperpar), value=TRUE )
    if (length(j) > 0) {

      summary_inv_prec = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(pmax(y,1e-12)), x) , silent=TRUE  )
      # summary_inv_prec_1024 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=1024L) , silent=TRUE  )
      # summary_inv_prec_512 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=512L) , silent=TRUE  )

      precs = try( list_simplify( apply_simplify( fit$marginals.hyperpar[j], FUN=summary_inv_prec ) ), silent=TRUE )  # prone to integration errors ..
      if (inherits(precs, "try-error")) precs = try( list_simplify( apply_simplify( fit$marginals.hyperpar[j], FUN=summary_inv_prec_1024 ) ), silent=TRUE )
      if (inherits(precs, "try-error")) {
        if (be_verbose)  message( "Model may be over parameterized. NAN and Inf values encountered. Try alt parameterizations or smaller number of n or masking negative values")
        precs = fit$summary.hyperpar[j,1:5]
        precs[,c(1,3:5)] = 1/sqrt( precs[,c(1,3:5)] )
        rownames(precs) = gsub("Precision for", "SD", rownames(precs) )
        colnames(precs) = tokeep
        O[["summary"]][["random_effects"]] = precs

      } else {
        # precs[,"mode"] =  1/sqrt( fit$summary.hyperpar[j,"mode"]  )
        toadd = setdiff( colnames(O[["summary"]]), colnames(precs) )
        precs[,toadd] = NA
        rownames(precs) = gsub("Precision for", "SD", rownames(precs) )
        O[["summary"]][["random_effects"]] = precs[, tokeep, drop =FALSE]
      }
    }
    precs = NULL

    j = grep( ".*Rho.*", rownames(fit$summary.hyperpar), value=TRUE )
    if (length(j) > 0) {
      rhos = try( list_simplify( apply_simplify( fit$marginals.hyperpar[j], FUN=function(x) inla.zmarginal( x, silent=TRUE  ) ) ), silent=TRUE )
      if (any( inherits(rhos, "try-error"))) {
        if (be_verbose)  message( "Model may be over parameterized. NAN and Inf values encountered. Try alt parameterizations or smaller number of n or masking negative values")
      } else {
        #  rhos[,"mode"] = apply_simplify( fit$marginals.hyperpar[j], FUN=function(x) inla.mmarginal( x ))
        O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], rhos[, tokeep, drop =FALSE] )
        rhos = NULL
      }
    }

    # update phi's
    j = grep( "^Phi.*", rownames(fit$summary.hyperpar), value=TRUE )
    if (length(j) > 0) {
      phis = fit$summary.hyperpar[j, ,drop=FALSE]
      colnames(phis) = c( tokeep, "mode")
      O[["summary"]][["random_effects"]] = rbind( O[["summary"]][["random_effects"]], phis[, tokeep, drop =FALSE] )
      phis = NULL
    }

    j = NULL
    gc()

    print(  O[["summary"]] )

  }


  if (any( grepl("random", toget) )) {
    if (!exists("random", O)) O[["random"]] = list()
    if (exists("marginals.random", fit)) {

      if ("random_other" %in% toget) {
        summary_inv_random = function(x) inla.zmarginal( inla.tmarginal( invlink_random, x) , silent=TRUE  )
        raneff = setdiff( names( fit$marginals.random ), c(vnS, vnST, vnSI, vnSTI ) )
        for (re in raneff) {
          if (be_verbose)  message("Extracting from marginals: random covariate ", re  )
          g = fit$marginals.random[[re]]
          O[["random"]] [[re]] = list_simplify ( apply_simplify( g, summary_inv_random ) )  [, tokeep, drop =FALSE]
          O[["random"]] [[re]]$ID = fit$summary.random[[re]]$ID
        }
        g = raneff = NULL
      }
      gc()

      if ("random_spatial" %in% toget) {
        # space only
        if (be_verbose)  message("Extracting from marginals: random spatial errors"  )
        summary_inv_random = function(x) inla.zmarginal( inla.tmarginal( invlink_random, x) , silent=TRUE  )

        bym = iid = NULL
        matchto = list( space=O[[vnS]] )

        if ( exists(vnSI, fit$marginals.random) ) {
          O[["random"]] [[vnSI]] = list()  # space as a main effect
          model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnSI) ]  # should be iid

          # single spatial effect (eg besag) .. indexing not needed but here in case more complex models ..
          Z = expand.grid( space=O[[vnSI]], type=model_name, stringsAsFactors =FALSE )
          iid =  which(Z$type==model_name)
          matchfrom0 = list( space=Z[["space"]][iid] )

          m = list_simplify ( apply_simplify( fit$marginals.random[[vnSI]], summary_inv_random ) )

          #  spatial effect ... besag, etc main effects
          W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) )
          names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
          matchfrom = list( space=Z[["space"]][iid] )
          for (k in 1:length(names(m))) {
            W[,k] = reformat_to_array( input = unlist(m[iid,k]), matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vnSI]] [[model_name]] = W [, tokeep, drop =FALSE]
          W = NULL
        }

        if ( exists(vnS, fit$marginals.random) ) {
          O[["random"]] [[vnS]] = list()  # space as a main effect
          model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnS) ]
          if ( model_name %in% c("bym", "bym2") ) {
            # bym2 effect: bym and iid
            Z = expand.grid( space=O[[vnS]], type = c("iid", model_name), stringsAsFactors =FALSE )
            bym = which(Z$type==model_name)
            iid = which(Z$type=="iid")
            matchfrom0 = list( space=Z[["space"]][iid] )
          } else {
            # single spatial effect (eg besag, etc)
            Z = expand.grid( space=O[[vnS]], type=model_name, stringsAsFactors =FALSE )
            bym = which(Z$type==model_name)
            matchfrom0 = list( space=Z[["space"]][bym] )
          }

          m = list_simplify ( apply_simplify( fit$marginals.random[[vnS]], summary_inv_random ) )

          #  spatial effect ... besag, etc main effects
          W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) )
          names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
          matchfrom = list( space=Z[["space"]][bym] )
          for (k in 1:length(names(m))) {
            W[,k] = reformat_to_array( input = unlist(m[bym,k]), matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vnS]] [[model_name]] = W [, tokeep, drop =FALSE]

          if  ( model_name %in% c("bym", "bym2") ) {
            #  iid main effects
            W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            matchfrom = list( space=Z[["space"]][iid] )
            for (k in 1:length(names(m))) {
              W[,k] = reformat_to_array( input = unlist(m[iid,k]), matchfrom=matchfrom, matchto=matchto )
            }
            O[["random"]] [[vnS]] [["iid"]] = W [, tokeep, drop =FALSE]
            W = NULL
          }
          m = NULL
        }


        # POSTERIOR SIMS
        if ( exists(vnSI, fit$marginals.random ) &  exists(vnS, fit$marginals.random ) )  {
          # sep besag + iid  
          selection=list()
          selection[vnSI] = 0  # 0 means everything matching space
          selection[vnS] = 0  # 0 means everything matching space
          aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
          aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
          iid = which(aa_rn==vnSI)
          bym = which(aa_rn==vnS)
          g = apply_simplify( aa, function(x) {invlink_random(x$latent[iid] + x$latent[bym] ) } )
        } else if ( exists(vnSI, fit$marginals.random ) & ! exists(vnS, fit$marginals.random ) )  {
          # iid  only
          selection=list()
          selection[vnSI] = 0  # 0 means everything matching space
          aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
          aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
          iid = which(aa_rn==vnSI)
          g = apply_simplify( aa, function(x) {invlink_random(x$latent[iid] ) } )
        } else if ( ! exists(vnSI, fit$marginals.random ) & exists(vnS, fit$marginals.random ) ) {
          # besag only or bym/bym2
          selection=list()
          selection[vnS] = 0  # 0 means everything matching space
          aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
          if ( model_name %in% c("bym", "bym2") ) {
            g = apply_simplify( aa, function(x) {invlink_random(x$latent[iid] + x$latent[bym] ) } )
          } else {
            aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
            bym = which(aa_rn==vnS)
            g = apply_simplify( aa, function(x) {invlink_random(x$latent[bym] ) } )
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

        if ( !is.null(exceedance_threshold)  | !is.null(deceedance_threshold) ) {
          if (!is.null(exceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( input = m, matchfrom=matchfrom0, matchto = matchto )
            names(dimnames(W))[1] = vnS
            dimnames( W )[[vnS]] = O[[vnS]]
            O[["random"]] [[vnS]] [["exceedance"]] = W
            W = m = NULL
          }

          if (!is.null(deceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( input = m, matchfrom=matchfrom0, matchto=matchto  )
            names(dimnames(W))[1] = vnS
            dimnames( W )[[vnS]] = O[[vnS]]
            O[["random"]] [[vnS]] [["deceedance"]] = W
            W = m = NULL
          }
        }
      }

      Z = g = NULL
      gc()


      if ("random_spatiotemporal"  %in% toget ) {
        # space-year
        if (be_verbose)  message("Extracting from marginals: random spatiotemporal errors"  )
        summary_inv_random = function(x) inla.zmarginal( inla.tmarginal( invlink_random, x) , silent=TRUE  )

        bym = iid = NULL
        matchto = list( space=O[[vnS]], time=O[[vnT]]  )
      
        if (exists(vnSTI, fit$marginals.random )) {
          O[["random"]] [[vnSTI]] = list()
          model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnSTI) ]  # should be iid

          Z = expand.grid( space=O[[vnS]], type=model_name, time=O[[vnT]], stringsAsFactors =FALSE )
          iid =  which(Z$type==model_name)
          matchfrom0 = list( space=Z[["space"]][iid], time=Z[["time"]][iid]  )

          m = list_simplify ( apply_simplify( fit$marginals.random[[vnSTI]], summary_inv_random ) )

          #  spatiotemporal interaction effects  iid
          W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
          names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
          names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
          matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  )
          for (k in 1:length(names(m))) {
            W[,,k] = reformat_to_array(  input = unlist(m[iid,k]), matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vnSTI]] [[model_name]] = W [,, tokeep, drop =FALSE]
          W = NULL
          iid = NULL
        }

        if (exists(vnST, fit$marginals.random )) {
          O[["random"]] [[vnST]] = list()
          model_name = fm$random_effects$model[ which(fm$random_effects$vn == vnST) ]

          if ( model_name %in% c("bym", "bym2") ) {
            # bym2 effect: bym and iid with annual results
            Z = expand.grid( space=O[[vnS]], type = c("iid", model_name), time=O[[vnT]], stringsAsFactors =FALSE )
            bym = which(Z$type==model_name)
            iid = which(Z$type=="iid")
            matchfrom0 = list( space=Z[["space"]][iid], time=Z[["time"]][iid]  )
          } else {
            # besag effect: with annual results
            Z = expand.grid( space=O[[vnS]], type =model_name, time=O[[vnT]], stringsAsFactors =FALSE )
            bym =  which(Z$type==model_name)
            matchfrom0 = list( space=Z[["space"]][bym], time=Z[["time"]][bym]  )
          }

          m = list_simplify ( apply_simplify( fit$marginals.random[[vnST]], summary_inv_random ) )

          #  spatiotemporal interaction effects  bym
          W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
          names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
          names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
          matchfrom = list( space=Z[["space"]][bym],  time=Z[["time"]][bym]  )
          for (k in 1:length(names(m))) {
            W[,,k] = reformat_to_array(  input = unlist(m[bym,k]), matchfrom = matchfrom, matchto = matchto  )
          }
          O[["random"]] [[vnST]] [[model_name]] = W [,, tokeep, drop =FALSE]
          W = NULL

          if ( model_name %in% c("bym", "bym2") ) {
            #  spatiotemporal interaction effects  iid
            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
            matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  )
            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = unlist(m[iid,k]), matchfrom = matchfrom, matchto = matchto )
            }
            O[["random"]] [[vnST]] [["iid"]] = W [,, tokeep, drop =FALSE]
            W = NULL
          }
          m = NULL

          if ( exists(vnSTI, fit$marginals.random ) &  exists(vnST, fit$marginals.random ) ) {
            # sep besag +iid  
              selection=list()
              selection[vnSTI] = 0  # 0 means everything matching space
              selection[vnST] = 0  # 0 means everything matching space
              aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
              aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
              iid = which(aa_rn==vnSTI)
              bym = which(aa_rn==vnST)
              g = apply_simplify( aa, function(x) {invlink_random(x$latent[iid] + x$latent[bym] ) } )
              matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  )
          } else if (exists(vnSTI, fit$marginals.random ) & ! exists(vnST, fit$marginals.random )) {
            # iid  
              selection=list()
              selection[vnSTI] = 0  # 0 means everything matching space
              aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
              aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
              iid = which(aa_rn==vnSTI)
              g = apply_simplify( aa, function(x) {invlink_random(x$latent[iid] ) } )
              matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  )
          } else if ( !exists(vnSTI, fit$marginals.random ) & exists(vnST, fit$marginals.random ) ) {
            # besag  only or bym/bym2
              selection=list()
              selection[vnST] = 0  # 0 means everything matching space
              aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
              if ( model_name %in% c("bym", "bym2") ) {
                g = apply_simplify( aa, function(x) {invlink_random(x$latent[iid] + x$latent[bym] ) } )
                matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  )
              } else {
                aa_rn = gsub( "[:].*$", "", rownames(aa[[1]]$latent) )
                bym = which(aa_rn==vnST)
                g = apply_simplify( aa, function(x) {invlink_random(x$latent[bym] ) } )
                matchfrom = list( space=Z[["space"]][bym],  time=Z[["time"]][bym]  )
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
            W[,,k] = reformat_to_array(  input = m[,k], matchfrom=matchfrom, matchto=matchto )
          }
          O[["random"]] [[vnST]] [["combined"]] = W [,, tokeep, drop =FALSE]
          W = NULL
        }


        if ( !is.null(exceedance_threshold)  | !is.null(deceedance_threshold) ) {
      
          if (!is.null(exceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( input=m, matchfrom=matchfrom0,  matchto=matchto )
            m = NULL
            names(dimnames(W))[1] = vnS
            names(dimnames(W))[2] = vnT
            dimnames( W )[[vnS]] = O[[vnS]]
            dimnames( W )[[vnT]] = O[[vnT]]
            O[["random"]] [[vnST]] [["exceedance"]] = W
            W = NULL
          }

          if (!is.null(deceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( input = m, matchfrom=matchfrom0,  matchto=matchto )
            m = NULL
            names(dimnames(W))[1] = vnS
            names(dimnames(W))[2] = vnT
            dimnames( W )[[vnS]] = O[[vnS]]
            dimnames( W )[[vnT]] = O[[vnT]]

            O[["random"]] [[vnST]] [["deceedance"]] = W
            W = NULL
          }  
          # ned space-year
        }
        
        Z = g = NULL
        gc()
      }
    }  # end random effects
  }


  if ("predictions"  %in% toget ) {

    if (!exists("predictions", O)) O[["predictions"]] = list()

    truncate_upperbound = function( b, upper_limit, eps=1e-12 ) {
      # not used
      k = which( b[,1] > upper_limit )
      if (length(k) > 0) b[k,2] = 0
      return( b )
    }

    if (exists("data_transformation", p))  {
      backtransform = function( b ) {
        b[,1] =  p$data_transformation$backward( b[,1]   )
        return( b )
      }
    } 

    summary_inv_predictions = function(x) inla.zmarginal( x, silent=TRUE  )
    
    if (!exists("tag", M)) M$tag="predictions" # force predictions for all data

    # adjusted by offset
    if (exists("marginals.fitted.values", fit)) {

      if (be_verbose)  message("Extracting from marginals: predictions"  )

      if (  p$aegis_dimensionality == "space" ) {
        ipred = which( M$tag=="predictions"  &  M[,vnS0] %in% O[[vnS]] )  # filter by S and T in case additional data in other areas and times are used in the input data
        g = fit$marginals.fitted.values[ipred]   
        if ( exists("offset_scale_revert", O) ) g = apply_generic( g, function(u) {inla.tmarginal( O$offset_scale_revert, u) } )    

        g = apply_generic( g, function(u) {inla.tmarginal( invlink_predictions, u) } )    
        if (exists("data_transformation", p)) g = apply_generic( g, backtransform )

        m = list_simplify ( apply_simplify( g, summary_inv_predictions ) )
        W = array( NA, dim=c( length(O[[vnS]]),  length(names(m)) ),  dimnames=list( space=O[[vnS]], stat=names(m) ) )
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        
        matchfrom = list( space=M[ ipred, vnS0] ) 
        matchto = list( space=O[[vnS]] )

        for (k in 1:length(names(m))) {
          W[,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=matchfrom, matchto=matchto )
        }
        O[["predictions"]] = W[, tokeep, drop =FALSE]
        W = m = NULL
      }

      if ( p$aegis_dimensionality == "space-year" ) {
        ipred = which( M$tag=="predictions" & M[,vnS0] %in% O[[vnS]] & M[,vnT0] %in% O[[vnT]] )
        g = fit$marginals.fitted.values[ipred]   
        if ( exists("offset_scale_revert", O) ) g = apply_generic( g, function(u) {inla.tmarginal( O$offset_scale_revert, u) } )    

        g = apply_generic( g, function(u) {inla.tmarginal( invlink_predictions, u) } )    

        if (exists("data_transformation", p)) g = apply_generic( g, backtransform )

        m = list_simplify ( apply_simplify( g, summary_inv_predictions ) )
        W = array( NA, dim=c( length(O[[vnS]]), length(O[[vnT]]), length(names(m)) ),  dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) )
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
      
        matchfrom = list( space=M[ipred,vnS0], time=M[ipred,vnT0] )
        matchto = list( space=O[[vnS]], time=O[[vnT]] )

        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=matchfrom, matchto=matchto)
        }
        O[["predictions"]] = W[,, tokeep, drop =FALSE]
        W = m = NULL
      }


      if ( p$aegis_dimensionality == "space-year-season" ) {
        ipred = which( M$tag=="predictions" & M[,vnS0] %in% O[[vnS]]  &  M[,vnT0] %in% O[[vnT]] &  M[,vnU0] %in% O[[vnU]])  # ignoring U == predict at all seassonal components ..
        g = fit$marginals.fitted.values[ipred]   
        if ( exists("offset_scale_revert", O) ) g = apply_generic( g, function(u) {inla.tmarginal( O$offset_scale_revert, u) } )    
 
 
        g = apply_generic( g, function(u) {inla.tmarginal( invlink_predictions, u) } )    

        if (exists("data_transformation", p)) g = apply_generic( g, backtransform )

        m = list_simplify ( apply_simplify( g, summary_inv_predictions ) )
        W = array( NA, dim=c( length(O[[vnS]]), length(O[[vnT]]), length(O[[vnU]]), length(names(m)) ),  dimnames=list( space=O[[vnS]], time=O[[vnT]], season=O[[vnU]], stat=names(m) ) )
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
        names(dimnames(W))[3] = vnU  # need to do this in a separate step ..

        matchfrom = list( space=M[ipred, vnS0], time=M[ipred, vnT0], season=M[ipred, vnU0] )
        matchto = list( space=O[[vnS]], time=O[[vnT]], season=O[[vnU]] )

        for (k in 1:length(names(m))) {
          W[,,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=matchfrom, matchto=matchto )
        }
        O[["predictions"]] = W[,,, tokeep, drop =FALSE]
        W = m = NULL
      }

      if (!is.null(exceedance_threshold)) {
        m = list_simplify ( apply_simplify( g, FUN=exceedance_prob, threshold=exceedance_threshold ) )
        W = reformat_to_array(  input = unlist(m ), matchfrom = matchfrom, matchto = matchto )
        names(dimnames(W))[1] = vnS
        dimnames( W )[[vnS]] = O[[vnS]]
        if (!is.null(vnT)) {
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnT]] = O[[vnT]]
        }
        if (!is.null(vnU)) {
          names(dimnames(W))[3] = vnU
          dimnames( W )[[vnU]] = O[[vnU]]
        }
        
        O[["exceedance_probability"]] = W
        W = m = NULL
      }

      if (!is.null(deceedance_threshold)) {
        m = list_simplify ( apply_simplify( g, FUN=deceedance_prob, threshold=deceedance_threshold ) )
        W = reformat_to_array(  input = unlist(m ), matchfrom = matchfrom, matchto = matchto )
        names(dimnames(W))[1] = vnS
        dimnames( W )[[vnS]] = O[[vnS]]
        if (!is.null(vnT)) {
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnT]] = O[[vnT]]
        }
        if (!is.null(vnU)) {
          names(dimnames(W))[3] = vnU
          dimnames( W )[[vnU]] = O[[vnU]]
        }
        
        O[["deceedance_probability"]] = W
        W = m = NULL
      }
    }
  }
  
  fit = NULL
  gc()

  # copy data in case needed for plotting ..
  if (!is.null(M)) O[["M"]] = M
  if (!is.null(sppoly)) O[["sppoly"]] = sppoly

  save( O, file=fn_res, compression_level=compression_level )

  if (be_verbose)  message( "Carstm summary saved as: ", fn_res )

  return(O)

}