carstm_extract_inla = function(p, M, fn_fit, toget="summary", file_compress_method=FALSE, exceedance_threshold=NULL, deceedance_threshold=NULL, nposteriors=NULL, improve.hyperparam.estimates=NULL, quantile_limit=NULL ) {

  # outputs
  O = list()
  O[["summary"]] = list()
  O[["random"]] = list()

  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )


  if (0) {

      # usual variable names used in aegis .. INLA requires these to be numeric
      p$vnS = "auid_main"  # "space"
      p$vnT = "yr"
      p$vnU = "dyri"  # sub annual time 
      
      # alt character descrptions of vars
      p$vnSn = "AUID"  # as character 
      p$vnTn = "year"  # as character 
      p$vnUn = "dyear"
      
      p$vnST = "auid"  # vnST = "space_time" (copy of vnS)
      p$vnTS = "year_factor"  # vnTS = "time_space" (copy of vnT)
      
      p$vnSf = "AUID"  # "space" as a factor
      p$vnTf = "year_factor"  # "time" as a factor
      p$vnUf = "season_factor"  # "time" as a factor

      # AUID is character; auid is factor -> numeric 

  }

  vnY = p$variabletomodel

  vnO = ifelse( exists("vnO", p), p$vnS, "data_offset" )  # input data is expected to be on user scale

  # these are numeric codes representing each 
  vnS = ifelse( exists("vnS", p), p$vnS, "space" )
  vnT = ifelse( exists("vnT", p), p$vnT, "time" )
  vnU = ifelse( exists("vnU", p), p$vnU, "season" )  # sub-annual time

  # as.factor
  vnSf = ifelse( exists("vnSf", p), p$vnSf, "space_factor" )
  vnTf = ifelse( exists("vnTf", p), p$vnTf, "time_factor" )
  vnUf = ifelse( exists("vnUf", p), p$vnUf, "season_factor" )

  # alt character descriptions (names)
  vnSn = ifelse( exists("vnSn", p), p$vnSn, "space_name" ) 
  vnTn = ifelse( exists("vnTn", p), p$vnTn, "time_name" ) 
  vnUn = ifelse( exists("vnUn", p), p$vnUn, "season_name" ) 


  vnST = ifelse( exists("vnST", p), p$vnST, "space_time" )
  vnTS = ifelse( exists("vnTS", p), p$vnTS, "time_space" )

  vnST = ifelse( exists("vnST", p), p$vnST, "space_time" )
  vnTS = ifelse( exists("vnTS", p), p$vnTS, "time_space" )

  # copies

  if (!exists(vnST, M))  M[,vnST] = M[,vnS] 
  if (!exists(vnTS, M))  M[,vnTS] = M[,vnT] 

  # these should already be provided
  if (!exists(vnSf, M))  M[,vnSf] = as.factor( M[,vnS] )
  if (!exists(vnTf, M))  M[,vnTf] = as.factor( M[,vnT] )

  if (!exists(vnSn, M))  M[,vnSn] = as.character( M[,vnS] )
  if (!exists(vnTn, M))  M[,vnTn] = as.character( M[,vnT] )
  if (!exists(vnUn, M))  M[,vnUn] = as.character( M[,vnU] )

  # covert chars to numeric for INLA
   if (exists(vnSf, M)) M[,vnS] = as.numeric(M[,vnSf])
   if (exists(vnTf, M)) M[,vnT] = as.numeric(M[,vnTf])
   if (exists(vnUf, M)) M[,vnU] = as.numeric(M[,vnUf])


  if (is.null(nposteriors)) {   
    nposteriors = ifelse( exists("nposteriors", p), p$nposteriors, 1000 )
  }  

  
  # TODO: add extraction for besag, bym, etc..  


  sppoly = areal_units( p=p )  # required by car fit
 
  if (exists("data_transformation", p)) M[, vnY]  = p$data_transformation$forward( M[, vnY] ) # make all positive


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
  j = which( is.finite(M[,vnY]) )
  if (any( grepl( vnO, p$carstm_model_formula ))) {
    m = lnk_function( M[ j, vnY ] / M[ j, vnO ] )
  } else {
    m = lnk_function( M[ j, vnY ] )
  }


  ii = which(is.finite(m))
  mrange = range( m[ii]  )  # on data scale not internal
  if (!is.null(quantile_limit)) mlimit = quantile( m[ii], probs=quantile_limit )
  H = carstm_hyperparameters( sd(m[ii] ), alpha=0.5, median(m[ii], na.rm=TRUE) )
  m = NULL
  gc()

  p = parameters_add_without_overwriting( p, options.control.family = inla.set.control.family.default() )

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
  if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.6, diff.logdens=12 )  # get improved estimates for the hyperparameters

  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")

  message( "Saving carstm fit (this can be slow): ", fn_fit )

  save( fit, file=fn_fit, compress=file_compress_method )


  # do the computations here as fit can be massive ... best not to copy, etc ..
  message( "Computing summaries ..." )



  O[["dimensionality"]] = p$aegis_dimensionality

  if ( O[["dimensionality"]] == "space") {
    O[[vnSn]] = as.character( slot( slot(sppoly, "nb"), "region.id" ) )
    O[[vnSf]] = as.factor( O[[ vnSn ]] )
    ipred = which( M$tag=="predictions" & M[,vnSn] %in% O[[vnSn]] )  # filter by S and T in case additional data in other areas and times are used in the input data
    mfrom = list( space=M[ ipred, vnSn] )
    mto   = list( space=O[[vnSn]] )
  }

  if ( O[["dimensionality"]] == "space-year") {
    O[[vnSn]] = as.character( slot( slot(sppoly, "nb"), "region.id" ) )
    O[[vnTn]] = as.character( p$yrs )
    O[[vnSf]] = as.factor( O[[ vnSn ]] )
    O[[vnTf]] = as.factor( O[[ vnTn ]] )

    M[,vnTn] = as.character(M[,vnTn])

    # filter by vnSn and years in case additional data in other areas and times are used in the input data
    ipred = which( M$tag=="predictions" & M[,vnSn] %in% O[[vnSn]] & M[,vnTn] %in% O[[vnTn]] )  
    mfrom = list( space=M[ipred,vnSn], time=M[ipred,vnTn] )
    mto   = list( space=O[[vnSn]], time=O[[vnTn]] )
  }

  if ( O[["dimensionality"]] == "space-year-season") {
    O[[vnSn]] = as.character( slot( slot(sppoly, "nb"), "region.id" ) )
    O[[vnTn]] = as.character( p$yrs )
    O[[vnUn]] = as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[[vnUn]] ) )
    O[[vnSf]] = as.factor( O[[ vnSn ]] )
    O[[vnTf]] = as.factor( O[[ vnTn ]] )
    O[[vnUf]] = as.factor( O[[ vnUn ]] )

    M[,vnTn] = as.character( M[,vnTn] )
    M[,vnUn] = as.character( discretize_data( M[,vnUn], p$discretization[[vnUn]] ) )

    ipred = which( M$tag=="predictions" & M[,vnSn] %in% O[[vnSn]]  &  M[,vnTn] %in% O[[vnTn]] )  # ignoring U == predict at all seassonal components ..
    mfrom = list( space=M[ipred,vnSn], time=M[ipred,vnTn], season=M[ipred,vnUn] )
    mto   = list( space=O[[vnSn]], time=O[[vnTn]], season=O[[vnUn]] )
  }




  invlink = function(x) lnk_function( x,  inverse=TRUE ) 

  summary_inv = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )

  summary_inv_prec = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(pmax(y,1e-12)), x) , silent=TRUE  )
  # summary_inv_prec_1024 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=1024L) , silent=TRUE  )
  # summary_inv_prec_512 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=512L) , silent=TRUE  )

  list_simplify = function(x) as.data.frame( t( as.data.frame( x )))

  exceedance_prob = function(x, threshold)  {1 - inla.pmarginal(q = threshold, x)}

  deceedance_prob = function(x, threshold)  { inla.pmarginal(q = threshold, x)}

  truncate_upperbound = function( b, mlimit, eps=1e-12 ) {
    k = which( b[,1] > mlimit )  
    if (length(k) > 0) b[k,2] = 0
    return( b )
  }

  tokeep =  c("mean", "sd", "quant0.025", "quant0.5", "quant0.975")

  # exceedance_threshold=1 
  # deceedance_threshold=1

  if (exists("deceedance_threshold", p)) deceedance_threshold=p[["deceedance_threshold"]]
  if (exists("exceedance_threshold", p)) exceedance_threshold=p[["exceedance_threshold"]]



  if ( "summary" %in% toget) {
    
      O[["summary"]][["direct"]] = summary(fit)
      
      # parameters
      # back-transform from marginals
      W = cbind ( t (sapply( fit$marginals.fixed, FUN=summary_inv ) ) )

      O[["summary"]][["fixed_effects"]] = W [, tokeep, drop =FALSE]

      # hyperpar (variance components)   
      j = grep( "^Precision.*", rownames(fit$summary.hyperpar), value=TRUE )
      if (length(j) > 0) {
        precs = try( list_simplify( sapply( fit$marginals.hyperpar[j], FUN=summary_inv_prec ) ), silent=TRUE )  # prone to integration errors ..
        if (inherits(precs, "try-error")) precs = try( list_simplify( sapply( fit$marginals.hyperpar[j], FUN=summary_inv_prec_1024 ) ), silent=TRUE )
        if (inherits(precs, "try-error")) {
          message( "Try alt parameterization for prec -> sd or smaller number of n or masking negative values, using direct cionversion of summaries instead")
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

  if ("random_other" %in% toget) {

    if (exists("marginals.random", fit)) {

      raneff = names( fit$marginals.random )
      raneff = setdiff( raneff, c(vnS, vnST) ) 
      for (re in raneff) {
        g = fit$marginals.random[[re]]
        if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, mlimit=mlimit )
        O[["random"]] [[re]] = list_simplify ( sapply( g, summary_inv ) )  [, tokeep, drop =FALSE]
        O[["random"]] [[re]]$ID = fit$summary.random[[re]]$ID
      }
    }

  }

  if ("random_spatial" %in% toget) {
  
    if (exists("marginals.random", fit)) {

      if ( exists(vnS, fit$marginals.random) ) {
        
        O[["random"]] [[vnS]] = list()  # space as a main effect

        # bym2 output index locations
        Z = expand.grid( space=O[[vnSn]], type = c("iid", "bym"), stringsAsFactors =FALSE )

        u = which(Z$type=="bym")
        v = which(Z$type=="iid")
        
        g = fit$marginals.random[[vnS]]
        if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, mlimit=mlimit )

        # m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  
        m = list_simplify ( sapply( g, summary_inv ) ) 
        
        #  iid main effects  
        W = array( NA,  dim=c( length( O[[vnSn]]), length(names(m)) ), dimnames=list( space=O[[vnSn]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        for (k in 1:length(names(m))) {
          W[,k] = reformat_to_array( input = unlist(m[v,k]), matchfrom = list( space=Z[["space"]][v] ), matchto = list( space=O[[vnSn]]  ) )
        }
        O[["random"]] [[vnS]] [["iid"]] = W [, tokeep, drop =FALSE]

        #  spatial main effects  
        W = array( NA, dim=c( length( O[[vnSn]]), length(names(m)) ), dimnames=list( space=O[[vnSn]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        for (k in 1:length(names(m))) {
          W[,k] = reformat_to_array(  input = unlist(m[u,k]), matchfrom = list( space=Z[["space"]][u]  ), matchto = list( space=O[[vnSn]]  ) )
        } 
        O[["random"]] [[vnS]] [["bym"]] = W [, tokeep, drop =FALSE]

        # space == iid + bym combined:
        selection=list()
        selection[vnS] = 0  # 0 means everything matching space
        aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
        g = sapply( aa, function(x) invlink(x$latent[v] + x$latent[u]) )
        mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
        mm = apply( g, 1, mean, na.rm =TRUE)
        ms = apply( g, 1, sd, na.rm =TRUE)
        W = cbind(mm, ms, mq)
        attr(W, "dimnames") = list( space=O[[vnSn]], stat=tokeep  ) 
        O[["random"]] [[vnS]] [["combined"]] = W
        

        if (!is.null(exceedance_threshold)) {
          m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
          W = reformat_to_array( 
            input = m, 
            matchfrom = list( space=Z[["space"]] ), 
            matchto = list( space=O[[vnSn]]  )
          )
          names(dimnames(W))[1] = vnS 
          dimnames( W )[[vnS]] = O[[vnSn]]
          O[["random"]] [[vnS]] [["exceedance"]] = W
        }

        if (!is.null(deceedance_threshold)) {
          m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
          W = reformat_to_array( 
            input = m, 
              matchfrom = list( space=Z[["space"]] ), 
            matchto = list( space=O[[vnSn]]  )
          )
          names(dimnames(W))[1] = vnS 
          dimnames( W )[[vnS]] = O[[vnSn]]
          O[["random"]] [[vnS]] [["deceedance"]] = W
        }

      }
    }

  }


  if ("random_spatiotemporal"  %in% toget ) {

    if (exists("marginals.random", fit)) {

      if ( exists( vnST, fit$marginals.random) ) {

        O[["random"]] [[vnST]] = list()

        Z = expand.grid( space=O[[vnSn]], type = c("iid", "bym"), time=O[[vnTn]], stringsAsFactors =FALSE )
    
        u = which(Z$type=="bym")
        v = which(Z$type=="iid")
        
        g = fit$marginals.random[[vnST]]
        if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, mlimit=mlimit )

        # m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  
        m = list_simplify ( sapply( g, summary_inv ) ) 
    
        #  spatiotemporal interaction effects  bym
        W = array( NA, dim=c( length( O[[vnSn]]), length(O[[vnTn]]), length(names(m)) ), dimnames=list( space=O[[vnSn]], time=O[[vnTn]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..

        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array(  input = unlist(m[u,k]), matchfrom = list( space=Z[["space"]][u],  time=Z[["time"]][u]  ), matchto = list( space=O[[vnSn]], time=O[[vnTn]]  ) )
        } 
        O[["random"]] [[vnST]] [["bym"]] = W [,, tokeep, drop =FALSE]


        #  spatiotemporal interaction effects  iid
        W = array( NA, dim=c( length( O[[vnSn]]), length(O[[vnTn]]), length(names(m)) ), dimnames=list( space=O[[vnSn]], time=O[[vnTn]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array(  input = unlist(m[v,k]), matchfrom = list( space=Z[["space"]][v],  time=Z[["time"]][v]  ), matchto = list( space=O[[vnSn]], time=O[[vnTn]]  ) )
        } 
        O[["random"]] [[vnST]] [["iid"]] = W [,, tokeep, drop =FALSE]


        # space == iid + bym combined:
        selection=list()
        selection[vnST] = 0  # 0 means everything matching space
        aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  
        g = sapply( aa, function(x) invlink(x$latent[v] + x$latent[u]) )
        mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
        mm = apply( g, 1, mean, na.rm =TRUE)
        ms = apply( g, 1, sd, na.rm =TRUE)
        m = data.frame( cbind(mm, ms, mq) )
        names(m) = tokeep
        W = array( NA, dim=c( length( O[[vnSn]]), length(O[[vnTn]]), length(names(m)) ), dimnames=list( space=O[[vnSn]], time=O[[vnTn]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
        
        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array(  input = m[,k], matchfrom = list( space=Z[["space"]][v],  time=Z[["time"]][v]  ), matchto = list( space=O[[vnSn]], time=O[[vnTn]]  ) )
        } 
        O[["random"]] [[vnST]] [["combined"]] = W [,, tokeep, drop =FALSE]
        

        if (!is.null(exceedance_threshold)) {
          m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
          W = reformat_to_array( 
            input = m, matchfrom = list( space=Z[["space"]][v],  time=Z[["time"]][v]  ), matchto = list( space=O[[vnSn]], time=O[[vnTn]]  )
          )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnSn]]
          dimnames( W )[[vnT]] = O[[vnTn]]

          O[["random"]] [[vnST]] [["exceedance"]] = W
        }

        if (!is.null(deceedance_threshold)) {
          m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
          W = reformat_to_array( 
            input = m, matchfrom = list( space=Z[["space"]][v],  time=Z[["time"]][v]  ), matchto = list( space=O[[vnSn]], time=O[[vnTn]]  )
          )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnSn]]
          dimnames( W )[[vnT]] = O[[vnTn]]

          O[["random"]] [[vnST]] [["deceedance"]] = W
        }
      }
    }
  }


  if (any(grepl("predictions", toget))) {
   # revert these to character as they are used for matching
   
    SP = levels( O[[vnSf]] ) [ M[,vnS] ]
    TP = levels( O[[vnTf]] ) [ M[,vnT] ]

    i = which(  M$tag=="predictions" & SP %in% O[[vnSn]] & TP %in% O[[vnTn]] )  

  }
 

  if ("predictions"  %in% toget ) {
    
     if (exists("marginals.fitted.values", fit)) {

       g = fit$marginals.fitted.values[i]  # already on user scale
       if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, mlimit=invlink(mlimit) )

      m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  # already backtransformed by link=1
      W = array( NA,  dim=c( length( O[[vnSn]]), length(O[[vnTn]]), length(names(m)) ),  dimnames=list( space=O[[vnSn]], time=O[[vnTn]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
     
      
      for (k in 1:length(names(m))) {
        W[,,k] = reformat_to_array( 
          input = unlist(m[,k]), 
          matchfrom = list( space=SP[ i ], time=TP[i] ), 
          matchto =   list( space=O[[vnSn]], time=O[[vnTn]] )
        )
      } 
      O[["predictions"]] = W[,, tokeep, drop =FALSE]

    }

  }

  if ("predictions_adjusted"  %in% toget ) {
      # adjusted by offset
    if (exists("marginals.fitted.values", fit)) {

      g = fit$marginals.fitted.values  # g is already on user scale
      if (any( grepl( vnO, p$carstm_model_formula ))) {
        o = M[, vnO ]  # offsets; M[,vnO] is on the user scale .. it is log transformed in formula
        for (e in i) g[[e]][,1] = g[[e]][,1] - o[e] 
      }
      
      if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, mlimit=invlink(mlimit) )
      m = list_simplify ( sapply( g[i], inla.zmarginal, silent=TRUE ) )  # already backtransformed by link=1
      W = array( NA, dim=c( length(O[[vnSn]]), length(O[[vnTn]]), length(names(m)) ),  dimnames=list( space=O[[vnSn]], time=O[[vnTn]], stat=names(m) ) ) 
       names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..

      for (k in 1:length(names(m))) {
        W[,,k] = reformat_to_array( 
          input = unlist(m[,k]), 
          matchfrom = list( space=SP[ i ], time=TP[i] ), 
          matchto =   list( space=O[[vnSn]], time=O[[vnTn]] )
        )
      } 
      dimnames( W )[[vnS]] = O[[vnSn]]
      O[["predictions_adjusted"]] = W[,, tokeep, drop =FALSE]

      # offsets on user scale
     if (any( grepl( vnO, p$carstm_model_formula ))) {
        o =  M[, vnO ]   # offsets; M[,vnO] is on the user scale .. it is log transformed in formula .. so this is untransformed
        W = reformat_to_array( 
          input = o[i], 
          matchfrom = list( space=SP[ i ], time=TP[i] ), 
          matchto =   list( space=O[[vnSn]], time=O[[vnTn]] )
        )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnSn]]
          dimnames( W )[[vnT]] = O[[vnTn]]

        O[["prediction_offsets"]] = W

        O[["predictions_adjusted_direct"]] = O[["predictions"]][,,"mean"] /  O[["prediction_offsets"]]   ## relative rate as a check
     }


      if (!is.null(exceedance_threshold)) {
        m = list_simplify ( sapply( g[i], FUN=exceedance_prob, threshold=exceedance_threshold ) )
        W = reformat_to_array( 
          input = unlist(m ), 
          matchfrom = list( space=SP[ i ], time=TP[i] ), 
          matchto =   list( space=O[[vnSn]], time=O[[vnTn]] )
        )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnSn]]
          dimnames( W )[[vnT]] = O[[vnTn]]

          O[["exceedance_probability"]] = W
      }

      if (!is.null(deceedance_threshold)) {
        m = list_simplify ( sapply( g[i], FUN=deceedance_prob, threshold=deceedance_threshold ) )
        W = reformat_to_array( 
          input = unlist(m ), 
          matchfrom = list( space=SP[ i ], time=TP[i] ), 
          matchto =   list( space=O[[vnSn]], time=O[[vnTn]] )
        )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnSn]]
          dimnames( W )[[vnT]] = O[[vnTn]]

          O[["deceedance_probability"]] = W
      }

      if ( length(O[[vnTn]]) > 2 ) {
        if ( exists("predictions_adjusted", O ) ) {
          ti = as.numeric( O[[vnTn]] )  
          lmslope = function( x ) summary( lm(x~ti) )$coefficients["ti",1:2] 
          W = t ( apply( O[["predictions_adjusted"]][,,"mean"], 1, lmslope ) )  # relative rate per year
          names(dimnames(W))[1] = vnS 
            dimnames( W )[[vnS]] = O[[vnSn]]
               O[["time_slope"]] = W
        }
      }
    }
  }
 


  # copy data in case needed for plotting ..

  O[["M"]] = M
  O[["sppoly"]] = sppoly
  O[["fn_res"]] = fn_res 

  return(O)
 
}