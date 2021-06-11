

carstm_model_inla = function(p, M, fn_fit, toget="summary", file_compress_method=FALSE, exceedance_threshold=NULL, deceedance_threshold=NULL, nposteriors=NULL, improve.hyperparam.estimates=NULL, quantile_limit=NULL ) {

  # outputs
  O = list()
  O[["summary"]] = list()
  O[["random"]] = list()

  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (class(M)=="character") assign("M", eval(parse(text=M) ) )

  if (0) {
      # usual variable names used in aegis .. char / num / factor
      p$vnS = "space"  # "space"
      p$vnT = "time"
      p$vnU = "season" # "dyri"  # sub annual time 


      # alt character descrptions of vars (if any, taken from main above if not)
      p$vnS0 = "AUID"  # as character 
      p$vnT0 = "yr"  # as character 
      p$vnU0 = "dyear"
      
      # copies of main effects for inla in formulae
      p$vnST = "space_time"  # vnST = "space_time" (copy of vnS)
      p$vnTS = "time_space"  # vnTS = "time_space" (copy of vnT)
      
      # p$vnT0 = "time0"  # "time" as a factor
      # p$vnU0 = "season0"  # "time" as a factor

      # AUID is character; space is factor -> numeric 
  }


  # fiddling of AU and TU inputs: for bym2, etc, they need to be numeric, matching numerics of polygon id ("region.id")
  # set as factor to carry both name and index
  sppoly = areal_units( p=p )  # required by car fit
  region.id = as.character( slot( slot(sppoly, "nb"), "region.id" ) ) # the master / correct sequence of the AU's and neighbourhood matrix index values
  nAUID = nrow(sppoly)
  
  vnST = vnTS = NULL  # this is also used as a flag for random st effects extraction

  if (grepl("space", p$aegis_dimensionality)) {
    # labels
    vnS = ifelse( exists("vnS", p), p$vnS, "space" )
    vnS0 = ifelse( exists("vnS0", p), p$vnS0, "space0" ) # local storage/copy as a character

    O[[vnS]] = region.id  # the sequence is key as index matches nb matrix values
    O[[vnS]] = as.character( O[[vnS]] ) # local copy
    M[,vnS0] = as.character( M[,vnS] ) # local copy
    M[,vnS] = match( M[,vnS0], O[[vnS]] )  # overwrite with numeric values that must match index of neighbourhood matrix
  }

  if (grepl("year", p$aegis_dimensionality)) {
    vnT = ifelse( exists("vnT", p), p$vnT, "time" )
    vnT0 = ifelse( exists("vnT0", p), p$vnT0, "time0" )

    O[[vnT]] = as.character( p$yrs )
    M[,vnT0] = as.character( M[,vnT] ) # copy
    M[,vnT] = match( M[,vnT0], O[[vnT]] )

      # internal vars, for inla
    vnST = ifelse( exists("vnST", p), p$vnST, "space_time" )
    vnTS = ifelse( exists("vnTS", p), p$vnTS, "time_space" )
    M[,vnST] = M[,vnS]   
    M[,vnTS] = M[,vnT] 
  }

  if (grepl("season", p$aegis_dimensionality)) {
    vnU = ifelse( exists("vnU", p), p$vnU, "season" )  # sub-annual time
    vnU0 = ifelse( exists("vnU0", p), p$vnU0, "season0" )

    O[[vnU]] = as.character( p$dyears + diff(p$dyears)[1]/2)   
    M[,vnU0] = as.character( M[,vnU] )
    M[,vnU] = match( M[,vnU0], O[[vnU]] )
  }
 

  if (is.null(nposteriors))  nposteriors = ifelse( exists("nposteriors", p), p$nposteriors, 1000 )

  vnY = p$variabletomodel
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
  m = M[ j, vnY ]
  
  vnO = ifelse( exists("vnO", p), p$vnO, "data_offset" )  # input data is expected to be on user scale
  if (any( grepl( vnO, p$carstm_model_formula )))  {
    if (exists( vnO, M)) m = m / M[ j, vnO ]
  } 
  

  # on user scale
  ii = which(is.finite(m))
  if (!is.null(quantile_limit)) upper_limit = quantile( m[ii], probs=quantile_limit ) 
  mq = quantile( m[ii], probs=c(0.025, 0.5, 0.975) )
  
  O[["data_range"]] = c( mean=mean(m[ii]), sd=sd(m[ii]), min=min(m[ii]), max=max(m[ii]),  lb=mq[1], median=mq[2], ub=mq[3]  )  # on data /user scale not internal link 
  
  # on link scale:
  ml = lnk_function( m ) # necessary in case of log(0)
  ll = which(is.finite(ml))
  H = carstm_hyperparameters( sd(ml[ll] ), alpha=0.5, median(ml[ll], na.rm=TRUE) )  # sd slightly biased due to 0's being dropped
  m = ml = ii = ll = j = NULL
  gc()

  p = parameters_add_without_overwriting( p, options.control.family = inla.set.control.family.default() )

  fit  = NULL

  if (!exists("options.control.inla", p )) p$options.control.inla = list(
    inla.set.control.inla.default(),  # first try defaults
    list( strategy="adaptive", improved.simplified.laplace=TRUE ), # default h=0.005
    list( stupid.search=FALSE, strategy="adaptive", h=0.05, cmin=0, tolerance=1e-9),
    list( stupid.search=FALSE, strategy="adaptive", h=0.1, cmin=0),
    list( stupid.search=FALSE, strategy="adaptive", h=0.001, cmin=0), # default h=0.005  
    list( stupid.search=TRUE, strategy="adaptive", h=0.2, cmin=0, optimiser="gsl" ), # default h=0.005
    list( stupid.search=TRUE, fast=FALSE, step.factor=0.1),
    list( stupid.search=TRUE, cmin=0, optimiser="gsl" )
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
    message("If you are using MSWindows and you get a popup complaining about 'inla stopped working',")
    message("try setting the flag in the following link to 1, using regedit. Be careful.")
    message("e.g., see: https://monitormyweb.com/guides/how-to-disable-stopped-working-message-in-windows")
    stop( "solution did not converge")
  }

  # to improve hyper param estimates..
  if (!is.null(improve.hyperparam.estimates)) if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.6, diff.logdens=9  )  # get improved estimates for the hyperparameters

  if (is.null(fit)) warning("model fit error")
  if ("try-error" %in% class(fit) ) warning("model fit error")

  message( "Saving carstm fit (this can be slow): ", fn_fit )

  save( fit, file=fn_fit, compress=file_compress_method )


  # do the computations here as fit can be massive ... best not to copy, etc ..
  message( "Computing summaries ..." )




  invlink = function(x) lnk_function( x,  inverse=TRUE ) 

  summary_inv = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )

  summary_inv_prec = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(pmax(y,1e-12)), x) , silent=TRUE  )
  # summary_inv_prec_1024 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=1024L) , silent=TRUE  )
  # summary_inv_prec_512 = function(x) inla.zmarginal( inla.tmarginal( function(y) 1/sqrt(y), x, n=512L) , silent=TRUE  )

  list_simplify = function(x) as.data.frame( t( as.data.frame( x )))

  exceedance_prob = function(x, threshold)  {1 - inla.pmarginal(q = threshold, x)}

  deceedance_prob = function(x, threshold)  { inla.pmarginal(q = threshold, x)}

  truncate_upperbound = function( b, upper_limit, eps=1e-12 ) {
    k = which( b[,1] > upper_limit )  
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
        O[["random"]] [[re]] = list_simplify ( sapply( g, summary_inv ) )  [, tokeep, drop =FALSE]
        O[["random"]] [[re]]$ID = fit$summary.random[[re]]$ID
      }
    }

  }

  if ("random_spatial" %in% toget) {
    # space only
    
    if (exists("marginals.random", fit)) {

      if ( exists(vnS, fit$marginals.random) ) {
        
        O[["random"]] [[vnS]] = list()  # space as a main effect

        n_S = length(fit$marginals.random[[vnS]]) 
        iid = NULL
        bym = NULL
        if ( n_S == nAUID) {
          # single spatial effect (eg besag)
          Z = expand.grid( space=O[[vnS]], type = c("iid"), stringsAsFactors =FALSE )
          iid =  which(Z$type=="iid")
        }

        if ( n_S == nAUID*2) {
          # bym2 effect: bym and iid
          Z = expand.grid( space=O[[vnS]], type = c("iid", "bym"), stringsAsFactors =FALSE )
          bym = which(Z$type=="bym")
          iid = which(Z$type=="iid")
        }

        g = fit$marginals.random[[vnS]]
        # m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  
        m = list_simplify ( sapply( g, summary_inv ) ) 
        
        if (!is.null(iid)) {
          #  iid main effects  
          W = array( NA,  dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) ) 
          names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
          for (k in 1:length(names(m))) {
            W[,k] = reformat_to_array( input = unlist(m[iid,k]), matchfrom = list( space=Z[["space"]][iid] ), matchto = list( space=O[[vnS]] ) )
          }
          O[["random"]] [[vnS]] [["iid"]] = W [, tokeep, drop =FALSE]
        }

        if (!is.null(bym)) {
          #  spatial main effects  
          W = array( NA, dim=c( length( O[[vnS]]), length(names(m)) ), dimnames=list( space=O[[vnS]], stat=names(m) ) ) 
          names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
          for (k in 1:length(names(m))) {
            W[,k] = reformat_to_array(  input = unlist(m[bym,k]), matchfrom = list( space=Z[["space"]][bym]  ), matchto = list( space=O[[vnS]] ) )
          } 
          O[["random"]] [[vnS]] [["bym"]] = W [, tokeep, drop =FALSE]
        }

        if (!is.null(iid)  & !is.null(bym) ) {
          # space == iid + bym combined:
          selection=list()
          selection[vnS] = 0  # 0 means everything matching space
          aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  # 0 means everything matching space
          g = sapply( aa, function(x) invlink(x$latent[iid] + x$latent[bym]) )
          mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
          mm = apply( g, 1, mean, na.rm =TRUE)
          ms = apply( g, 1, sd, na.rm =TRUE)
          W = cbind(mm, ms, mq)
          attr(W, "dimnames") = list( space=O[[vnS]], stat=tokeep  ) 
          O[["random"]] [[vnS]] [["combined"]] = W
        }

        if (!is.null(exceedance_threshold)) {
          m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
          W = reformat_to_array( input = m, matchfrom = list( space=Z[["space"]][iid] ), matchto = list( space=O[[vnS]] ) )
          names(dimnames(W))[1] = vnS           
          dimnames( W )[[vnS]] = O[[vnS]]
          O[["random"]] [[vnS]] [["exceedance"]] = W
        }

        if (!is.null(deceedance_threshold)) {
          m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
          W = reformat_to_array( input = m, matchfrom = list( space=Z[["space"]][iid] ), matchto = list( space=O[[vnS]] )  )
          names(dimnames(W))[1] = vnS           
          dimnames( W )[[vnS]] = O[[vnS]]
          O[["random"]] [[vnS]] [["deceedance"]] = W
        }

      }
    }

  }


  if ("random_spatiotemporal"  %in% toget ) {
    # space-year

    if (exists("marginals.random", fit)) {

      if (!is.null(vnST)) {

        if ( exists( vnST, fit$marginals.random) ) {

          O[["random"]] [[vnST]] = list()

          n_ST = length(fit$marginals.random[[vnST]]) 
          iid = NULL
          bym = NULL
    
          if ( n_ST == nAUID* p$ny) {
            # besag effect: with annual results
            Z = expand.grid( space=O[[vnS]], type = c("iid"), time=O[[vnT]], stringsAsFactors =FALSE )
            iid =  which(Z$type=="iid")
          }

          if ( n_ST == nAUID*2 * p$ny) {
            # bym2 effect: bym and iid with annual results
            Z = expand.grid( space=O[[vnS]], type = c("iid", "bym"), time=O[[vnT]], stringsAsFactors =FALSE )
            bym = which(Z$type=="bym")
            iid = which(Z$type=="iid")
          }
  
          g = fit$marginals.random[[vnST]]
          # m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  
          m = list_simplify ( sapply( g, summary_inv ) ) 

          if (!is.null(iid)) {
            #  spatiotemporal interaction effects  iid
            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) ) 
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = unlist(m[iid,k]), matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  ), matchto = list( space=O[[vnS]], time=O[[vnT]]  ) )
            } 
            O[["random"]] [[vnST]] [["iid"]] = W [,, tokeep, drop =FALSE]
          }

          if (!is.null(bym)) {
            #  spatiotemporal interaction effects  bym
            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) ) 
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..

            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = unlist(m[bym,k]), matchfrom = list( space=Z[["space"]][bym],  time=Z[["time"]][bym]  ), matchto = list( space=O[[vnS]], time=O[[vnT]]  ) )
            } 
            O[["random"]] [[vnST]] [["bym"]] = W [,, tokeep, drop =FALSE]
          }

          if (!is.null(iid)  & !is.null(bym) ) {
            # space == iid + bym combined:
            selection=list()
            selection[vnST] = 0  # 0 means everything matching space
            aa = inla.posterior.sample( nposteriors, fit, selection=selection, add.names =FALSE )  
            g = sapply( aa, function(x) invlink(x$latent[iid] + x$latent[bym]) )
            mq = t( apply( g, 1, quantile, probs =c(0.025, 0.5, 0.975), na.rm =TRUE) )
            mm = apply( g, 1, mean, na.rm =TRUE)
            ms = apply( g, 1, sd, na.rm =TRUE)
            m = data.frame( cbind(mm, ms, mq) )
            names(m) = tokeep
            W = array( NA, dim=c( length( O[[vnS]]), length(O[[vnT]]), length(names(m)) ), dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) ) 
            names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
            names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
            
            for (k in 1:length(names(m))) {
              W[,,k] = reformat_to_array(  input = m[,k], matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  ), matchto = list( space=O[[vnS]], time=O[[vnT]]  ) )
            } 
            O[["random"]] [[vnST]] [["combined"]] = W [,, tokeep, drop =FALSE]
          }         

          if (!is.null(exceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x > exceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( 
              input = m, matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  ), matchto = list( space=O[[vnS]], time=O[[vnT]]  )
            )
            names(dimnames(W))[1] = vnS 
            names(dimnames(W))[2] = vnT
            dimnames( W )[[vnS]] = O[[vnS]]
            dimnames( W )[[vnT]] = O[[vnT]]
            O[["random"]] [[vnST]] [["exceedance"]] = W
          }

          if (!is.null(deceedance_threshold)) {
            m = apply ( g, 1, FUN=function(x) length( which(x < deceedance_threshold) ) ) / nposteriors
            W = reformat_to_array( 
              input = m, matchfrom = list( space=Z[["space"]][iid],  time=Z[["time"]][iid]  ), matchto = list( space=O[[vnS]], time=O[[vnT]]  )
            )
            names(dimnames(W))[1] = vnS 
            names(dimnames(W))[2] = vnT
            dimnames( W )[[vnS]] = O[[vnS]]
            dimnames( W )[[vnT]] = O[[vnT]]

            O[["random"]] [[vnST]] [["deceedance"]] = W
          }
        }  # ned space-year
      }
    }
  }

   

  if ("predictions"  %in% toget ) {
      # adjusted by offset
    if (exists("marginals.fitted.values", fit)) {
 

      if (  "space"==p$aegis_dimensionality ) {
        ipred = which( M$tag=="predictions"  &  M[,vnS0] %in% O[[vnS]] )  # filter by S and T in case additional data in other areas and times are used in the input data

        g = fit$marginals.fitted.values[ipred]  # g is already on user scale
        if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, upper_limit=upper_limit )
        m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  # already backtransformed by link=1

        W = array( NA, dim=c( length(O[[vnS]]),  length(names(m)) ),  dimnames=list( space=O[[vnS]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..

        for (k in 1:length(names(m))) {
          W[,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=list( space=M[ ipred, vnS0] ), matchto=list( space=O[[vnS]] ))
        } 
        O[["predictions"]] = W[,, tokeep, drop =FALSE]
      }

      if (  "space-year"==p$aegis_dimensionality ) {
        ipred = which( M$tag=="predictions" & M[,vnS0] %in% O[[vnS]] & M[,vnT0] %in% O[[vnT]] )  
        g = fit$marginals.fitted.values[ipred]  # g is already on user scale
        if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, upper_limit=upper_limit )
        m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  # already backtransformed by link=1

        W = array( NA, dim=c( length(O[[vnS]]), length(O[[vnT]]), length(names(m)) ),  dimnames=list( space=O[[vnS]], time=O[[vnT]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..

        for (k in 1:length(names(m))) {
          W[,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=list( space=M[ipred,vnS0], time=M[ipred,vnT0] ), matchto= list( space=O[[vnS]], time=O[[vnT]] ))
        } 
        O[["predictions"]] = W[,, tokeep, drop =FALSE]

      }

      if (  "space-year-season"==p$aegis_dimensionality ) {
        ipred = which( M$tag=="predictions" & M[,vnS0] %in% O[[vnS]]  &  M[,vnT0] %in% O[[vnT]] &  M[,vnU0] %in% O[[vnU]])  # ignoring U == predict at all seassonal components ..
        g = fit$marginals.fitted.values[ipred]  # g is already on user scale
        if (!is.null(quantile_limit)) g = lapply( g, truncate_upperbound, upper_limit=upper_limit )
        m = list_simplify ( sapply( g, inla.zmarginal, silent=TRUE ) )  # already backtransformed by link=1

        W = array( NA, dim=c( length(O[[vnS]]), length(O[[vnT]]), length(O[[vnU]]), length(names(m)) ),  dimnames=list( space=O[[vnS]], time=O[[vnT]], season=O[[vnU]], stat=names(m) ) ) 
        names(dimnames(W))[1] = vnS  # need to do this in a separate step ..
        names(dimnames(W))[2] = vnT  # need to do this in a separate step ..
        names(dimnames(W))[3] = vnU  # need to do this in a separate step ..

        for (k in 1:length(names(m))) {
          W[,,,k] = reformat_to_array( input=unlist(m[,k]), matchfrom=list( space=M[ipred, vnS0], time=M[ipred, vnT0], season=M[ipred, vnU0] ), matchto=list( space=O[[vnS]], time=O[[vnT]], season=O[[vnU]] ))
        } 
        O[["predictions"]] = W[,, tokeep, drop =FALSE]
      }


      if (exists("data_transformation", p) ) O[["predictions"]] = p$data_transformation$backward( O[["predictions"]] ) # make all positive


      if (!is.null(exceedance_threshold)) {
        m = list_simplify ( sapply( g, FUN=exceedance_prob, threshold=exceedance_threshold ) )
        W = reformat_to_array( 
          input = unlist(m ), 
          matchfrom = list( space=M[ipred,vnS], time=M[ipred,vnT]  ), 
          matchto =   m_to
        )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnS]]
          dimnames( W )[[vnT]] = O[[vnT]]

          O[["exceedance_probability"]] = W
      }

      if (!is.null(deceedance_threshold)) {
        m = list_simplify ( sapply( g, FUN=deceedance_prob, threshold=deceedance_threshold ) )
        W = reformat_to_array( 
          input = unlist(m ), 
          matchfrom = list( space=M[ipred,vnS], time=M[ipred,vnT] ), 
          matchto =   m_to
        )
          names(dimnames(W))[1] = vnS 
          names(dimnames(W))[2] = vnT
          dimnames( W )[[vnS]] = O[[vnS]]
          dimnames( W )[[vnT]] = O[[vnT]]

          O[["deceedance_probability"]] = W
      }

      if ( length(O[[vnT]]) > 2 ) {
        if ( exists("predictions", O ) ) {
          ti = as.numeric( O[[vnT]] )  
          lmslope = function( x ) summary( lm(x~ti) )$coefficients["ti",1:2] 
          W = t ( apply( O[["predictions"]][,,"mean"], 1, lmslope ) )  # relative rate per year
          names(dimnames(W))[1] = vnS 
            dimnames( W )[[vnS]] = O[[vnS]]
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