
carstm_posterior_predictions_simulations = function(
  sppoly =NULL,
  fit = NULL,
  space_id=NULL, time_id=NULL, cyclic_id=NULL, 
  fn_fit=NULL, 
  num.threads = "1:1"
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
    
  if (is.null(fit)) {
    if (file.exists(fn_fit)) {
      if (grepl("\\.RDS$", fn_fit)) {
        fit = readRDS(fn_fit)
      } else {
        load( fn_fit )
      }
    }
  }

  
 
  inla.setOption(num.threads=num.threads)
  
  mc.cores = as.numeric( unlist(strsplit(num.threads, ":") )[1] )

  # do the computations here as fit can be massive ... best not to copy, etc ..
  if (be_verbose)  message( "\nComputing summaries and extracting posterior simulations ..." )

  if (!exists("modelinfo", fit)) stop("modelinfo not in fit, fit object needs to be re-run")
  O = fit$modelinfo
  
  fit$modelinfo = NULL
  gc()


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
      Z = try( apply_generic( Z, inla.tmarginal, fun=invlink) )
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
 
    S = inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores ) 

    message( "Sampling complete ... now reformatting and extracting required components" )

    for (z in c("tag", "start", "length") ) assign(z, attributes(S)[[".contents"]][[z]] )  # index info

      if (0) {
        # not used at the moment .. 
        nlogdens = length(S[[1]]$logdens)
        logdens = array(NA, dim=c( nlogdens, nposteriors  ) )
        for (i in 1:nposteriors) {
          logdens[,i] = unlist(S[[i]]$logdens)
        }
        logdens_names =  names(S[[1]]$logdens)
        logdens = format_results( logdens, labels=logdens_names  )# same seq as space_id ( == attributes(space)$row.names )
      }
  }

  names.fixed = fit$names.fixed
  summary.fixed  = fit$summary.fixed
  summary.random = fit$summary.random
  summary.hyperpar = fit$summary.hyperpar 
  
  dic = fit$dic[c("dic", "p.eff", "dic.sat", "mean.deviance")]
  waic = fit$waic[c("waic", "p.eff")]
  mlik = fit$mlik[2]


  
  # separate out random spatial and randomm spatiotemporal (as they can be large arrays)
  if ("random_spatial" %in% toget) {
    # space only
    Z = NULL
    iSP = which( re$dimensionality=="s" & re$level=="main")
    if (length(iSP) > 0 ) {

      if (be_verbose)  message("Extracting random spatial errors"  )
      if (exists("debug")) if (is.character(debug)) if ( debug =="random_spatial") browser()

      matchto = list( space=O[["space_id"]] )

      W = array( NA, dim=c( O[["space_n"]], length(tokeep) ), dimnames=list( space=O[["space_id"]], stat=tokeep ) )
      names(dimnames(W))[1] = vS  # need to do this in a separate step ..

      O[["random"]] [[vS]] = list()  # space as a main effect  vS==vS

      if (length(iSP) == 1) {

        # vS = re$vn[ iSP ]  # == vS as it is a single spatial effect
        
        model_name = re$model[ iSP ]  # should be iid

        m = fit$marginals.random[[vS]]
 
        m = try( apply_generic( m, inla.tmarginal, fun=invlink) , silent=TRUE )
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
         
          m = try( apply_generic( m, inla.tmarginal, fun=invlink), silent=TRUE  )
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
        slabels = O[["space_id"]]

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
            dimnames( m )[[vS]] = O[["space_id"]]
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
            dimnames( m )[[vS]] = O[["space_id"]]
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
        dimnames=list( space=O[["space_id"]], time=O[["time_id"]], stat=tokeep ) )
      names(dimnames(W))[1] = vS  # need to do this in a separate step ..
      names(dimnames(W))[2] = vT  # need to do this in a separate step ..

      if (length(iST) == 1) {

        vnST = re$vn[ iST ]
        model_name = re$model[ iST ]   

        if (exists(vnST, fit$marginals.random )) {
 
          m = fit$marginals.random[[vnST]]
          m = try( apply_generic( m, inla.tmarginal, fun=invlink) , silent=TRUE )

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
          m = try( apply_generic( m, inla.tmarginal, fun=invlink), silent=TRUE )
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
            dimnames( m )[[vS]] = O[["space_id"]]
            dimnames( m )[[vT]] = O[["time_id"]]
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
            names(dimnames(m))[2] = O[["time_id"]]
            dimnames( m )[[vS]] = O[["space_id"]]
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
    npredictions = diff(O[["predictions_range"]])+1 
 
    preds = O[["predictions_range"]][1]:O[["predictions_range"]][2]  
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
          dimnames=list( space=O[["space_id"]], stat=names(m) ) )
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
            dimnames=list( space=O[["space_id"]], sim=1:nposteriors ) )
  
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
          dimnames=list( space=O[["space_id"]], time=O[["time_id"]], stat=names(m) ) 
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
            dimnames=list( space=O[["space_id"]], time=O[["time_id"]], sim=1:nposteriors ) )

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
          dimnames=list( space=O[["space_id"]], time=O[["time_id"]], cyclic=O[["cyclic_id"]], stat=names(m) ) )
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
            dimnames=list( space=O[["space_id"]], time=O[["time_id"]], cyclic=O[["cyclic_id"]], sim=1:nposteriors ) )

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
              dimnames( W )[[vS]] = O[["space_id"]]
              m = NULL
              if ( length(O[["matchfrom"]]) > 1 ) {
                names(dimnames(W))[2] = vT
                dimnames( W )[[vT]] = O[["time_id"]]
              }
              if (length(O[["matchfrom"]]) > 2  ) {
                names(dimnames(W))[3] = vU
                dimnames( W )[[vU]] = O[["cyclic_id"]]
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
              dimnames( W )[[vS]] = O[["space_id"]]
              m = NULL
              if ( length(O[["matchfrom"]]) > 1 ) {
                names(dimnames(W))[2] = vT
                dimnames( W )[[vT]] = O[["time_id"]]
              }
              if ( length(O[["matchfrom"]]) > 2 ) {
                names(dimnames(W))[3] = vU
                dimnames( W )[[vU]] = O[["cyclic_id"]]
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
