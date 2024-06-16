
carstm_posterior_simulations = function( p=NULL, pN=NULL, pW=NULL, pH=NULL, pB=NULL, qmax=NULL,
  wgts_min=0, wgts_max=NULL, N_min=0, N_max=NULL, B_min=0, B_max=NULL, max_value=NULL, degree_day=FALSE, 
  pa_threshold=0.05, hurdle_direct=FALSE, denoise="", carstm_directory=NULL ) {
  
  # N and B are actually densities
  
  operation = NULL
  
  if (!is.null(p)) operation = c( operation, "generic" )
  if (!is.null(pN)) operation = c( operation, "number" )
  if (!is.null(pW)) operation = c( operation, "meansize" )
  if (!is.null(pH)) operation = c( operation, "presence_absence" )
  if (!is.null(pB)) operation = c( operation, "biomass" )
 
  if (is.null(operation)) {
    stop( "The operation is ambiguous, check parameters." )
  }

  gen = NULL
  nums = NULL
  biom = NULL
  pa = NULL
  wgts = NULL

  if ( "generic" %in% operation) {
    # simple wrapper to data .. direct access would give more control if a single set is wanted .. 
    # this is more useful when complex transformations are required using multiple sims (e.g. 2 or 3-operations below)
    gen = carstm_model( p=pW, DS="carstm_modelled_summary", carstm_directory=carstm_directory  )

    if ( denoise == "" ) {
      gen =  gen[["sims"]][["predictions"]]
    } else {
      if (grepl( "spatial", denoise)) vnS = gen[["fm"]][["vn"]][["S"]]
      if (grepl( "spatiotemporal", denoise)) vnST = gen[["fm"]][["vn"]][["S2"]]
      if ("spatial" %in% denoise )  gen = gen[["sims"]][["predictions"]] - gen[["sims"]] [[vnS]] [["re_total"]] 
      if ("spatial_iid" %in% denoise )  gen = gen[["sims"]][["predictions"]] - gen[["sims"]] [[vnS]] [["re_unstructured"]]
      if ("spatial_bym2" %in% denoise )  gen = gen[["sims"]][["predictions"]] - gen[["sims"]] [[vnS]] [["re_neighbourhood"]]
      if ("spatiotemporal" %in% denoise ) gen = gen[["sims"]][["predictions"]] - gen[["sims"]] [[vnST]] [["re_total"]] 
      if ("spatiotemporal_iid" %in% denoise ) gen = gen[["sims"]][["predictions"]] - gen[["sims"]] [[vnST]] [["re_unstructured"]]
      if ("spatiotemporal_bym2" %in% denoise ) gen = gen[["sims"]][["predictions"]] - gen[["sims"]] [[vnST]] [["re_neighbourhood"]]
    }

    gen[!is.finite(gen)] = NA
    if ( !is.null(min_value) ) {
      i = which( gen < min_value )
      if (length(i) > 0 ) gen[ i ] = min_value
    }
    if ( !is.null(max_value) ) {
      i = which( gen > max_value )
      if (length(i) > 0 ) gen[ i ] = max_value
    }

    if (degree_day) {
      # compute degree days from sim ... todo 
      # collapse seasonal cycle into a cummulative distribution

    }

    if (length(operation)==1) return(gen)
    if (length(operation)>1) stop( "More than one generic operation specified ... carstm_posterior_simulations needs more rules ...")
  }


  if ( "presence_absence" %in% operation ) {
    pa = carstm_model( p=pH, DS="carstm_modelled_summary", carstm_directory=carstm_directory    )
    if ( denoise == "" ) {
      pa = pa[["sims"]][["predictions"]]
    } else {
      if (grepl( "spatial", denoise)) vnS = pa[["fm"]][["vn"]][["S"]]
      if (grepl( "spatiotemporal", denoise)) vnST = pa[["fm"]][["vn"]][["S2"]]
      if ("spatial" %in% denoise )  pa = pa[["sims"]][["predictions"]] - pa[["sims"]] [[vnS]] [["re_total"]] 
      if ("spatial_iid" %in% denoise )  pa = pa[["sims"]][["predictions"]] - pa[["sims"]] [[vnS]] [["re_unstructured"]]
      if ("spatial_bym2" %in% denoise )  pa = pa[["sims"]][["predictions"]] - pa[["sims"]] [[vnS]] [["re_neighbourhood"]]
      if ("spatiotemporal" %in% denoise ) pa = pa[["sims"]][["predictions"]] - pa[["sims"]] [[vnST]] [["re_total"]] 
      if ("spatiotemporal_iid" %in% denoise ) pa = pa[["sims"]][["predictions"]] - pa[["sims"]] [[vnST]] [["re_unstructured"]]
      if ("spatiotemporal_bym2" %in% denoise ) pa = pa[["sims"]][["predictions"]] - pa[["sims"]] [[vnST]] [["re_neighbourhood"]]
    }

    j = which( !is.finite(pa) )
    if (length(j) > 0 ) pa[j] = NA
    # pa = inverse.logit(pa)
    attr( pa, "unit") = "probability"
    if (length(operation)==1) return(pa)
  }


  # construct meansizes matrix used to convert number to weight
  if ( "meansize" %in% operation) {
    wgts = carstm_model( p=pW, DS="carstm_modelled_summary", carstm_directory=carstm_directory  )
    if ( denoise == "" ) {
      wgts = wgts[["sims"]][["predictions"]]
    } else {
      if (grepl( "spatial", denoise)) vnS = wgts[["fm"]][["vn"]][["S"]]
      if (grepl( "spatiotemporal", denoise)) vnST = wgts[["fm"]][["vn"]][["S2"]]
      if ("spatial" %in% denoise )  wgts = wgts[["sims"]][["predictions"]] - wgts[["sims"]] [[vnS]] [["re_total"]] 
      if ("spatial_iid" %in% denoise )  wgts = wgts[["sims"]][["predictions"]] - wgts[["sims"]] [[vnS]] [["re_unstructured"]]
      if ("spatial_bym2" %in% denoise )  wgts = wgts[["sims"]][["predictions"]] - wgts[["sims"]] [[vnS]] [["re_neighbourhood"]]
      if ("spatiotemporal" %in% denoise ) wgts = wgts[["sims"]][["predictions"]] - wgts[["sims"]] [[vnST]] [["re_total"]] 
      if ("spatiotemporal_iid" %in% denoise ) wgts = wgts[["sims"]][["predictions"]] - wgts[["sims"]] [[vnST]] [["re_unstructured"]]
      if ("spatiotemporal_bym2" %in% denoise ) wgts = wgts[["sims"]][["predictions"]] - wgts[["sims"]] [[vnST]] [["re_neighbourhood"]]
    }

    j = which( !is.finite(wgts) )
    if (length(j) > 0 ) wgts[j] = NA

    if (!is.null(qmax)) if ( is.null(wgts_max) )  wgts_max = quantile( wgts, probs=qmax, na.rm=TRUE )

    if ( !is.null(wgts_min) ) {
      i = which( wgts < wgts_min )
      if (length(i) > 0 ) wgts[ i ] = wgts_min
    }
    if ( !is.null(wgts_max) ) {
      i = which( wgts > wgts_max )
      if (length(i) > 0 ) wgts[ i ] = wgts_max
    }
    

    attr( wgts, "unit") = "kg"
    if (length(operation)==1) return(wgts)
  }


  if ( "biomass" %in% operation) {
    biom = carstm_model( p=pB, DS="carstm_modelled_summary", carstm_directory=carstm_directory   )
    if ( denoise == "" ) {
      biom = biom[["sims"]][["predictions"]]
    } else {
      if (grepl( "spatial", denoise)) vnS = biom[["fm"]][["vn"]][["S"]]
      if (grepl( "spatiotemporal", denoise)) vnST = biom[["fm"]][["vn"]][["S2"]]
      if ("spatial" %in% denoise )  biom = biom[["sims"]][["predictions"]] - biom[["sims"]] [[vnS]] [["re_total"]] 
      if ("spatial_iid" %in% denoise )  biom = biom[["sims"]][["predictions"]] - biom[["sims"]] [[vnS]] [["re_unstructured"]]
      if ("spatial_bym2" %in% denoise )  biom = biom[["sims"]][["predictions"]] - biom[["sims"]] [[vnS]] [["re_neighbourhood_neighbourhood"]]
      if ("spatiotemporal" %in% denoise ) biom = biom[["sims"]][["predictions"]] - biom[["sims"]] [[vnST]] [["re_total"]] 
      if ("spatiotemporal_iid" %in% denoise ) biom = biom[["sims"]][["predictions"]] - biom[["sims"]] [[vnST]] [["re_unstructured"]]
      if ("spatiotemporal_bym2" %in% denoise ) biom = biom[["sims"]][["predictions"]] - biom[["sims"]] [[vnST]] [["re_neighbourhood"]]
    }

    j = which( !is.finite(biom) )
    if (length(j) > 0 ) biom[j] = NA

    if (!is.null(qmax)) if ( is.null(B_max) ) B_max = quantile( biom, probs=qmax, na.rm=TRUE )

    if ( !is.null( B_min ) ) {
      i = which( biom < B_min )
      if (length(i) > 0 ) biom[ i ] = B_min
    }
    if ( !is.null( B_max ) ) {
      i = which( biom > B_max )
      if (length(i) > 0 ) biom[ i ] = B_max
    }
    # biom = biom / 10^6  # kg / km^2 -> kt / km^2
    attr(biom, "unit") = "kg/km^2"
    if (length(operation)==1) return(biom)
  }


  if ( "number" %in% operation ) {
    nums = carstm_model( p=pN, DS="carstm_modelled_summary", carstm_directory=carstm_directory  )
        
    if ( denoise == "" ) {
      nums = nums[["sims"]][["predictions"]]
    } else {
      if (grepl( "spatial", denoise)) vnS = nums[["fm"]][["vn"]][["S"]]
      if (grepl( "spatiotemporal", denoise)) vnST = nums[["fm"]][["vn"]][["S2"]]
      if ("spatial" %in% denoise )  nums = nums[["sims"]][["predictions"]] - nums[["sims"]] [[vnS]] [["re_total"]] 
      if ("spatial_iid" %in% denoise )  nums = nums[["sims"]][["predictions"]] - nums[["sims"]] [[vnS]] [["re_unstructured"]]
      if ("spatial_bym2" %in% denoise )  nums = nums[["sims"]][["predictions"]] - nums[["sims"]] [[vnS]] [["re_neighbourhood"]]
      if ("spatiotemporal" %in% denoise ) nums = nums[["sims"]][["predictions"]] - nums[["sims"]] [[vnST]] [["re_total"]] 
      if ("spatiotemporal_iid" %in% denoise ) nums = nums[["sims"]][["predictions"]] - nums[["sims"]] [[vnST]] [["re_unstructured"]]
      if ("spatiotemporal_bym2" %in% denoise ) nums = nums[["sims"]][["predictions"]] - nums[["sims"]] [[vnST]] [["re_neighbourhood"]]
    }

    j = which( !is.finite(nums) )
    if (length(j) > 0 ) nums[j] = NA

    if (!is.null(qmax)) if ( is.null(N_max) )  N_max = quantile( nums, probs=qmax, na.rm=TRUE )

    if ( !is.null(N_min) ) {
      i = which( nums < N_min )
      if (length(i) > 0 ) nums[ i ] = N_min
    }
    if ( !is.null(N_max) ) {
      i = which( nums > N_max )
      if (length(i) > 0 ) nums[ i ] = N_max
    }
    attr(nums, "unit") = "n/km^2"
    if (length(operation)==1) return(nums)
  }

 
  if (length(operation) == 2 ) {
    if ("meansize" %in% operation ) {
      if ( "biomass" %in% operation ) {
        nums = biom / wgts   # (n * 10^6) / km^2
        attr(nums, "unit") = "n/km^2"
        return(nums)
      }
      if ( "number" %in% operation ) {
        # numerical: density no / m^2  -->>  (no. km^2)
        biom = nums * wgts # * 10^6 / 10^6  # cancels out .. kg / km^2 -> kt / km^2
        attr(biom, "unit") = "kg/km^2"
        return(biom)
      }
    } 
    

    # NOTE: r = binom0 / pois0 is the Hurdle correction for the truncation (censored mean)
    # binom0 = pa  # prob > pa_threshold
    # pois0 = ppois( pa_threshold, lambda = mu, lower.tail = FALSE)  # prob > pa_threshold
 
    if ("presence_absence" %in% operation) {
      if ( "biomass" %in% operation ) {
        if (hurdle_direct) {
          if (pB$family=="poisson") {
            pois0 = ppois( pa_threshold, lambda = biom, lower.tail = FALSE)  #mu is expected 
            biom = biom * ( pa / pois0) 
          } else if ( pB$family=="nbinomial") {
            # need to check
          }
        } else  {
          # simple weighted sum -- this is identical to the hurdle -- but more general
          if (!is.null(pa_threshold)) pa = ifelse( pa< pa_threshold, 0, 1)
          biom = biom * pa
        }
        attr(biom, "unit") = "kg/km^2"
        return(biom)
      }
      if ("number" %in% operation) {
        if (hurdle_direct) {
          if (pN$family=="poisson") {
            pois0 = ppois( pa_threshold, lambda = nums, lower.tail = FALSE)  #mu is expected 
            nums = nums * (pa / pois0)
          } else if ( pB$family=="nbinomial") {
            # need to check
          }
        } else {
          #  simple weighted sum-- this is identical to the hurdle 
          if (!is.null(pa_threshold)) pa = ifelse( pa< pa_threshold, 0, 1)
          nums = nums * pa
        }
        attr(nums, "unit") = "n/km^2"
        return(nums)
      }
    }
  }

  if (length(operation) == 3 ) {
    if ("meansize" %in% operation ) {
      if ( "biomass" %in% operation ) {
        if ( "presence_absence" %in% operation ) {
          if (hurdle_direct) {
            if (pB$family=="poisson") {
              nums = biom / wgts
              pois0 = ppois( pa_threshold, lambda = nums, lower.tail = FALSE)  #mu is expected 
              nums = nums * (pa / pois0)  
            } else if ( pB$family=="nbinomial") {
              # need to check
            }
          } else {
          #  simple weighted sum-- this is identical to the hurdle 
            if (!is.null(pa_threshold)) pa = ifelse( pa< pa_threshold, 0, 1)
            nums = biom*pa / wgts
          }
          return(nums)
        }
      }
      if ( "number" %in% operation ) {
        if ( "presence_absence" %in% operation ) {
          if (hurdle_direct) {
            if (pN$family=="poisson") {
              biom = nums * wgts
              pois0 = ppois( pa_threshold, lambda=biom, lower.tail = FALSE)  #mu is expected 
              biom = biom * (pa / pois0)
            } else if ( pB$family=="nbinomial") {
              # need to check
            }
          } else {
          #  simple weighted sum-- this is identical to the hurdle 
            if (!is.null(pa_threshold)) pa = ifelse( pa< pa_threshold, 0, 1)
            biom = nums*pa * wgts
          }
          return(biom)
        }
      }
    }
  }

  return(operation)


}
