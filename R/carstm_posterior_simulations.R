
carstm_posterior_simulations = function( p=NULL, pN=NULL, pW=NULL, pH=NULL, pB=NULL, sppoly=NULL, 
  wgts_max=NULL, N_max=NULL, B_max=NULL, max_value=NULL, degree_day=FALSE  ) {

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
    gen = carstm_model( p=pW, DS="carstm_modelled_summary", sppoly=sppoly  )
    gen = gen[[ "predictions_posterior_simulations"  ]]
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


  # construct meansizes matrix used to convert number to weight
  if ( "presence_absence" %in% operation ) {
    pa = carstm_model( p=pH, DS="carstm_modelled_summary", sppoly=sppoly  )
    pa = pa[[ "predictions_posterior_simulations" ]]
    j = which( !is.finite(pa) )
    if (length(j) > 0 ) pa[j] = NA
    # pa = inverse.logit(pa)
    attr( pa, "unit") = "probability"
    if (length(operation)==1) return(pa)
  }


  if ( "meansize" %in% operation) {
    wgts = carstm_model( p=pW, DS="carstm_modelled_summary", sppoly=sppoly  )
    wgts = wgts[[ "predictions_posterior_simulations"  ]]
    j = which( !is.finite(wgts) )
    if (length(j) > 0 ) wgts[j] = NA

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
    biom = carstm_model( p=pB, DS="carstm_modelled_summary", sppoly=sppoly  )
    biom = biom[[ "predictions_posterior_simulations" ]] 
    j = which( !is.finite(biom) )
    if (length(j) > 0 ) biom[j] = NA
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
    nums = carstm_model( p=pN, DS="carstm_modelled_summary", sppoly=sppoly  )
    nums = nums[[ "predictions_posterior_simulations" ]]    
    j = which( !is.finite(nums) )
    if (length(j) > 0 ) nums[j] = NA
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
    
    if ("pa" %in% operation) {
      if ( "biomass" %in% operation ) {
        biom = biom * pa
        attr(biom, "unit") = "kg/km^2"
        return(biom)
      }
      if ("number" %in% operation) {
        nums = nums * pa
        attr(nums, "unit") = "n/km^2"
        return(nums)
      }
    }
  }

  if (length(operation) == 3 ) {
    if ("meansize" %in% operation ) {
      if ( "biomass" %in% operation ) {
        if ( "presence_absence" %in% operation ) {
          nums = biom*pa / wgts
          return(nums)
        }
      }
      if ( "number" %in% operation ) {
        if ( "presence_absence" %in% operation ) {
          biom = nums*pa*wgts
          return(biom)
        }
      }
    }
  }

  return(operation)


}
