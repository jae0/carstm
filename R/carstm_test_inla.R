carstm_test_inla = function( family = "poisson" ) {

  library(INLA)
  data(Seeds)
  
  Seeds$rate = Seeds$r/ Seeds$n

  # NOTE CARSTM log-transforms offset so do not transform itsss

  if (family=="poisson") {

    # compare with GLM
    fitc = glm ( r ~ x1 + offset(log(n)), data=Seeds, family="poisson" )

    fit1 = glm ( r ~ x1, offset = log(n), data=Seeds, family="poisson" )

#### form used in carstm  "classical"
    fit2 = inla ( r ~ x1 + offset(log(n)), data=Seeds, family="poisson" ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1)
    )
####

    fit3 = inla ( r ~ x1, offset=log(n), data=Seeds,  family="poisson" ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1)
    )

    fit4 = inla ( r ~ x1, E=n, data=Seeds,  family="poisson",
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1)
    )


#### form used in carstm "compact"
    fit2e = inla ( r ~ x1 + offset(log(n)), data=Seeds, family="poisson"  ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1),
      inla.mode="compact"
    )

####

    fit3e = inla ( r ~ x1, offset=log(n), data=Seeds,  family="poisson" ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1),
      inla.mode="compact"
    )

    fit4e = inla ( r ~ x1, E=n, data=Seeds,  family="poisson",
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1),
      inla.mode="compact"
    )



    if (0) {

      # parameter estimates are stable ... the issue is that predictions are on different scales depeng upon run type

      summary(fitc)

      # Coefficients:
      #             Estimate Std. Error   z value Pr(>|z|)
      # (Intercept) -0.636577   0.057735 -11.02584  < 2e-16
      # x1          -0.119091   0.106761  -1.11549  0.26464

      summary(fit1)

      # Coefficients:
      #             Estimate Std. Error   z value Pr(>|z|)
      # (Intercept) -0.636577   0.057735 -11.02584  < 2e-16
      # x1          -0.119091   0.106761  -1.11549  0.26464

      summary(fit2)

      # Fixed effects:
      #               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
      # (Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
      # x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0

      summary(fit3)

      # Fixed effects:
      #               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
      # (Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
      # x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0

      summary(fit4)

      # Fixed effects:
      #               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
      # (Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
      # x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0
    
      summary(fit2e)

      # Fixed effects:
      #               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
      # (Intercept) -0.638 0.058     -0.752   -0.638     -0.525 -0.638   0
      # x1          -0.121 0.107     -0.331   -0.121      0.088 -0.121   0
      # Marginal log-Likelihood:  -77.93 

      summary(fit3e)

      # Fixed effects:
      #               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
      # (Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
      # x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0

      # Marginal log-Likelihood:  -77.93 


      summary(fit4e)

      # Fixed effects:
      #               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
      # (Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
      # x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0

      # Marginal log-Likelihood:  -77.93 

    }


    # predictions from marginals
    # both experimental and classical marginals.fitted.values are on response scale
    pfit2  = unlist( sapply( fit2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )  
    pfit2e = unlist( sapply( fit2e$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] ) 

    pfit3 = unlist( sapply( fit3$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
    pfit3e = unlist( sapply( fit3e$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    pfit4 = unlist( sapply( fit4$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
    pfit4e = unlist( sapply( fit4e$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    # changed behaviour 2023:
    pfit4 = pfit4 * Seeds$n  # ***
    pfit4e = pfit4e * Seeds$n  # ***


    # predictions from posterior distributions
    posterior_means = function( x, n=1000 ) {
      pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
      return(rowMeans( exp( pp ) ))  # <<< note conversion here
    }

    posterior_means2 = function( x, n=1000 ) {
      # method used in carstm  .. posterior is on link scale ((for both experimental and classical))
      pp = inla.posterior.sample( n, x, selection=list(Predictor=0)  )
      oo = apply( simplify2array( lapply( pp, function(x) { exp(x$latent) } ) ), 1, mean)
      return( oo )
    }

    
    # all on response scale:

    p0 = predict( fitc, type="response" )
    p1 = predict( fit1, type="response" )
    p2c = posterior_means2(fit2) 
    p2 = posterior_means(fit2) 

    # method used in carstm  .. 
    # posterior sims is on link scale ((for both experimental and classical))
    #  marginals.fitted.values are on response scale (for both experimental and classical)
    # summary.fitted.values on response scale ((for both experimental and classical))

    p2e = posterior_means(fit2e)  # captures offset information
    p2e2 = posterior_means2(fit2e) # captures offset information  -- method used in carstm
    p3 = posterior_means(fit3) 
    p3e = posterior_means(fit3e) 
    p4 = posterior_means(fit4) 
    p4e = posterior_means(fit4e) 

    # changed behaviour 2023:
    p4 = p4 * Seeds$n  # ***
    p4e = p4e * Seeds$n #***


    # summary.fitted.values on response scale ((for both experimental and classical))

    # all predictions
    u = cbind(
    
      p0,
      p1,
      # "classic"  .. .. incorporates offset in all computations
      fit2$summary.fitted.values$mean,  # with offset
      pfit2, # with offset
      p2,    # with offset

#### "experimental" ... used by default in carstm
#      fit2e$summary.fitted.values$mean * Seeds$n, # no offset
      fit2e$summary.fitted.values$mean  , # changed in 2023! ***
#      pfit2e * Seeds$n ,#  no offset needed any longer
      pfit2e  ,#  changed in 2023! ***
      p2e,  # with offset
      p2e2, # with offset
####

      fit3$summary.fitted.values$mean,
      pfit3,
      p3,

#      fit3e$summary.fitted.values$mean * Seeds$n, # experimental mode requires offset
      fit3e$summary.fitted.values$mean  , # experimental mode no longer requires offset *** 2023
#       pfit3e * Seeds$n ,
      pfit3e ,  # experimental mode no longer requires offset *** 2023 
      p3e,

      fit4$summary.fitted.values$mean  * Seeds$n,  # using E=... notation requires offset
#      pfit4 * Seeds$n,
      pfit4 ,  # no longer needs offsets  *** 2023
#      p4 * Seeds$n,
      p4 ,  # no longer needs offsets  *** 2023

      fit4e$summary.fitted.values$mean * Seeds$n,   # using E=... notation requires offset
#      pfit4e * Seeds$n,
      pfit4e ,  # no longer needs offsets  *** 2023
 #      p4e * Seeds$n   # no longer needs offsets  *** 2023
      p4e     # no longer needs offsets  *** 2023
      )

    plot( as.data.frame(u) )
    o = cor(u)  # should all be 1
    o [ abs(o-1) <= 0.001 ] = 1
    o [ abs(o-1) > 0.001 ] = 99999
    print(o)
    print( round( (u - rowMeans(u))/ rowMeans(u) , 4 ) )  
    if (any(o != 1)) {
      warning("\nThere has been a change in INLA's internals in terms of predcition scale!\n")
    } else {
      message("\nAll predictions look good!\n" )
    } 
  }


  if (family=="gaussian") {
    # simply convert to density/rates 


    # compare with GLM
    fit1 = glm ( r ~ n + as.factor(x1), data=Seeds, family="gaussian" )

    fit2 = inla ( r ~ n + as.factor(x1), data=Seeds, family="gaussian" ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1)
    )

####
    fit2e = inla ( r ~ n + as.factor(x1), data=Seeds, family="gaussian" ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1),
      inla.mode="compact"
    )
#####

    if (0) {

      # parameter estimates are stable ... the issue is that predictions are on different scales depeng upon run type

      summary(fit1)
 
        # Deviance Residuals: 
        #       Min          1Q      Median          3Q         Max  
        # -20.682338   -3.283449    0.594381    3.431895   16.331710  

        # Coefficients:
        #                 Estimate Std. Error  t value   Pr(>|t|)
        # (Intercept)    -1.444091   5.183555 -0.27859    0.78373
        # n               0.557116   0.088029  6.32878 5.7884e-06
        # as.factor(x1)1 -0.863782   4.253060 -0.20310    0.84134

        # (Dispersion parameter for gaussian family taken to be 69.0842404)

        #     Null deviance: 5169.238  on 20  degrees of freedom
        # Residual deviance: 1243.516  on 18  degrees of freedom
        # AIC: 153.3001



      summary(fit2)
        # Fixed effects:
        #                 mean    sd 0.025quant 0.5quant 0.975quant   mode kld
        # (Intercept)    -1.458 5.174    -11.713   -1.457      8.776 -1.456   0
        # n               0.557 0.088      0.383    0.557      0.732  0.557   0
        # as.factor(x1)1 -0.848 4.226     -9.220   -0.849      7.513 -0.850   0

        # Model hyperparameters:
        #                                         mean    sd 0.025quant 0.5quant 0.975quant  mode
        # Precision for the Gaussian observations 0.016 0.005      0.008    0.016      0.027 0.014

        # Marginal log-Likelihood:  -93.63 

      summary(fit2e)
        # Fixed effects:
        #                 mean    sd 0.025quant 0.5quant 0.975quant   mode kld
        # (Intercept)    -1.457 4.998    -11.337   -1.457      8.406 -1.456   0
        # n               0.557 0.085      0.389    0.557      0.725  0.557   0
        # as.factor(x1)1 -0.849 4.084     -8.920   -0.850      7.211 -0.851   0

        # Model hyperparameters:
        #                                         mean    sd 0.025quant 0.5quant 0.975quant  mode
        # Precision for the Gaussian observations 0.016 0.005      0.008    0.015      0.028 0.014

        # Marginal log-Likelihood:  -93.62 

    }


    # predictions from marginals

    pfit2  = unlist( sapply( fit2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    pfit2e  = unlist( sapply( fit2e$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
   
    # predictions from posterior distributions  -- identity link
    posterior_means = function( x, n=100 ) {
      pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
      return(rowMeans( ( pp ) ))
    }
  
    p1 = predict( fit1, type="response" )
    p2 = posterior_means(fit2) 
    p2e = posterior_means(fit2e) 
  
       # all predictions
    u = cbind(
      p1,
      
      fit2$summary.fitted.values$mean,
      pfit2,
      p2,

      fit2e$summary.fitted.values$mean,
      pfit2e,
      p2e

     )

    plot( as.data.frame(u) )
    o = cor(u)  # should all be 1
    o [ abs(o-1) <= 0.001 ] = 1
    o [ abs(o-1) > 0.001 ] = 99999
    print(o)
    print( round( (u - rowMeans(u))/ rowMeans(u) , 4 ) )  
    if (any(o != 1)) {
      warning("\nThere has been a change in INLA's internals in terms of predcition scale!\n")
    } else {
      message("\nAll predictions look good!\n" )
    } 
  }


  if (family=="binomial") {

    Seeds$Y = ifelse( (Seeds$r / Seeds$n) > 0.5, 1, 0 )

    # compare with GLM
    fit1 = glm ( Y ~ as.factor(x1) , data=Seeds, family="binomial" )

    fit2 = inla ( Y ~ as.factor(x1), data=Seeds, family="binomial" ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1)
    )

####
    fit2e = inla ( Y ~ as.factor(x1), data=Seeds, family="binomial"  ,
      control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
      control.predictor = list(compute=TRUE, link=1),
      inla.mode="compact"
    )
####

    if (0) {

      # parameter estimates are stable ... the issue is that predictions are on different scales depeng upon run type
 
      summary(fit1)
        # Coefficients:
        #                 Estimate Std. Error  z value Pr(>|z|)
        # (Intercept)     0.559616   0.626783  0.89284  0.37194
        # as.factor(x1)1 -1.406914   0.932227 -1.50920  0.13125

        # (Dispersion parameter for binomial family taken to be 1)

        #     Null deviance: 29.06454  on 20  degrees of freedom
        # Residual deviance: 26.63789  on 19  degrees of freedom
        # AIC: 30.63789



      summary(fit2)
        # Fixed effects:
        #                 mean    sd 0.025quant 0.5quant 0.975quant   mode kld
        # (Intercept)     0.560 0.626     -0.614     0.54      1.848  0.499   0
        # as.factor(x1)1 -1.479 0.931     -3.389    -1.45      0.267 -1.393   0

        # Marginal log-Likelihood:  -16.69 


 
      summary(fit2e)
        # Fixed effects:
        #                 mean    sd 0.025quant 0.5quant 0.975quant   mode kld
        # (Intercept)     0.608 0.627     -0.622    0.608      1.838  0.608   0
        # as.factor(x1)1 -1.542 0.932     -3.372   -1.542      0.285 -1.542   0

        # Marginal log-Likelihood:  -16.69 

    }


    # predictions from marginals

    inverse.logit = function( x ) {
      # x should be the log odds ratio
      oddsratio = exp(x)
      prob = oddsratio / (1 + oddsratio )
      return (prob)
    }


    pfit2  = unlist( sapply( fit2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
    pfit2e = unlist( sapply( fit2e$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    # predictions from posterior distributions -- logit link 
    posterior_means = function( x, n=100 ) {
      pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
      return(rowMeans( inverse.logit ( pp ) ))
    }

    p1 = predict( fit1, type="response" )
    p2 = posterior_means(fit2) 
    p2e = posterior_means(fit2e) 

    # all predictions .. on probability scale ... no need for scaling
    u = cbind(
   
      p1,
      
      fit2$summary.fitted.values$mean,
      pfit2,
      p2,

      fit2e$summary.fitted.values$mean , 
      pfit2e ,
      p2e
    )

    plot( as.data.frame(u) )
    o = cor(u)  # should all be 1
    o [ abs(o-1) <= 0.001 ] = 1
    o [ abs(o-1) > 0.001 ] = 99999
    print(o)
    print( round( (u - rowMeans(u))/ rowMeans(u) , 4 ) )  
    if (any(o != 1)) {
      warning("\nThere was been a change in INLA's internals in terms of predcition scale!\n")
    } else {
      message("\nAll predictions look good!\n" )
    } 
  }


}