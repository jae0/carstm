carstm_test_inla = function( family = "poisson" ) {

  library(INLA)
  data(Seeds)
  
  Seeds$rate = Seeds$r/ Seeds$n

  # NOTE CARSTM log-transforms offset so do not transform it
     # predictions from posterior distributions
    posterior_means = function( x, n=1000, func=function(x){x} ) {
      pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
      return(rowMeans( func( pp ) ))  # <<< note conversion here
    }
 
    posterior_means2 = function( x, n=1000, func=function(x){x} ) {
      # method used in carstm  .. posterior is on link scale ((for both experimental and classical))
      pp = inla.posterior.sample( n, x, selection=list(Predictor=0)  )
      oo = apply( simplify2array( lapply( pp, function(x) { func(x$latent) } ) ), 1, mean)
      return( oo )
    }

    
  if (family=="besag_poisson_nbinomial") {
 
    source(system.file("demodata/Bym-map.R", package="INLA"))  # load polygons "germany.graph", etc
    g = system.file("demodata/germany.graph", package = "INLA")
    
    Germany = INLA::Germany  # data
    Germany$region.copy = Germany$region
    Germany$logE = log(Germany$E)

    obsrate = Germany$Y / Germany$E

    
    # NOTE: all models give same parameter estimates

    # besag model .. offset not in formula but as a parameter for inla
    formula1 = Y ~ f(region.copy, model="besag", graph=g) + f(region,model="iid")
    

    # no offsets in formula, but use "offset log transformed" ( note: results are counts and not rates)
    m1a = inla(formula1, family="poisson", data=Germany, 
        offset=Germany$logE,
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    )

    # no offsets in formula, but use "E" (eta) not log transformed
    m1 = inla(formula1, family="poisson", data=Germany, E=Germany$E, 
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    )

    # not same results: one is a rate the the other is a count
    plot(m1a$summary.fitted.values$mean, m1$summary.fitted.values$mean)   
    PREDICTEDCOUNTS = m1$summary.fitted.values$mean * Germany$E  # rate * offset
    plot(m1a$summary.fitted.values$mean, PREDICTEDCOUNTS ) 


    # NOTE: using offset in formula is the approach used in carstm: it is consistent 
    # besag model .. offset in formula
    formula2 = Y ~ offset(logE) + f(region.copy, model="besag",graph=g) + f(region,model="iid")

    # offset in formula
    m2 = inla(formula2, family="poisson", data=Germany, 
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    )


    p1 = posterior_means(m1, func=exp)
    p1marginals  = unlist( sapply( m1$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    # results are rates and not counts: also, varying offset and scale issues no longer exist:
    if (0) {
      plot( obsrate ~ m1$summary.fitted.values$mean ) 
      plot( obsrate ~ p1 )   
      plot( obsrate ~ p1marginals )
    }  

    # however, using the "offset" argument returns counts
    p1a = posterior_means(m1a, func=exp)
    p1amarginals  = unlist( sapply( m1a$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    # results are rates and not counts: also, varying offset and scale issues no longer exist:
    if (0) {
      plot( Germany$Y ~ m1a$summary.fitted.values$mean ) # fitted.values on user scale and ignores offsets for prediction p.. prediction is rate
      plot( Germany$Y ~ p1a )  # posterior samples on link scale that are inverted (log/exp) and ignores offsets for prediction .. prediction is rate
      plot( Germany$Y ~ p1amarginals ) # posterior marginals also on user scale and ignores offsets for prediction  .. prediction is rate
    }    

 
    # results are counts and not rates:
    # NOTE: the varying offset and scale issues now fixed ... all are same
    p2 = posterior_means(m2, func=exp)  # range of p2 larger
    p2marginals  = unlist( sapply( m2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    if (0) {
      plot( Germany$Y ~ m2$summary.fitted.values$mean  ) #  summary.fitted.values$mean is on user scale with offset .. prediction is count
      plot( Germany$Y ~ p2 ) # posterior samples are on link scale .. which are then exponentiated and with offsets .. prediction is count
      plot( Germany$Y ~ p2marginals ) # marginals are again on user scale with offsets .. prediction is count
    }



    # negative binomial .. has offsets in formula
    nb2 = inla(formula2, family="nbinomial", data=Germany, 
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    )

    pnb2 = posterior_means(nb2, func=exp)
    pnb2marginals  = unlist( sapply( nb2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    if (0) {
      plot( Germany$Y ~ nb2$summary.fitted.values$mean ) # fitted.values on user scale and ignores offsets for prediction .. prediction is rate even when offsets given
      plot( Germany$Y ~ pnb2 )  # posterior means incorporates given offsets and predicts a count 
      plot( Germany$Y ~ pnb2marginals ) #   posterior marginals also on user scale and ignores offsets for prediction  .. prediction is rate
    } 
   
    # bym2 model .. offset in formula  
    formula3 = Y ~ offset(logE) + f(region, model="bym2",graph=g)  

    m3 = inla(formula3, family="poisson", data=Germany, 
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    )
   
    p3 = posterior_means(m3, func=exp)
    p3marginals  = unlist( sapply( m3$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

    # note: results are in counts
    if (0) {
      plot( Germany$Y ~ m3$summary.fitted.values$mean ) # fitted.values on user scale and ignores offsets for prediction p.. prediction is rate
      plot( Germany$Y ~ p3 )  # posterior samples on link scale that are inverted (log/exp) and ignores offsets for prediction .. prediction is rate
      plot( Germany$Y ~ p3marginals ) # posterior marginals also on user scale and ignores offsets for prediction  .. prediction is rate

      # recall marginals.fitted.values and summary.fitted.values are on user scale (counts) 
      plot( m3$marginals.fitted.values[[1]] )
      abline(v=m3$summary.fitted.values[1,"mean"])
  
      # however, marginals.random$region and summary.random$region are on internal link scale ... need to be inverted
      x = m3$marginals.random$region[[1]]
      plot( x[,1], x[,2] )
      abline(v=m3$summary.random$region[1,"mean"]) # (e.g., counts cannot be negative values)
      
      plot( exp(x[,1]), x[,2] )
      abline(v=exp(m3$summary.random$region[1,"mean"]))

      # note: posterior samples would also need an inverse transform

    }    
    
      
    u =  cbind(Germany$Y, 
      m1$summary.fitted.values$mean * Germany$E, p1 * Germany$E, p1marginals * Germany$E, 
      m1a$summary.fitted.values$mean, p1a, p1amarginals, 
      m2$summary.fitted.values$mean, p2, p2marginals, 
      m3$summary.fitted.values$mean, p3, p3marginals, 
      nb2$summary.fitted.values$mean, pnb2, pnb2marginals ) 
    plot( as.data.frame(u) )
    o = cor(u)  # should all be 1
    o [ abs(o-1) <= 0.05 ] = 1
    o [ abs(o-1) > 0.05 ] = 0
    print(o)
    print( round( (u - rowMeans(u))/ rowMeans(u) , 4 ) )  
    if (any(o != 1)) {
      warning("\nThere has been a change in INLA's internals in terms of predcition scale!\n")
    } else {
      message("\nAll predictions for poisson/nbinomial besag + bym2 look good!\n" )
    } 

    # note: lognormal is handled as log Y ~ normal() ... be careful .. manual back transforms are required

  }



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



    # all on response scale:

    p0 = predict( fitc, type="response" )
    p1 = predict( fit1, type="response" )
    p2c = posterior_means2(fit2, func=exp) 
    p2 = posterior_means(fit2, func=exp ) 

    # method used in carstm  .. 
    # posterior sims is on link scale ((for both experimental and classical))
    #  marginals.fitted.values are on response scale (for both experimental and classical)
    # summary.fitted.values on response scale ((for both experimental and classical))

    p2e = posterior_means(fit2e, func=exp)  # captures offset information
    p2e2 = posterior_means2(fit2e, func=exp) # captures offset information  -- method used in carstm
    p3 = posterior_means(fit3, func=exp) 
    p3e = posterior_means(fit3e, func=exp) 
    p4 = posterior_means(fit4, func=exp) 
    p4e = posterior_means(fit4e, func=exp) 

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
    o [ abs(o-1) > 0.001 ] = 0
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
    
  
    p1 = predict( fit1, type="response" )
    p2 = posterior_means(fit2)  # identity link
    p2e = posterior_means(fit2e) # identity link
   
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
    o [ abs(o-1) > 0.001 ] = 0
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
 
    p1 = predict( fit1, type="response" )
    p2 = posterior_means(fit2, func=inverse.logit) 
    p2e = posterior_means(fit2e, func=inverse.logit) 

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
    o [ abs(o-1) > 0.001 ] = 0
    print(o)
    print( round( (u - rowMeans(u))/ rowMeans(u) , 4 ) )  
    if (any(o != 1)) {
      warning("\nThere was been a change in INLA's internals in terms of predcition scale!\n")
    } else {
      message("\nAll predictions look good!\n" )
    } 
  }


  if (family=="bym2_poisson_nbinomial_random_effect") {
 
    source(system.file("demodata/Bym-map.R", package="INLA"))  # load polygons "germany.graph", etc
    g = system.file("demodata/germany.graph", package = "INLA")
    
    Germany = INLA::Germany  # data
    Germany$region.copy = Germany$region
    Germany$logE = log(Germany$E)

    obsrate = Germany$Y / Germany$E

     
    # bym2 model .. offset not in formula but as a parameter for inla
    formula1 = Y ~ f(region, model="bym2", graph=g)  
    

    # no offsets in formula, but use "offset log transformed" ( note: results are counts and not rates)
    m1a = inla(formula1, family="poisson", data=Germany, 
        offset=Germany$logE,
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    ) 
    
    m1a$summary.fitted.values$mean
    

    # NOTE: using offset in formula is the approach used in carstm:
    formula2 = Y ~ offset(logE) + f(region, model="bym2",graph=g)  

    m2 = inla(formula2, family="poisson", data=Germany, 
        control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
        control.fixed=list(prec.intercept=1),
        control.predictor = list(compute=TRUE, link=1)
    )

    m2$summary.fitted.values$mean

    cor(m1a$summary.fitted.values$mean, m2$summary.fitted.values$mean)
    plot(m1a$summary.fitted.values$mean, m2$summary.fitted.values$mean)

    # parameter estimates are good

    # random effects:
    re_m1a_summ = m1a$summary.random[["region"]]
    re_m1a_marg = m1a$marginals.random[["region"]]

    re_m2_summ = m2$summary.random[["region"]]
    re_m2_marg = m2$marginals.random[["region"]]

    abs(re_m2_summ[1,"mean"] - re_m1a_summ[1,"mean"]) < 1e6
    
    abs(inla.zmarginal(re_m2_marg[[1]])[[1]] - inla.zmarginal(re_m1a_marg[[1]])[[1]] ) < 1e6
 
    abs( inla.tmarginal(function(x) (x), re_m2_marg[[1]]) [[1]] - re_m2_summ[1,"mean"] ) < 1e6

    o1a = lapply( re_m1a_marg, inla.tmarginal, fun=exp )

    o2  = lapply( re_m2_marg,  inla.tmarginal, fun=exp )

    plot(o2[1,], o1a[1,])
    cor(c(o2), c(o1a))
    
    o2t = sapply( lapply( o2, inla.zmarginal), function(x) x$mean )
    o2z  = sapply( lapply( re_m2_marg,  inla.zmarginal ), function(x) x$mean )
    
    plot(o2z~o2t )
 
    plot(o2z~log(o2t) )
 
    plot(o2z~re_m2_summ$mean)
    cor(o2z,re_m2_summ$mean)

  }



}





### NOTE:: Bottom line 
###     when using offsets in formulae, marginals and summaries are on user scale.
###     using posterior samples requires inverse link    
### 