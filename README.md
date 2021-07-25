# carstm
carstm provides a sequence of examples with supporting functions that examines abundance estimation of snow crab and groundfish (Atlantic cod) in
Maritimes groundfish survey strata. The approach is easily generalizable to all species. The main sequence of example models and scripts are found in inst\scripts\* for the projects: bio.snowcrab and aegis.surveys, and leverages the https::\github.com\jae0\aegis data and GIS handing routines.

For Atlantoc cod, carstm replicates the standard analysis which is known as "stratanal", a basic stratified average estimate. This is shown to be equivalent to a Gaussian linear fixed effects model. Thereafter, a model-based approach is used to incrementally improve upon the assumptions of the model, focussing upon the distributional model (Poisson, overdispersed Poisson), adding environmental covariates and then employing an INLA-based ICAR (intrinsic conditionally autoregressive models; "bym2") approach towards accounting for areal unit modelling and an AR1 temporal autocorrelation assuming separability of the spacetime autocorrelation.


## Fast example of a pure space model:
  # Prepare data:

  require(INLA)
  data(Germany)
  g = system.file("demodata/germany.graph", package="INLA")
  source(system.file("demodata/Bym-map.R", package="INLA"))
  summary(Germany)


  # Basic BYM model using INLA directly:

  fm = Y ~ f(region, model="besag",graph.file=g) + f(region, model="iid") + f(x, model="rw2")
  fit =  inla( fm, family="poisson", data=Germany, E=E, verbose=TRUE,
      control.predictor = list( compute=TRUE),
    control.compute = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
 )


  # Same example using INLA's experiemental mode (faster and more efficient RAM usage) with marginals, and bym2 model

  fm2 = Y ~ f( region, model="bym2",graph.file=g) + f(x, model="rw2")
  fit2 =  inla( fm2, family="poisson", data=Germany, E=E, verbose=TRUE,
    # control.inla = list( strategy='adaptive' ) , 
    # control.inla = list( strategy='adaptive', int.strategy='eb' ),
    control.predictor = list( compute=TRUE  ),
    control.compute = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
    inla.mode="experimental"
  )

  summary(fit)
  summary(fit2)

  # random effects
  dev.new()
  par( mfrow=c(1,2))
  Bym.map(fit$summary.random$region$mean)
  Bym.map(fit2$summary.random$region$mean[1:nrow(Germany)])  # bym2 concatenates iid and bym effects 

  dev.new()
  plot(fit$summary.random$region$mean ~ fit2$summary.random$region$mean[1:544] )

  dev.new()
  hist(fit$summary.random$region$mean)

  # predictions
  dev.new()
  par( mfrow=c(1,2))
  Bym.map(fit$summary.fitted.values$mean)
  Bym.map(exp( fit2$summary.fitted.values$mean) )  # bym2 concatenates iid and bym effects 

  dev.new()
  plot(fit$summary.fitted.values$mean ~ exp( fit2$summary.fitted.values$mean)  )

  dev.new()
  hist(fit$summary.fitted.values$mean)


  # Same example using carstm .. extra steps as a data structure and options need to be specified and posterior sims to compute cobined spatial and spatiotemporal effects
  require(carstm)  # carstm options are sent via one controlling parameter list  
  
  Germany$tag = "predictions"  # predict on all locations
  Germany$log_E = log( Germany$E )  # offset on log scale
  Germany$region = as.character(Germany$region)

  sppoly = Germany  # construct "sppoly" with required attributes
  attributes(sppoly)[["areal_units_fn"]] = g  
  attributes(sppoly)[["nb"]] = inla.read.graph(g) 

  p = list(
    project_name = "test",  # basic run
    modeldir = tempdir(),
    carstm_model_label = "testlabel",
    variabletomodel = "Y",
    vnS = "region", 
    aegis_dimensionality = "space",   # a pure space model
    carstm_modelengine = "inla",
    carstm_model_formula = Y ~ offset(log_E) +  f(region, model="bym2", graph.file=g)  + f(x, model="rw2"),
    carstm_model_family = "poisson",
    nposteriors=5000
  ) 

  fn_fit = tempfile( pattern="fit", tmpdir=p$modeldir )
  fn_res = tempfile( pattern="res", tmpdir=p$modeldir )

  res = carstm_model( 
    p=p,
    M = Germany, 
    sppoly = sppoly,
    region.id = as.character(Germany$region),
    fn_fit = fn_fit,
    fn_res = fn_res,
    verbose=TRUE
  ) 

   or # res = carstm_model( p=p, fn_res = fn_res, DS="carstm_modelled_summary"  ) # to load currently saved results

     # extract results
  fit3 = carstm_model( p=p, fn_fit = fn_fit, DS="carstm_modelled_fit" )  # extract currently saved model fit
  fit3$summary$dic$dic
  fit3$summary$dic$p.eff
  summary(fit3)  # identical to fit2

  plot(fit3)
  plot(fit3, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

  plot(fit$summary.random$region$mean ~ fit3$summary.random$region$mean[1:544] )




  # random effects
  dev.new()
  par( mfrow=c(1,3))
  Bym.map(fit$summary.random$region$mean)
  Bym.map(fit2$summary.random$region$mean[1:nrow(Germany)])  # bym2 concatenates iid and bym effects 
  Bym.map(fit3$summary.random$region$mean[1:nrow(Germany)])  # bym2 concatenates iid and bym effects 

  dev.new()
  par( mfrow=c(1,2))
  plot(fit$summary.random$region$mean ~ fit2$summary.random$region$mean[1:544] )
  plot(fit$summary.random$region$mean ~ fit3$summary.random$region$mean[1:544] )

  dev.new()
  par( mfrow=c(1,3))
  hist(fit$summary.random$region$mean)
  hist(fit2$summary.random$region$mean)
  hist(fit3$summary.random$region$mean)


  # predictions
  dev.new()
  par( mfrow=c(1,4))
  Bym.map(fit$summary.fitted.values$mean)
  Bym.map(exp( fit2$summary.fitted.values$mean) )  # bym2 concatenates iid and bym effects 
  Bym.map(exp( fit3$summary.fitted.values$mean) )  # bym2 concatenates iid and bym effects 
  Bym.map( res$predictions[,"mean"] )  # bym2 concatenates iid and bym effects 

  dev.new()
  par( mfrow=c(1,3))
  plot(fit$summary.fitted.values$mean ~ exp( fit2$summary.fitted.values$mean)  )
  plot(fit$summary.fitted.values$mean ~ exp( fit3$summary.fitted.values$mean)  )
  plot(exp(fit3$summary.fitted.values$mean) ~ res$predictions[,"mean"]  )

  dev.new()
  par( mfrow=c(1,4))
  hist(fit$summary.fitted.values$mean)
  hist(exp(fit2$summary.fitted.values$mean))
  hist(exp(fit3$summary.fitted.values$mean))
  hist(res$predictions[,"mean"])


  plot(exp(fit3$summary.fitted.values$mean) ~ res$predictions[,"mean"]  )

  w = which(res$predictions[,"mean"] > 5) [1] 

  res$predictions[w,]
  fit3$summary.fitted.values[w,]


  plot( fit$marginals.fitted.values[[w]] )
  plot( fit2$marginals.fitted.values[[w]] )
  plot( exp(fit3$marginals.fitted.values[[w]] ))


  invlink = function(x) inla.link.log( x,  inverse=TRUE )
  invlinktransf = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )
  list_simplify = function(x) as.data.frame( t( as.data.frame( x )))
  summary_inv_predictions = function(x) inla.zmarginal( x, silent=TRUE  )

  m = fit3$marginals.fitted.values[[w]]


      o0 = inla.zmarginal( m)
      o = lapply( o0, FUN=invlink )
      t(o)  # sd is incorrect .. 

      ot = invlinktransf( m)
      t(ot)

      exp(fit3$summary.fitted.values[w,])


      # step-by-step

      plot(exp(m))

      g = list( invlink( m ) )[[1]]
      plot(g)

      zmarg = function(x) inla.zmarginal( x, silent=TRUE  )
  
      h = sapply( g, zmarg )
      h = list_simplify ( h )
      h

      h2 = list_simplify ( sapply( g, summary_inv_predictions ) )
      h2


# ----
  m = fit3$marginals.fitted.values[[w]]
  ot = invlinktransf( m)
  t(ot)

  zmarg = function(x) inla.zmarginal( x, silent=TRUE  )
  h = list_simplify ( sapply( list( invlink( m ) ), zmarg ) )
  h

  h2 = list_simplify ( sapply( g, summary_inv_predictions ) )


  inla.mmarginal(g)
  exp( inla.mmarginal(m) )

  inla.zmarginal(g)
  inla.zmarginal(m) 


  if (0) {

      data(warpbreaks)


      warpbreaks$breaks [sample(1:nrow(warpbreaks), 20) ] = 0
   
      str(warpbreaks)
      
      fit0 = glm(breaks ~ wool + tension, warpbreaks, family = poisson(link = "log"))
      summary(fit0)

      fit1 = glm(breaks ~ wool + tension, warpbreaks, family = quasipoisson(link = "log"))
      summary(fit1)

      fit2 = inla(breaks ~ wool + tension, data=warpbreaks, family="poisson" )
      summary(fit2)

      fit3 = inla(breaks ~ wool + tension, data=warpbreaks, family="zeroinflatedpoisson0" )
      summary(fit3)

      ## chose a marginal and compare the with the results computed by the

      res = fit2

      r = res$summary.fixed["woolB",]
      r

      m = res$marginals.fixed$woolB


      ## compute the 95% HPD interval
      inla.hpdmarginal(0.95, m)
     
      invlink = function(x) inla.link.log( x,  inverse=TRUE )
      invlinktransf = function(x) inla.zmarginal( inla.tmarginal( invlink, x) , silent=TRUE  )
      list_simplify = function(x) as.data.frame( t( as.data.frame( x )))
      summary_inv_predictions = function(x) inla.zmarginal( x, silent=TRUE  )

      o0 = inla.zmarginal( m)
      o = lapply( o0, FUN=invlink )
      t(o)  # sd is incorrect .. 

      ot = invlinktransf( m)
      t(ot)


      # step-by-step
      
      # g = list( invlink( m ) )
      g = inla.tmarginal( invlink, m)    
      
      zmarg = function(yy) inla.zmarginal( yy, silent=TRUE  )
  
      list_simplify ( sapply( list(m), zmarg ) )

      list_simplify ( sapply( list(g), zmarg ) )


  }



    # carstm_maps of some of the results  .. can't map unitl polygons are coverted to sf format (TODO)
  tmout = carstm_map(  res=res, vn=c( "predictions" ), 
    plot_elements=c( "compass", "scale_bar", "legend" ),
    main=paste( "Predictions")  
  )

  tmout = carstm_map(  res=res, vn=c( "random", "spatial", "combined" ), 
    plot_elements=c( "compass", "scale_bar", "legend" ),
    main=paste( "Spatial effects")  
  )

