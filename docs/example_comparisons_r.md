
# Example comparing a simple spatial only model (BYM) formulation with INLA, and checking they are the same as with CARSTM 

```r
  # Prepare data:
  require(INLA)

  # inla.pardiso.check()  # get it if you can and turn it on with:
  # inla.setOption(pardiso.license="~/paradiso.license" )  # point "pardiso.license" to the physical location of your license


  data(Germany)
  g = system.file("demodata/germany.graph", package="INLA")
  source(system.file("demodata/Bym-map.R", package="INLA"))
  summary(Germany)


  # Basic BYM model using INLA directly:
  Germany$region.iid = Germany$region 
  fm = Y ~ f(region, model="besag",graph.file=g) + f(region.iid, model="iid") + f(x, model="rw2")
  fit =  inla( fm, family="poisson", data=Germany, E=E, verbose=TRUE,
      control.predictor = list( compute=TRUE),
    control.compute = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
    inla.mode="classic"
 )
  summary(fit)


  # Same analysis using INLA's experiemental mode (faster and more efficient RAM usage) with marginals, and bym2 model
  fm2 = Y ~ f( region, model="bym2",graph.file=g) + f(x, model="rw2")
  fit2 =  inla( fm2, family="poisson", data=Germany, E=E, verbose=TRUE,
    # control.inla = list( strategy='adaptive' ) , 
    # control.inla = list( strategy='adaptive', int.strategy='eb' ),
    control.predictor = list( compute=TRUE  ),
    control.compute = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
    inla.mode="compact"
  )
  summary(fit2)
```

Checking some plots

```r
  
  plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=TRUE, single=TRUE )
   
  # random effects are nearly identical 
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
```


```r 
  # Same analysis, now using carstm 
  # there are extra steps as a data structure and options need to be specified 
  # and posterior sims to compute combined spatial (and spatiotemporal) effects "re_total"
  
  require(carstm)  # carstm options are sent via one controlling parameter list  
  # loadfunctions("carstm")

  Germany$tag = "predictions"  # predict on all locations
  
  Germany$region = as.character(Germany$region)

  sppoly = Germany  # construct "sppoly" with required attributes (though it is not a polygon)
  attributes(sppoly)[["areal_units_fn"]] = g  
  attributes(sppoly)$nb = inla.read.graph(g)

  p = list(
    modeldir = tempdir(),
    carstm_model_label = "testlabel",
    carstm_modelengine = "inla",
    vn = list(Y="Y", S="region", SI="region.iid", O = "E"),  # instructions on which are spatial and iid and offset terms. parsing will try to figure it out but this is most secure
    dimensionality = "space",   # a pure space model
    family = "poisson",
    nposteriors=5000
  ) 

  p$formula = formula( Y ~  1 + offset( E ) 
    + f(region, model="besag", graph=slot(sppoly, "nb"), scale.model=TRUE ) 
    + f(region.iid, model="iid" ) 
    + f(x, model="rw2", scale.model=TRUE)  )

  
  p$formula = formula( Y ~  1 + offset( E ) 
    + f(region, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE ) 
    + f(region.iid, model="iid") 
    + f(x, model="rw2", scale.model=TRUE)  )

  # Leroux model
  p$formula = formula( Y ~  1 + offset( E ) 
    + f(region, model="besagproper2", graph=slot(sppoly, "nb"), scale.model=TRUE ) 
    + f(region.iid, model="iid" ) 
    + f(x, model="rw2", scale.model=TRUE)   )



  # note offset is not logged ... link function handles it 

  fn_fit = tempfile( pattern="fit", tmpdir=p$modeldir )

  res = carstm_model( 
    p=p,
    data = Germany, 
    sppoly = sppoly,
    space.id = as.character(Germany$region),
    fn_fit = fn_fit, 
    num.threads="4:2",
    verbose=TRUE
  ) 

  if (0) {
    # Or, to load currently saved results
    res = carstm_model( p=p, DS="carstm_summary"  ) 
    
    # extract currently saved model fit
    fit3 = carstm_model( p=p, fn_fit = fn_fit, DS="carstm_modelled_fit" )  
    
    fit3$summary$dic$dic
    fit3$summary$dic$p.eff
    summary(fit3)  # identical to fit2

    plot(fit3)
    plot(fit3, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  }

  plot(fit$summary.random$region$mean ~ fit3$summary.random$region$mean[1:544] )
  
  plot(res$predictions[,"mean"] ~ fit$summary.fitted.values[,"mean"] )

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


  # carstm_maps of some of the results  .. can't map unitl polygons are coverted to sf format (TODO)
  plt = carstm_map(  res=res, vn=c( "predictions" ), 
    space="region",
    plot_elements=c( "compass", "scale_bar", "legend" ),
    main=paste( "Predictions")  
  )

  plt = carstm_map(  res=res, vn=c( "random", "region", "re_total" ), 
    space="region",
    plot_elements=c( "compass", "scale_bar", "legend" ),
    main=paste( "Spatial effects")  
  )



  if (0) {

    # this section compares: mean(exp(x)) vs exp(mean(x))
    # in this example, the differences are small .. but can be large deponding upon data distribtion
    
      data(warpbreaks)
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
      list_simplify = function(x) as.data.frame( t( as.data.frame( x )))

      t( lapply( inla.zmarginal( m), FUN=invlink )) # NOTE, sd is incorrect .. just looking at means 

      # this method recovers the the correct SD
      t( inla.zmarginal( inla.tmarginal( invlink, m) , silent=TRUE  )  )

      # or, step-by-step
      g = inla.tmarginal( invlink, m)    
      
      zmarg = function(yy) inla.zmarginal( yy, silent=TRUE  )

      invlink(unlist( list_simplify ( sapply( list(m), zmarg ) ) ))
      
      list_simplify ( sapply( list(g), zmarg ) )


  }

```

