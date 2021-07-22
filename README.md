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
  Germany = cbind(Germany,region.struct=Germany$region)


  # Basic BYM model using INLA directly:

  fm = Y ~ f(region.struct, model="besag",graph.file=g) + f(region, model="iid") + f(x, model="rw2")
  fit =  inla( fm, family="poisson", data=Germany, E=E, verbose=TRUE,
      control.predictor = list( compute=TRUE),
    control.compute = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
 )


  # Same example using INLA's experiemental mode (faster and more efficient RAM usage) with marginals, and bym2 model

  fm2 = Y ~ f( region.struct, model="bym2",graph.file=g) + f(x, model="rw2")
  fit2 =  inla( fm2, family="poisson", data=Germany, E=E, verbose=TRUE,
    # control.inla = list( strategy='adaptive' ) , 
    # control.inla = list( strategy='adaptive', int.strategy='eb' ),
    control.predictor = list( compute=TRUE ),
    control.compute = list(dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
    inla.mode="experimental"
  )

  summary(fit)
  summary(fit2)

  dev.new()
  par( mfrow=c(1,2))
  Bym.map(fit$summary.random$region.struct$mean)
  Bym.map(fit2$summary.random$region.struct$mean[1:nrow(Germany)])  # bym2 concatenates iid and bym effects 

  plot(fit$summary.random$region.struct$mean ~ fit2$summary.random$region.struct$mean[1:544] )


  # Same example using carstm .. extra steps as a data structure and options need to be specified 
  require(carstm)  # carstm options are sent via one controlling parameter list  
  
  Germany$tag = "predictions"  # predict on all locations
  Germany$log_E = log( Germany$E )  # offset on log scale

  sppoly = Germany  # construct "sppoly" with required attributes
  attributes(sppoly)[["areal_units_fn"]] = g  
  attributes(sppoly)[["nb"]] = inla.read.graph(g) 

  p = list(
    project_name = "test",  # basic run
    modeldir = tempdir(),
    carstm_model_label = "testlabel",
    variabletomodel = "Y",
    vnS = "region.struct", 
    aegis_dimensionality = "space",   # a pure space model
    carstm_modelengine = "inla",
    carstm_model_formula = Y ~   f(region.struct, model="bym2", graph.file=g)  + f(x, model="rw2"),
    carstm_model_family = "poisson",
    nposteriors=5000
  ) 

  fn_fit = tempfile( pattern="fit", tmpdir=p$modeldir )
  fn_res = tempfile( pattern="res", tmpdir=p$modeldir )

  res = carstm_model( 
    p=p,
    M = Germany, 
    E = Germany$E,
    sppoly = sppoly,
    region.id = as.character(Germany$region),
    fn_fit = fn_fit,
    fn_res = fn_res,
    verbose=TRUE
  ) 


     # extract results
  fit3 = carstm_model( p=p, fn_fit = fn_fit, DS="carstm_modelled_fit" )  # extract currently saved model fit
  fit3$summary$dic$dic
  fit3$summary$dic$p.eff
  summary(fit3)  # identical to fit2

  

  plot(fit3)
  plot(fit3, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

  plot(fit$summary.random$region.struct$mean ~ fit3$summary.random$region.struct$mean[1:544] )

  res = carstm_model( p=p, fn_res = fn_res, DS="carstm_modelled_summary"  ) # to load currently saved results

    # plot random effects
  dev.new()
  par( mfrow=c(2,2))
  Bym.map(fit$summary.random$region.struct$mean)
  Bym.map(fit2$summary.random$region.struct$mean[1:nrow(Germany)])  # bym2 concatenates iid and bym effects 
  Bym.map(fit3$summary.random$region.struct$mean[1:nrow(Germany)])  # bym2 concatenates iid and bym effects 

    # plot predictions
  dev.new()
  par( mfrow=c(1,2))
  Bym.map(fit3$summary.fitted.values$mean)  # bym2 concatenates iid and bym effects 
  Bym.map(res$predictions[, "mean"])  # bym2 concatenates iid and bym effects 

  plot( fit3$summary.fitted.values$mean ~ res$predictions[, "mean"] )

    # maps of some of the results  .. can't map unitl polygons are coverted to sf format (TODO)
  tmout = carstm_map(  res=res, vn=c( "predictions" ), 
    plot_elements=c( "compass", "scale_bar", "legend" ),
    main=paste( "Predictions")  
  )

  tmout = carstm_map(  res=res, vn=c( "random", "spatial", "combined" ), 
    plot_elements=c( "compass", "scale_bar", "legend" ),
    main=paste( "Spatial effects")  
  )

