# carstm
carstm provides a sequence of examples with supporting functions that examines abundance estimation of groundfish (Atlantic cod) in
Maritimes groundfish survey strata. Though focussing upon cod, the approach is easily generalizable to all species. The main sequence of example models and scripts are found in inst\scripts\* and leverages the https::\github.com\jae0\aegis data and GIS handing routines.

Initially, carstm replicates the standard analysis which is known as "stratanal", a basic stratified average estimate. This is shown to be equivalent to a Gaussian linear fixed effects model. Thereafter, a model-based approach is used to incrementally improve upon the assumptions of the model, focussing upon the distributional model (Poisson, overdispersed Poisson), adding environmental covariates and then employing an INLA-based ICAR (intrinsic conditionally autoregressive models; "bym2") approach towards accounting for areal unit modelling and an AR1 temporal autocorrelation assuming separability of the spacetime autocorrelation.

## Example:
  # Prepare data:

  require(INLA)
  data(Germany)
  g = system.file("demodata/germany.graph", package="INLA")
  source(system.file("demodata/Bym-map.R", package="INLA"))
  summary(Germany)
  Germany = cbind(Germany,region.struct=Germany$region)

  fm = Y ~ f(region.struct,model="besag",graph.file=g) + f(region,model="iid") + f(x, model="rw2")

  # Basic BYM model using INLA directly:

  fit =  inla( fm, family="poisson", data=Germany, E=E, verbose=TRUE )


  # Same example using INLA's experiemental mode and empirical Bayes, with marginals

  fit2 =  inla( fm, family="poisson", data=Germany, E=E, verbose=TRUE,
    # control.inla = list( strategy='adaptive' ) , 
    control.inla = list( strategy='adaptive', int.strategy='eb' ),
    control.predictor = list( compute=TRUE),
    control.compute = list( return.marginals.predictor=TRUE ),
    inla.mode="experimental"
  )

  summary(fit)
  summary(fit2)

  dev.new()
  par( mfrow=c(1,2))
  Bym.map(fit$summary.random$region.struct$mean)
  Bym.map(fit2$summary.random$region.struct$mean)


  # Same example using carstm
  require(carstm)

