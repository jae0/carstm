
# some tests of INLA and CARSTM

# how INLA treats offsets are not always clear or consistent
# CARSTM tries to handle them: the following tests should succeed so long as INLA does not change base behaviour

library(carstm)
# loadfunctions("carstm")

# these are simple regression tests
carstm_test_inla("gaussian")

carstm_test_inla("poisson")

carstm_test_inla("binomial")

# carstm_test_inla("nbinomial")


# ---- 
# testing of slightly more complex, "bym" models

# NOTE key difference:  .. trying to account for this divergent behaviour causes a lot of the complexity within carstm
# classical mode: predictions on user scale, incorporating offects 
# experimental mode:  predictions are on user scale, NOT incorporating offects  


library(INLA)

data(Germany)
g = system.file("demodata/germany.graph", package="INLA")
source(system.file("demodata/Bym-map.R", package="INLA"))
Germany = cbind(Germany,region.struct=Germany$region)
Germany$logE = log(Germany$E)

obsrate = Germany$Y / Germany$E

formula1 = Y ~ f(region.struct,model="besag",graph=g) + f(region,model="iid")
formula2 = Y ~ offset(logE) + f(region.struct,model="besag",graph=g) + f(region,model="iid")
 
# NOTE: all models give same parameter estimates

posterior_means = function( x, n=100 ) {
  pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
  return(rowMeans( exp(pp) ))  # marginals link is log/exp
}
 

# mC2 = "classic" mode for model with offsets
# mE0 = "experimental" mode for model with NO offsets
# mE2 = "experimental" mode for model with offsets

# eta* = A ( eta + offset )
mC2 = inla(formula2, family="poisson", data=Germany, 
    inla.mode="classic", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pC2 = posterior_means(mC2)  # range of pC2 larger
pC2marg  = unlist( sapply( mC2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

# NOTE the varying offset and scale issues:
  plot( Germany$Y ~ mC2$summary.fitted.values$mean  ) #  summary.fitted.values$mean is on user scale with offset .. prediction is count
  plot( Germany$Y ~ pC2 ) # posterior samples are on link scale .. which are then exponentiated and with offsets .. prediction is count
  plot( Germany$Y ~ pC2marg ) # marginals are again on user scale with offsets .. prediction is count
  

# no offets in formula
mE0 = inla(formula1, family="poisson", data=Germany, E=E, 
    inla.mode="experimental", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pE0 = posterior_means(mE0)
pE0marg  = unlist( sapply( mE0$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

# NOTE the varying offset and sacle issues:
  plot( obsrate ~ mE0$summary.fitted.values$mean ) # fitted.values on user scale and ignores offsets for prediction p.. prediction is rate
  plot( obsrate ~ pE0 )  # posterior samples on link scale that are inverted (log/exp) and ignores offsets for prediction .. prediction is rate
  plot( obsrate ~ pE0marg ) # posterior marginals also on user scale and ignores offsets for prediction  .. prediction is rate
  
 
 # has offsets in formula
mE2 = inla(formula2, family="poisson", data=Germany, 
    inla.mode="experimental", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pE2 = posterior_means(mE2)
pE2marg  = unlist( sapply( mE2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

# NOTE the varying offset and sacle issues:
  plot( obsrate ~ mE2$summary.fitted.values$mean ) # fitted.values on user scale and ignores offsets for prediction .. prediction is rate even when offsets given
  plot( Germany$Y ~ pE2 )  # posterior means incorporates given offsets and predicts a count 
  plot( obsrate ~ pE2marg ) #   posterior marginals also on user scale and ignores offsets for prediction  .. prediction is rate
   

  
 # negative binomial .. has offsets in formula
nbin0 = inla(formula2, family="nbinomial", data=Germany, 
    # inla.mode="experimental", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)

pnbin0 = posterior_means(nbin0)
pnbinmarg  = unlist( sapply( nbin0$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )

# NOTE the varying offset and sacle issues:
  plot( obsrate ~ nbin0$summary.fitted.values$mean ) # fitted.values on user scale and ignores offsets for prediction .. prediction is rate even when offsets given
  plot( Germany$Y ~ pnbin0 )  # posterior means incorporates given offsets and predicts a count 
  plot( obsrate ~ pnbinmarg ) #   posterior marginals also on user scale and ignores offsets for prediction  .. prediction is rate
   

  

  # eta* = A ( eta + offset )
mm = inla(formula1, family="poisson", data=Germany, 
    inla.mode="classic",
    offset=Germany$logE,
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)

m2 = inla(formula1, family="poisson", data=Germany, 
    inla.mode="classic",
    E=Germany$E, 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)

> plot(mC2$summary.fitted.values$mean, m2$summary.fitted.values$mean)

> plot(mC2$summary.fitted.values$mean, m2$summary.fitted.values$mean*Germany$E)

