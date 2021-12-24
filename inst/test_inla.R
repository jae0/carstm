
# some tests of INLA and CARSTM

# how INLA treats offsets are not always clear or consistent
# CARSTM tries to handle them: the following tests should succeed so long as INLA does not change base behaviour

library(carstm)
# loadfunctions("carstm")

carstm_test_inla("gaussian")

carstm_test_inla("poisson")

carstm_test_inla("binomial")



# ---- 
# NOTE key difference: 
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
  return(rowMeans( exp(pp) ))
}
 

# eta* = A ( eta + offset )
mC2 = inla(formula2, family="poisson", data=Germany, 
    inla.mode="classic", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pC2 = posterior_means(mC2)  # range of pC2 larger
pC2marg  = unlist( sapply( mC2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
plot( Germany$Y ~ mC2$summary.fitted.values$mean  ) #  ??? not sure what is going on
plot( Germany$Y ~ pC2 ) # ??? not sure what is going on
plot( Germany$Y ~ pC2marg ) # ??? not sure what is going on


mE0 = inla(formula1, family="poisson", data=Germany, E=E, 
    inla.mode="experimental", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pE0 = posterior_means(mE0)
pE0marg  = unlist( sapply( mE0$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
plot( obsrate ~ mE0$summary.fitted.values$mean ) # fitted.values  on link scale
plot( obsrate ~ pE0 )  # posterior samples on link scale
plot( obsrate ~ pE0marg ) # ??? not sure what is going on

 
mE2 = inla(formula2, family="poisson", data=Germany, 
    inla.mode="experimental", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pE2 = posterior_means(mE2)
pE2marg  = unlist( sapply( mE2$marginals.fitted.values, function(u) inla.zmarginal(u) )["mean",] )
plot( obsrate ~ mE2$summary.fitted.values$mean ) # fitted.values on link scale
plot( Germany$Y ~ pE2 )  # ??? not sure what is going on 
plot( obsrate ~ pE2marg ) # ??? not sure what is going on
 


data(Germany)
g = system.file("demodata/germany.graph", package="INLA")
source(system.file("demodata/Bym-map.R", package="INLA"))
Germany = cbind(Germany,region.struct=Germany$region)
Germany$logE = log(Germany$E)


# specifying offsets directly in the formula gives strange predictions
mC = inla( Y ~ f(region.struct,model="besag",graph=g) + f(region,model="iid") + offset(logE), 
  family="poisson",   
  data=Germany, 
  control.compute = list(config = TRUE, return.marginals.predictor=TRUE) , 
  control.fixed=list(prec.intercept=1),
  control.predictor = list(compute=TRUE, link=1),
  inla.mode="classic"
)

# specifying offsets directly in the formula gives strange predictions
mE = inla( Y ~ f(region.struct,model="besag",graph=g) + f(region,model="iid") + offset(logE), 
  family="poisson",   
  data=Germany, 
  control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 
  control.fixed=list(prec.intercept=1),
  control.predictor = list(compute=TRUE, link=1),
  inla.mode="experimental"
)

plot( mC$summary.fitted.values$mean ~ mE$summary.fitted.values$mean    ) # issue: link scale with no offset
plot( mC$summary.fitted.values$mean ~ I(exp(mE$summary.fitted.values$mean)*Germany$E )   )



posterior_means = function( x, n=1000 ) {
  pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
  return(rowMeans( exp(pp) ))
}

pC = posterior_means(mC)
pE = posterior_means(mE)
plot( pC ~ pE  )
plot( mC$summary.fitted.values$mean ~ exp(pC)   )
plot( I(exp(mE$summary.fitted.values$mean)*Germany$E )  ~ exp(pE) )

plot( mC$summary.random$region$mean ~ mE$summary.random$region$mean    ) # issue: link scale with no offset

plot( mC$summary.random$region.struct$mean ~ mE$summary.random$region.struct$mean    ) # issue: link scale with no offset



 