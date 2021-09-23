library(INLA)

data(Seeds)


fitc = glm ( r ~ x1 + offset(log(n)), 
  data=Seeds, 
  family="poisson" 
)
summary(fitc)
Coefficients:
             Estimate Std. Error   z value Pr(>|z|)
(Intercept) -0.636577   0.057735 -11.02584  < 2e-16
x1          -0.119091   0.106761  -1.11549  0.26464

fit1 = glm ( r ~ x1, offset = log(n),
  data=Seeds, 
  family="poisson" 
)
summary(fit1)
Coefficients:
             Estimate Std. Error   z value Pr(>|z|)
(Intercept) -0.636577   0.057735 -11.02584  < 2e-16
x1          -0.119091   0.106761  -1.11549  0.26464

fit2 = inla ( r ~ x1 + offset(log(n)), 
  data=Seeds, 
  family="poisson"  ,
  control.compute = list(config = TRUE), 
  control.predictor = list(compute=TRUE, link=1)
)
summary(fit2)
Fixed effects:
              mean    sd 0.025quant 0.5quant 0.975quant   mode kld
(Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0


fit3 = inla ( r ~ x1, offset=log(n),
  data=Seeds, 
  family="poisson" ,
  control.compute = list(config = TRUE), 
  control.predictor = list(compute=TRUE, link=1)
)
summary(fit3)
Fixed effects:
              mean    sd 0.025quant 0.5quant 0.975quant   mode kld
(Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0


fit4 = inla ( r ~ x1, E=n,
  data=Seeds, 
  family="poisson",
  control.compute = list(config = TRUE), 
  control.predictor = list(compute=TRUE, link=1)
)
summary(fit4)
Fixed effects:
              mean    sd 0.025quant 0.5quant 0.975quant   mode kld
(Intercept) -0.637 0.058     -0.752   -0.636     -0.525 -0.635   0
x1          -0.119 0.107     -0.331   -0.118      0.088 -0.116   0


# so far, all good!

fit2e = inla ( r ~ x1 + offset(log(n)), 
  data=Seeds, 
  family="poisson"  ,
  control.compute = list(config = TRUE), 
  control.predictor = list(compute=TRUE, link=1),
  inla.mode="experimental"
)
summary(fit2e)


fit3e = inla ( r ~ x1, offset=log(n),
  data=Seeds, 
  family="poisson" ,
  control.compute = list(config = TRUE), 
  control.predictor = list(compute=TRUE, link=1)
)
summary(fit3e)


fit4e = inla ( r ~ x1, E=n,
  data=Seeds, 
  family="poisson",
  control.compute = list(config = TRUE), 
  control.predictor = list(compute=TRUE, link=1)
)
summary(fit4e)

---

plot( predict( fitc, type="response" ) ~ predict( fit1, type="response" ) ) # all good
plot( predict( fitc, type="response" ) ~ fit2$summary.fitted.values$mean ) # all good
plot( predict( fitc, type="response" ) ~ fit3$summary.fitted.values$mean ) # all good
plot( predict( fitc, type="response" ) ~ fit3e$summary.fitted.values$mean ) # all good
plot( predict( fitc, type="response" ) ~ fit4$summary.fitted.values$mean ) # not so good .. why?
plot( predict( fitc, type="response" ) ~ fit4e$summary.fitted.values$mean ) # not so good .. why?



posterior_means = function( x, n=1000 ) {
  pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
  return(rowMeans( pp ))
}

p0 = predict( fitc, type="link" )
p1 = predict( fit1, type="link" )
p2 = posterior_means(fit2) 
p3 = posterior_means(fit3) 
p3e = posterior_means(fit3e) 
p4 = posterior_means(fit4) 
p4e = posterior_means(fit4e) 

plot(p0 ~ p1)
plot(p0 ~ p2)
plot(p0 ~ p3)
plot(p0 ~ p3e)
plot(p0 ~ p4)
plot(p0 ~ p4e)

# ---- 

data(Germany)
g = system.file("demodata/germany.graph", package="INLA")
source(system.file("demodata/Bym-map.R", package="INLA"))
Germany = cbind(Germany,region.struct=Germany$region)
Germany$logE = log(Germany$E)

formula1 = Y ~ f(region.struct,model="besag",graph=g) + f(region,model="iid")
formula2 = Y ~ offset(logE) + f(region.struct,model="besag",graph=g) + f(region,model="iid")
 
# NOTE: all models give same parameter estimates

posterior_means = function( x, n=1000 ) {
  pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
  return(rowMeans( pp ))
}
 

# eta* = A ( eta + offset )
mC2 = inla(formula2, family="poisson", data=Germany, 
    inla.mode="classic", 
    control.compute = list(config = TRUE, return.marginals.predictor=TRUE), 

    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pC2 = posterior_means(mC2)  # range of pC2 larger
plot( pC0_fitted  ~ mC2$summary.fitted.values$mean ) #  ??? not sure what is going on
plot( pC0_fitted  ~ pC2 ) # ??? not sure what is going on

 

mE0 = inla(formula1, family="poisson", data=Germany, E=E, 
    inla.mode="experimental", 
    control.compute = list(config = TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pE0 = posterior_means(mE0)
plot( pC0_fitted ~ mE0$summary.fitted.values$mean ) # fitted.values  on link scale
plot( pC0_fitted ~ pE0 )  # posterior samples on link scale

 
mE2 = inla(formula2, family="poisson", data=Germany, 
    inla.mode="experimental", 
    control.compute = list(config = TRUE), 
    control.fixed=list(prec.intercept=1),
    control.predictor = list(compute=TRUE, link=1)
)
pE2 = posterior_means(mE2)
plot( pC0_fitted ~ mE2$summary.fitted.values$mean ) # fitted.values on link scale
plot( pC0_fitted ~ pE2 )  # ??? not sure what is going on 








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
  control.compute = list(config = TRUE), 
  control.fixed=list(prec.intercept=1),
  control.predictor = list(compute=TRUE, link=1),
  inla.mode="experimental"
)

plot( mC$summary.fitted.values$mean ~ mE$summary.fitted.values$mean    ) # issue: link scale with no offset
plot( mC$summary.fitted.values$mean ~ I(exp(mE$summary.fitted.values$mean)*Germany$E )   )



posterior_means = function( x, n=1000 ) {
  pp = inla.posterior.sample.eval( function() Predictor, inla.posterior.sample( n, x ) )
  return(rowMeans( pp ))
}

pC = posterior_means(mC)
pE = posterior_means(mE)
plot( pC ~ pE  )
plot( mC$summary.fitted.values$mean ~ exp(pC)   )
plot( I(exp(mE$summary.fitted.values$mean)*Germany$E )  ~ exp(pE) )

plot( mC$summary.random$region$mean ~ mE$summary.random$region$mean    ) # issue: link scale with no offset

plot( mC$summary.random$region.struct$mean ~ mE$summary.random$region.struct$mean    ) # issue: link scale with no offset

