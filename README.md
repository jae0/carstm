# carstm
carstm provides a sequence of examples with supporting functions that examines abundance estimation of groundfish (Atlantic cod) in
Maritimes groundfish survey strata. Though focussing upon cod, the approach is easily generalizable to all species. The main sequence of example models and scripts are found in inst\scripts\* and leverages the https::\github.com\jae0\aegis and https::\github.com\jae0\aegis.env data and GIS handing routines.

Initially, carstm replicates the standard analysis which is known as "stratanal", a basic stratified average estimate. This is shown to be equivalent to a Gaussian linear fixed effects model. Thereafter, a model-based approach is used to incrementally improve upon the assumptions of the model, focussing upon the distributional model (Poisson, overdispersed Poisson), adding environmental covariates and then employing an INLA-based ICAR (intrinsic conditionally autoregressive models; "bym2") approach towards accounting for areal unit modelling and an AR1 temporal autocorrelation assuming separability of the spacetime autocorrelation.

