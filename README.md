# carstm

(C)onditional (A)uto(R)egressive (S)pace(T)ime (M)odels
Conditional AutoRegressive ("CAR") models are just about th simplest possible way of accounting for spatial autorcorrelation. It is essentially the analogue of the temporal ARIMA type models but in space. It is arguably even simpler as the units of space are discrete areal units (arbitrary polygons). Originally it was formulated in lattice form and related to the Ising models in physics. Its use is most prevalent in epidemiology made accessible by the original Besag-York-Mollie paper, Leroux and Riebler, and many others.

This project is basically an accounting front-end to the computational engines (INLA, sf) to facillitate the computing, storage and access to the results of these models. Relatively simple to do for a simple spatial model without carstm, however, in spatio-temporal models, a bit of a challenge. This tool is mostly written to facilitate my work flow to estimate relationships and predict in a Bayesian spato-temporal context, variables such as depth (aegis.bathymetry), substrate grain size (aegis.substrate), bottom temperature (aegis.temperature), species composition (aegis.speciescomposition), snow crab numerical abundance, size and probability of observation to estimate viable habitat (bio.snowcrab), and extentions to cod abundance, number, size and habitat (aegis.survey), etc.. Many other variables require spatiotemporal modelling (e.g., physiological condition of fish, system-level metabolism) and are planned.  

The dynamic modelling is limited to basic statistical modelling and is not able to compute latent state space models (yet). Currently use of Julia/Turing/Differential Equations/SciML poses a serious risk of making these really interesting questions viable ...   


In the project aegis.survey, you will find a sequence of examples using *carstm* with supporting functions that examines abundance estimation of snow crab and groundfish (Atlantic cod) in
Maritimes groundfish survey strata. The approach is easily generalizable to all species. The main sequence of example models and scripts are found in inst\scripts\* for the projects: bio.snowcrab and aegis.surveys, and leverages the https::\github.com\jae0\aegis data and GIS handing routines.

For Atlantic cod, see the preprints at: 

[Use in estimating an abundance index of Atlantic Cod](https://doi.org/10.1101/2022.05.05.490753)
 
[Use in decomposing environmental effects for Atlantic Cod](https://doi.org/10.1101/2022.04.21.488963)

[Use in decomposing environmental effects for Snow Crab](https://doi.org/10.1101/2022.12.20.520893)


*carstm* can also replicate the "standard" analysis used by DFO Maritimes (known as "stratanal"), a really basic stratified average estimate. This is shown to be equivalent to a Gaussian linear fixed effects model. Thereafter, a model-based approach is used to incrementally improve upon the assumptions of the model, focussing upon the distributional model (Poisson, overdispersed Poisson), adding environmental covariates and then employing an INLA-based ICAR (intrinsic conditionally autoregressive models; "bym2") approach towards accounting for areal unit modelling and an AR1 temporal autocorrelation assuming separability of the spacetime autocorrelation.


The remainder here is an example of a real use case that fully shows the space-time-seasonal ("cyclic") approach: modelling temperature near Halifax, Nova Scotia, with season discretized to 10 units (to improve data density). Some example data are found in carstm::extdata/aegis_spacetime_test/RDS and is used in the example file: [carstm::scripts/example_temperature_carstm.md](inst/scripts/example_temperature_carstm.md)  
