# CARSTM

Conditional AutoRegressive ("CAR") models are just about th simplest possible way of accounting for spatial autorcorrelation. It is essentially the analogue of the temporal ARIMA type models but in space. It is arguably even simpler as the units of space are discrete areal units (arbitrary polygons). Originally it was formulated in lattice form and related to the Ising models in physics. Its use is most prevalent in epidemiology made accessible by the original Besag-York-Mollie paper, Leroux and Riebler, and many others.

For the Julia implementation, see examples and progression of model building in:

    - [notes_carstm_julia.qmd](./docs/notes_carstm_julia.qmd)
    - [notes_gaussian_process_comparisons.qmd](./docs/notes_gaussian_process_comparisons.qmd)
    - [notes_carstm_julia.qmd](/docs/notes_carstm_julia.qmd)
    
For the R implementation, they are structured to be installed as a package (see below) and documentation in:
    
    - [](./docs/carstm_methods.md)
    - [](./inst/scripts/example_temperature_carstm.md]
    - and various aegis*. and bio.snowcrab packages make heavy use of it. 

## Conditional AutoRegressive Space-Time Models (in R)

The R project is basically an accounting front-end to the computational engines (INLA, sf) to facillitate the computing, storage and access to the results of these models. Relatively simple to do for a simple spatial model without CARSTM, however, in spatio-temporal models, a bit of a challenge. This tool is mostly written to facilitate my work flow to estimate relationships and predict in a Bayesian spatiotemporal context.

[Here is an example of a real use case that fully shows the space-time-seasonal approach: modelling ocean bottom temperature near Halifax, Nova Scotia, with season discretized to 10 units (to improve data density) is here.](inst/scripts/example_temperature_carstm.md)  

As it operates with other aegis.* projects, each of the following are modelled and becomes accessible as covariates: ocean depth (aegis.bathymetry), ocean substrate grain size (aegis.substrate), ocean bottom temperature (aegis.temperature), demersal species composition (aegis.speciescomposition), Snow crab numerical abundance, size and probability of observation to estimate viable habitat (bio.snowcrab), and extentions to cod abundance, number, size and habitat (aegis.survey), etc. Many other variables require spatiotemporal modelling (e.g., physiological condition of fish, system-level metabolism) and are planned for one of these days.  

The dynamic modelling is limited to basic statistical modelling and is not able to compute latent state space models (yet). Currently use of Julia/Turing/Differential Equations/SciML poses a serious risk of making these really interesting questions viable ...   

---

In the project aegis.survey, you will find a sequence of examples using *CARSTM* with supporting functions that examines abundance estimation of snow crab and groundfish (Atlantic cod) in Maritimes groundfish survey strata. The approach is easily generalizable to all species. The main sequence of example models and scripts are found in inst\scripts\* for the projects: bio.snowcrab and aegis.surveys, and leverages the https::\github.com\jae0\aegis data and GIS handing routines.

For Atlantic cod, see the preprints at: 

[Use in estimating an abundance index of Atlantic Cod](https://doi.org/10.1101/2022.05.05.490753)
 
[Use in decomposing environmental effects for Atlantic Cod](https://doi.org/10.1101/2022.04.21.488963)

[Use in decomposing environmental effects for Snow Crab](https://doi.org/10.1101/2022.12.20.520893)

---

Another (less interesting) usage of *CARSTM* is to replicate the "standard" analysis used by DFO Maritimes (known as "stratanal"): a really basic stratified average estimate. This is shown to be equivalent to a Gaussian linear fixed effects model. The utility becomes evident as this model-based approach can be used to incrementally improve upon the assumptions of the basic model, focussing upon the distributional assumptions of the model (Poisson, overdispersed Poisson), adding environmental covariates and then employing an INLA-based ICAR (intrinsic conditionally autoregressive models; "bym2") approach towards accounting for areal unit modelling and an AR1 temporal autocorrelation assuming separability of the spacetime autocorrelation. [See this document for more details and notes.](docs/carstm_methods.pdf)


### Installation in R

To install, run the following:

```r
  remotes::install_github( "jae0/aegis")  # helper functions
  remotes::install_github( "jae0/carstm")
``` 

You probably will want to have an Rprofile set up properly such as:

```r
homedir = path.expand("~")
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

require( aegis )
require( carstm )

```


## Conditional AutoRegressive Space-Time Models (in Julia) and other approaches: 

- CAR/ICAR/BYM areal unit models (space)
- GP models in time or space or covariates ("features")
- Temporal (dynamical) models  
- Spacetime models that combine all of the above using Julia/Turing

Most didactic are the examples that use STAN (see info copied from various sources in "docs/example_car_stan.md"). 

The INLA/R, BRMS/R, Turing/Julia (subdirectory examples) implementations are confirming parameter estimates.

Ultimately, we use Julia/Turing due to simplicity and speed and attempt to develop a useful statistical and dynamical model of snow crab spatiotemporal processes.

It is similar in scope to the github.com/jae0/carstm, however, with no reliance upon INLA for computation. 

### Installation ... clone the repository or copy/download the Julia files


## Useful references

Mitzi Morris: https://mc-stan.org/users/documentation/case-studies/icar_stan.html (very thorough)

Max Joseph: Exact sparse CAR models in Stan: https://github.com/mbjoseph/CARstan 

https://github.com/ConnorDonegan/Stan-IAR
https://github.com/ConnorDonegan/survey-HBM#car-models-in-stan

https://www.mdpi.com/1660-4601/18/13/6856

Besag, Julian, Jeremy York, and Annie Mollié. "Bayesian image restoration, with two applications in spatial statistics." Annals of the institute of statistical mathematics 43.1 (1991): 1-20.

Gelfand, Alan E., and Penelope Vounatsou. "Proper multivariate conditional autoregressive models for spatial data analysis." Biostatistics 4.1 (2003): 11-15.

Jin, Xiaoping, Bradley P. Carlin, and Sudipto Banerjee. "Generalized hierarchical multivariate CAR models for areal data." Biometrics 61.4 (2005): 950-961.

