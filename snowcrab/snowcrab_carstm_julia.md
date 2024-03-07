---
title: "CARSTM in Julia for snow crab"
header: "CARSTM Julia Snow crab"
date: \today
journal: 'biorxiv'
author:
- name: Jae S. Choi
  email: choi.jae.seok@gmail.com
#  footnote: 1
  orcid: 0000-0003-3632-5723 
#  corresponding: true
affiliation:
#- number: 1
#  name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
keyword: |
	Keywords - Guassian Process / CAR 
abstract: |
	Compare Guassian Process / CARSTM with snow crab data.
unnumberedsections: true
citeproc: true
citeproc-method: biblatex
biblatex: true
biblio-style: authoryear
biblatexoptions: |
 \usepackage[authordate, maxcitenames=1, uniquename=mininit, backend=biber, natbib]{biblatex-chicago}
# csl-refs: true
csl: media/chicago-author-date.csl
# csl: "https://raw.githubusercontent.com/citation-style-language/styles/master/harvard-anglia-ruskin-university.csl"
# csl from: Zotero Style Repository https://www.zotero.org/styles
# https://github.com/citation-style-language/styles  # to get more csl's
# https://www.overleaf.com/learn/latex/Biblatex_citation_styles
bibliography: references.bib
#standalone: true
acknowledgements: |
additionalinformation: |
documentclass: paper
papersize: letterpaper
fontsize: 11pt
output:
	pdf_document:
# remainder are Quarto options
toc: true
number-sections: true
highlight-style: pygments
editor:
  render-on-save: false
format:
  html: 
    code-fold: true
    html-math-method: katex
    embed-resources: true
  pdf:
    pdf-engine: lualatex
  docx: default 
---



<!-- This is a Markdown/Quarto/pandoc-latex document -->

<!-- 
As a Quarto doc.. copy this file to a work directory (e.g., ~/tmp/ )  and run Quarto from there:

# quarto render *.qmd --to html 

Can add "--to docx --to pdf" as additional documents, but their formatting is awkward and will require more work.  

## Examples

```r
#| eval: true
#| output: false
#| warning: false
#| error: false
# code evaluation
```


```r
#| eval: true
#| output: true
#| label: predator_list
#| tbl-cap: "All predators of snow crab on Scotian Shelf of Atlantic Canada. Of 58287 finfish stomach samples, 159 had snow crab (0.28%). There is no indormation on snow crab diet in the database."
# table
```

```{r}
#| eval: false
#| output: false
#| label: temp_depth
#| fig-cap: "Temperature and depths of snow crab predation on the Scotian Shelf of Atlantic Canada. Grey is all species observations in diet data base. Red is snow crab as prey."
#| fig-dpi: 144
#| fig-height: 4
# figure
```

-->


<!--
pandoc/latex template related ...to get started:

Write paper in as *.md, refs in references.bib and use my_template.tex as a pandoc template
Modify makefile as required

citations: are cited as \cite{mittner2014brain} 

equations: 
	referenced as:  eq. (\ref{1}) 
	tagged as:  $$ \pi=\pi \tag{1a} \label{1}$$

pictures:
	referenced as: figure \ref{fig2}
    tagged as: ![my figure captions. \label{fig2}](Figures\Tolmukapea.jpg){width=100px height=50px} 
    
    no working: fignos ..  Fig. @fig:dummy{#fig:dummy width=40% height=20%}

tables: 
	do same as pictures or
	direct latex Tables: 
	
	\begin{table}[ht]
	\centering
	\caption{Probability to observe Bayes Factors of a certain magnitude or above for the used sample-size of $N=60$ assuming the original and the null-hypothesis.}
	\begin{tabular}{llrrr}
	  & & \multicolumn{3}{l}{$P(\text{BF}\ge\theta)$}\\
	  Hypothesis & BF Type & $\theta=3$ & $\theta=10$ & $\theta=20$ \\
	  \hline
	  $d\sim \mathcal{N}(1.57, 0.51)$ & JZS BF$_{10}$ & 0.98 & 0.97 & 0.96 \\
		 & Replication BF$_{10}$ & 0.98 & 0.96 & 0.96 \\
		 & Meta-Analysis BF$_{10}$ & 0.99 & 0.99 & 0.99 \\\cline{2-5}
		$d=0$ & JZS BF$_{01}$ & 0.81 & 0.00 & 0.00 \\
	   & Replication BF$_{01}$ & 0.98 & 0.95 & 0.91 \\
		 & Meta-Analysis BF$_{01}$ & 0.63 & 0.27 & 0.06 \\
	   \hline
	\end{tabular}
	\label{tab:probbf}
	\end{table}

-->



<!--
# Example Makefile for operating with pandoc/latex:
FILENAME=thermodynamics_onsager
TEMPLATE=my_template.tex
PDFENGINE=lualatex

# export TEXINPUTS=.:media//:
# export BIBINPUTS=.:media//:
# export BSTINPUTS=.:media//:

all: latex pdf 

latex: 
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FILENAME).tex \
	--pdf-engine=$(PDFENGINE)

pdf: latex  
	pdflatex $(FILENAME).tex
	biber $(FILENAME) 
	pdflatex $(FILENAME).tex
	pdflatex $(FILENAME).tex

view: 
	zathura $(FILENAME).pdf 

test: 
	pandoc -s $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FILENAME).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FILENAME).tex

html:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FILENAME).html \
	--mathjax \
	--self-contained

epub:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FILENAME).epub \
	--epub-cover-image=<cover-image> \
	--toc

docx:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FILENAME).docx \


odt:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=odt \
	--output=$(FILENAME).odt \

  
git:
	git commit -m"update"
	git checkout master
	git merge develop
	git push
	git pull
	git checkout develop
	git merge master
	git status
	

 
watch: $(FILENAME).md 
	fswatch -o $^ | xargs -n1 -I{} make

.PHONY: clean all

clean:
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out  *.run.xml *.spl  *.docx *.odt *.epub *.html
	

-->


<!-- This is a Markdown/Quarto document -->

<!-- 
Copy this file to a work directory (e.g., ~/tmp/ ) 
and run Quarto from there:

# quarto render *.qmd --to html 

Can add "--to docx --to pdf" as additional documents, but their formatting is awkward and will require more work.  
-->


# Space-time model with snow crab

## NOTES::
## car method works nicely
## see: regression_functions.jl: example_data() and featurize_poly() to work withing model design matrix 
## GP code works in example_3b_gaussian_processes_julia_test.md ... use this for covariates ..
## need to add method for fixed effects and AR1 .. then CARSTM should be complete



```R
# create data from snow crab survey

  source( file.path( "~", ".Rprofile" )  )
  source( file.path( code_root, "bio_startup.R" )  )
   
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 

  year.assessment = 2023
  
  yrs = 1999:year.assessment
  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
  
  runlabel= paste( "1999_present", snowcrab_filter_class, sep="_" )

  # params for number
  pN = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    carstm_model_label= runlabel,  
    selection = list(
      type = "number",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  # params for mean size .. mostly the same as pN
  pW = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,   
    areal_units_type="tesselation",
    carstm_model_label= runlabel,  
    selection = list(
      type = "meansize",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )

  )

  # params for probability of observation
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    carstm_model_label= runlabel,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )

  redo_data = FALSE
  
  if (redo_data) {
    xydata = snowcrab.db( p=pN, DS="areal_units_input", redo=TRUE )
    sppoly = areal_units( p=pN, xydata=xydata, 
        spbuffer=5,  n_iter_drop=5, redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
    
    plot(sppoly["AUID"])

    sppoly$dummyvar = ""
    xydata = st_as_sf( xydata, coords=c("lon","lat") )
    st_crs(xydata) = st_crs( projection_proj4string("lonlat_wgs84") )

    additional_features = snowcrab_mapping_features(pN)  # for mapping below
  
    tmap_mode("plot")
    
    plt = 
      tm_shape(sppoly) +
        tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5) + 
        tm_shape( xydata ) + tm_sf() +
        additional_features +
        tm_compass(position = c("right", "TOP"), size = 1.5) +
        tm_scale_bar(position = c("RIGHT", "BOTTOM"), width =0.1, text.size = 0.5) +
        tm_layout(frame = FALSE, scale = 2) +
        tm_shape( st_transform(polygons_rnaturalearth(), st_crs(sppoly) )) + 
        tm_borders(col = "slategray", alpha = 0.5, lwd = 0.5)

    dev.new(width=14, height=8, pointsize=20)
    plt
    
    M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly, redo=TRUE )  # will redo if not found
  }



  sppoly = areal_units( p=pN )

  M = snowcrab.db( p=pN, DS="carstm_inputs", sppoly=sppoly  )  # will redo if not found
  M = M[ which(is.finite( M$t + M$z +M$pca1)), ]  # drop missings (72)

  save(M, file="/home/jae/bio/carstm//snowcrab/snowcrab_data.rdata" )

  nb = attributes(sppoly)$nb
  save(nb, file="/home/jae/bio/carstm//snowcrab/snowcrab_nb.rdata" )
  
  sppoly = st_drop_geometry(sppoly)
  save(sppoly, file="/home/jae/bio/carstm//snowcrab/snowcrab_sppoly.rdata" )

#   iq = unique( c( which( M$totno > 0), ip ) )
#   iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
 
#   space_id = sppoly$AUID,
#   time_id =  pN$yrs,
#   cyclic_id = pN$cyclic_levels,
 

```


# Snow crab example: 

## Prepare julia environment and import Rdata files

```julia
 
# source_directory = joinpath(  dirname(@__DIR__()) ) #  same folder as the current file
source_directory = "/home/jae/bio/carstm/"
cd( source_directory )

# basic setup 
include( "startup.jl" )  

# load support functions
include( joinpath( source_directory, "car_functions.jl"  ))  

fndat1 = joinpath( source_directory, "snowcrab", "snowcrab_data.rdata" )
fndat2 = joinpath( source_directory, "snowcrab", "snowcrab_nb.rdata" )
fndat3 = joinpath( source_directory, "snowcrab", "snowcrab_sppoly.rdata" )

# load and unwrap containers
M  = load( fndat1, convert=true)["M"]
nb = load( fndat2, convert=true)["nb"]["nbs"]
sp = load( fndat3, convert=true)["sppoly"]

# indexes for identify preds and obs
ip = findall(M.tag .== "predictions")
io = findall(M.tag .== "observations" )
# io = findall( (M.tag .== "observations") .& (M.year .> 2016) )


nAU = length(nb)
ndata = size(M)[1]
nobs = length(io)
npred = length(ip)
 
tuid = M.year
auid = M.space

# mf = ModelFrame(@formula(totno ~ yr), M)  # DataFrames.jl
# mm = ModelMatrix(mf)
# X = mm.m 

# GP vars
Gvars = ["z", "t", "pca1", "pca2"]
G = Matrix( M[:, Gvars] )
nG = length(Gvars)

# Z-scores or covariates:
G_means = mean(G, dims=1)
G_sds = std(G, dims=1)
for i in 1:nG 
  G[:,i] = ( G[:,i] .- G_means[i] ) ./ G_sds[i]
end
G = G[io,:]
nG = size(G)[2]

# fixed effects:
# using StatsModels

X = modelmatrix( ModelFrame( @formula(totno ~ 1 + year ), M[io,:], contrasts = Dict( :year => StatsModels.EffectsCoding() )))
 
nX = size(X)[2]
nY = length(io)

# covars for GP
density = M.totno[io] ./ M.data_offset[io]
dens = log.( density )
dens = (dens .- mean(dens[isfinite.(dens)]) ) ./ std(dens[isfinite.(dens)])
dens[ findall(x->!isfinite(x), dens)] .= 0 # assume 1 / km^2

log_offset = log.(M.data_offset[io])
y = floor.(Int, M.totno[io])
pa = floor.(Int, M.pa[io])
wt = M.meansize[io]


# defined good habitat aprori == 1
good = findall(x -> x==1, pa)

minimum( y[findall(x -> x==1, pa)]  ) # detection limit == 1
maximum( y[findall(x -> x==0, pa)]  ) # detection limit == 1

# so truncate at [2, Inf]

auid = auid[io]

# adjacency_matrix
node1, node2, scaling_factor = nodes(nb)
W = nb_to_adjacency_matrix(nb)
D = diagm(vec( sum(W, dims=2) ))


```

## run model: CAR in space with covariates (linear base model)  

```julia

# simple spatial form .. year as fixed effect covariates (no GP covars)
 
m = turing_car(D, W, X, log_offset, y, nX )

# testing convergence and timing
o = sample(m, Turing.SMC(), 10)   
 
m = turing_car_prec(D, W, X, log_offset, y, nX )

# testing convergence and timing
o = sample(m, Turing.SMC(), 10)   
 
# simple spatial form .. year as fixed effect covariates (no GP covars)
m = turing_icar_latent_bym2( X, log_offset, y, auid, nX, nY, nAU, node1, node2, scaling_factor) 

rand(m)  # check a sample

# testing convergence and timing
o = sample(m, Turing.MH(), 5000)   
  # 6yr subset: 0.65 sec; sum_phi   -0.3638; sigma    0.5645; rho    0.0729
  # 23yr full: 11.84 sec: 
  parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 

     sum_phi   -0.1618    0.3086    0.0661    24.4525    16.0435    1.2542        2.0644
       sigma    0.6558    0.1954    0.0317    30.7535    23.0754    1.2531        2.5963
         rho    0.4147    0.1782    0.0313    35.0456    17.0560    1.2520        2.9587
     beta[1]    3.5643    0.6587    0.1959    11.5897        NaN    1.2497        0.9785

o = sample(m, Turing.SMC(), 5000)  # 6.5 sec! sum_phi    0.3752; sigma    0.2356; rho    0.2265
  # 23yr full: 14.9 sec: 
5000 samples:
     sum_phi   -0.0076    0.6653    0.0094   5012.1524   4899.3818    0.9999      336.8156
       sigma    0.7957    0.6018    0.0086   4900.8344   4885.7508    1.0002      329.3350
         rho    0.4975    0.3527    0.0051   4728.1707   4679.5472    1.0000      317.7321
     beta[1]   -0.0089    1.0019    0.0141   5072.4783   4915.3341    1.0001      340.8695

10000 samples:
     sum_phi   -0.0016    0.6596    0.0065   10175.0758    9766.0377    0.9999      306.5151
       sigma    0.7914    0.6020    0.0062    9723.5347    9555.2985    1.0007      292.9128
         rho    0.5029    0.3537    0.0036    9368.3785    9656.1857    1.0002      282.2141
     beta[1]   -0.0196    1.0116    0.0101    9999.4833    9714.0759    1.0000      301.2255


# NUTS:: 9.5 hrs for 100 samples of full dataset
     sum_phi   -0.0126    0.6171    0.0474   170.0691   112.1358    1.0511        0.0051
       sigma    1.4281    0.0487    0.0071    49.7466    51.5338    1.0027        0.0015
         rho    0.9880    0.0112    0.0023    22.3970    44.5810    1.0238        0.0007
     beta[1]    6.1779    0.1258    0.0121   106.2940    77.3213    0.9902        0.0032

o = optimize(m, MLE() )  
  # Optim.Options(iterations=5_000, allow_f_increases=true): 
  # does not converge 
  # ...  lp of -11645.64 for 6 yr subset, 
  # full data: sumphi: -0.04039559020440821, sigma 0.7094199293530146, rho 0.00018084535378351768; lp -28326.52

o = optimize(m, MAP() )  
  # 6yr test:  ~ 20 sec; lp of -12885.31; 
  # full data: ~ 5 min: -29629.15
  # 0.004562825353062191, 7.221872596866684, 0.9534728688028509

# ~ 14 min; 
num_elbo_samples, max_iters = (10, 1000)
q = vi(m, ADVI(num_elbo_samples, max_iters) );  # MultivariateDistribution
o = rand(q, 5000);
mean( o, dims=2) #  0.011947792078646548 (sigma_sum); 0.7747005353673827 (sigma);   0.3407988254614214 (rho) ... much better
std( o, dims=2) #  0.7028373853942271 (sigma_sum); 0.17148272838070977 (sigma);  0.04214342974149995  (rho) ... much better

n_samples, n_adapts, n_chains = 100, 100, 4
target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.001   
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

o = sample( m, turing_sampler,  n_samples ) # to see progress
# o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress
 
showall( summarize(o) )

6 yr subset (500, 500): ~ 1 hr (n_samples, n_adapts, n_chains = 100, 100, 4)

  parameters      mean       std   naive_se      mcse         ess      rhat   ess_per_sec 
     sum_phi   -0.0175    0.6777    0.0232    926.1497   364.2469    1.0000        0.2398
       sigma    1.6974    0.0794    0.0079     97.6563   104.1238    0.9999        0.0253
         rho    0.9913    0.0103    0.0033     17.7665    23.9218    1.0999        0.0046
     beta[1]    5.6389    0.1616    0.0542     10.9671    11.6697    1.0702        0.0028

 
full data, (10,10) 191.7 sec

     sum_phi   -0.5064    0.0500    0.0158    10.0000    10.0000    1.1929        0.0521
       sigma    0.4872    0.0075    0.0024     7.4799    10.0000    1.0765        0.0390
         rho    0.5733    0.0131    0.0041     7.6634    10.0000    1.0765        0.0400




# ------------------------
## GP on covars 
 

# one kernel for all GP
m = turing_icar_latent_bym2_gp1( X, G, log_offset, y, dens, auid, nY, nX, nG, nAU, node1, node2, scaling_factor ) 

# product of separate kernels for all GP
m = turing_icar_latent_bym2_gp2( X, G, log_offset, y, dens, auid, nY, nX, nG, nAU, node1, node2, scaling_factor ) 

# product of separate kernels for all GP  snowcrab HURDLE model
# see: https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
m = turing_icar_latent_bym2_gp_hurdle( X, G, log_offset, y, pa, wt, dens, auid, nY, nX, nG, nAU, node1, node2, scaling_factor, good ) 


rand(m)  # check a sample

# testing
o = sample(m, Turing.MH(), 10)

n_samples, n_chains = 10, 2
o = sample(m, Turing.SMC(), n_samples  )  # ~ 25 min (1000 samples) 
# o = sample(m, Turing.SMC(),  MCMCThreads(), n_samples, n_chains )  # ~ 25 min (1000 samples) 

    parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
        Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 

100 samples: subset

       beta[1]    0.0332    0.9460    0.0778   137.7638   116.7774    1.0056        4.5827
       sum_phi   -0.1180    0.7153    0.0754    89.3331   112.1358    0.9937        2.9716
         sigma    0.8525    0.5762    0.0555   100.8281    77.3213    0.9985        3.3540
           rho    0.4835    0.3787    0.0355   103.0927    78.3393    0.9998        3.4293
    kernel_var    0.9693    0.6329    0.0607    98.9947    49.8569    1.0017        3.2930
  kernel_scale    0.2283    0.1656    0.0211    86.8437    51.5563    1.0128        2.8888
         l2reg    0.0021    0.0015    0.0001   132.0515   112.1358    1.0059        4.3926

  full 936.29 sec
       beta[1]   -0.0501    1.0014    0.0926   111.0253   112.1358    1.0041        0.1186
       sum_phi   -0.0055    0.6258    0.0601   103.3219   116.7774    0.9913        0.1104
         sigma    0.8229    0.6264    0.0616   107.6554    78.3393    0.9905        0.1150
           rho    0.4759    0.3472    0.0411    84.9645    78.8353    1.0165        0.0907
    kernel_var    0.8929    0.6772    0.0740    77.1534   101.7705    0.9970        0.0824
  kernel_scale    0.1862    0.1281    0.0123    95.1022    94.4971    1.0028        0.1016
         l2reg    0.0020    0.0013    0.0001    70.3354    57.7252    0.9902        0.0751

1000 samples: subset

       beta[1]   -0.0209    0.9710    0.0328    873.5733    983.2831    1.0022        0.5600
       sum_phi    0.0022    0.6559    0.0218    901.2611    941.3886    0.9993        0.5778
         sigma    0.7674    0.6035    0.0206    861.6209    983.2831    0.9992        0.5524
           rho    0.4979    0.3611    0.0113   1021.1132    942.0177    1.0023        0.6546
    kernel_var    0.9712    0.6793    0.0238    945.6097    841.9459    0.9995        0.6062
  kernel_scale    0.1993    0.1384    0.0047    879.8438    982.7992    0.9993        0.5640
         l2reg    0.0019    0.0014    0.0000   1028.0057   1021.4676    0.9995        0.6590

o = optimize(m, MLE())
# o = optimize(m, MLE(), SimulatedAnnealing(), Optim.Options(iterations=10_000, allow_f_increases=true))

o = optimize(m, MAP())

num_elbo_samples, max_iters = (10, 1000)
# num_elbo_samples, max_iters = (2, 10)  # testing
q = vi(m, ADVI(num_elbo_samples, max_iters) );  # MultivariateDistribution
o = rand(q, 5000)


# n_samples, n_adapts, n_chains = 5, 5, 4
# target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.001   
# turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)





# o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress
# o = sample( m, turing_sampler, init_params = o.values.array) # Sample with the MAP estimate as the starting point.
# o = sample( m, turing_sampler,  n_samples ) # to see progress
 
# o = sample(m, Turing.HMC(0.001, 10), 10)

# o = sample(m, Turing.HMCDA( 5, 0.65, 0.3; ϵ=0.001), 10)

# o = sample(m, Turing.Gibbs( HMC(0.2, 3, :v1), SMC(20, :v2) ), 10)

    SMC: number of particles.
    PG: number of particles, number of iterations.
    HMC: leapfrog step size, leapfrog step numbers.
    Gibbs: component sampler 1, component sampler 2, ...
    HMCDA: total leapfrog length, target accept ratio.
    NUTS: number of adaptation steps (optional), target accept ratio.


    parameters      mean       std   naive_se      mcse          ess      rhat   ess_per_sec 
       beta[1]    0.2189    0.9864     0.3119    0.1836     -13.6352    0.9177       -0.2413
       sum_phi   -0.0747    0.5817     0.1839    0.2148    -449.3469    0.8968       -7.9521
         sigma    0.6706    0.4297     0.1359    0.2152       6.8018    1.1394        0.1204
           rho    0.2819    0.2682     0.0848    0.0954       5.6222    1.1564        0.0995
    kernel_var    1.0288    0.5491     0.1736    0.0600      85.5647    0.9276        1.5142
  kernel_scale    0.1519    0.1194     0.0377    0.0406     -40.4384    0.8944       -0.7156
        lambda    0.0018    0.0012     0.0004    0.0008       8.8172    1.0525        0.1560
           eta    0.2150    1.1391     0.3602    0.3007     -59.0233    0.9895       -1.0445
 
using ReverseDiff, Memoization 
Turing.setadbackend(:reversediff)

n_samples, n_adapts, n_chains = 5, 5, 1

# Morris uses 0.97 for target_acceptance, stan default is 0.95; such high acceptance rate does not work well -- divergent chains
target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.001   

turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

o = sample( m, turing_sampler, n_samples   ) # to see progress

# o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress

showall( summarize(o) )


# p = turing_icar_latent_bym2_predict( o, Xp; Gp=Gp, scaling_factor=scaling_factor, n_sample=10, nAU=length(auid) )






# production run  
n_samples, n_adapts, n_chains = 500, 100, 4

# Morris uses 0.97 for target_acceptance, stan default is 0.95; such high acceptance rate does not work well -- divergent chains
target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.1   


# seems to want:
init_ϵ = 0.001

turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress
# o = sample(m, turing_sampler, n_samples) 

# if on windows and threads are still not working, use single processor mode:
# o = mapreduce(c -> sample(m, turing_sampler, n_samples), chainscat, 1:n_chains)

showall( summarize(o) )

p = turing_icar_latent_bym2_predict( o, Xp; scaling_factor=scaling_fa
# seems to want:
init_ϵ = 0.001

turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress



```
 

 
