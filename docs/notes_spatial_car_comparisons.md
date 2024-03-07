---
title: "Comparison of CAR methods"
header: "CAR methods"
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
	Keywords - Conditional Autoregressive models, Julia, iNLA, brms, STAN
abstract: |
	Compare CAR methods using simple data.
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

## Introduction

Good write up here: 

https://www.multibugs.org/documentation/latest/spatial/SpatialDistributions.html#Appendix


stan implementation:

https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html

https://mc-stan.org/users/documentation/case-studies/icar_stan.html#auto-regressive-models-for-areal-data

https://www.sciencedirect.com/science/article/abs/pii/S1877584518301175

https://github.com/ConnorDonegan/Stan-IAR


Compare across various methods. Space only, using a toy data set:

## INLA (in R)

```r
require(carstm)

# also save to file to have same random seed / numbers as above for Julia
toydata = toy_spatial_data(seed=123, nx=10, fn="~/tmp/toy_spatial_data.RData" )

attach(toydata)  # (data, wb)

str(data)
str(wb)
  
 
require(INLA)

formula = y~ x1 + x2 + f(ID, model="bym2", graph=W, scale.model = TRUE, constr = TRUE, 
  # priors
  hyper = list(theta1 = list("PCprior", c(1, 0.01)),  
  # Pr(sd<1) = 0.01, unlikely to have rr>3just based on the spatial confounding
              theta2 = list("PCprior", c(0.5, 0.5))) 
  # Pr(phi<0.5)=0.5, we believe that the unmeasured spatial confounding
  # is driven 50% from the structured and 50% from the unstructured random effect
)

fit <- inla(formula, data =dat, family="binomial", Ntrials=size,
      control.compute = list( waic=TRUE, dic=TRUE, config=TRUE), # waic=TRUE
      control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
      verbose = F)


summary(fit)
```

~~~
Time used:
    Pre = 4.21, Running = 0.265, Post = 0.114, Total = 4.59 
Fixed effects:
              mean    sd 0.025quant 0.5quant 0.975quant   mode kld
(Intercept) -0.241 0.039     -0.317   -0.241     -0.164 -0.241   0
x1           1.075 0.056      0.966    1.074      1.185  1.074   0
x2           1.047 0.051      0.948    1.047      1.148  1.047   0

Random effects:
  Name	  Model
    ID BYM2 model

Model hyperparameters:
                  mean    sd 0.025quant 0.5quant 0.975quant  mode
Precision for ID 7.938 2.385      4.327    7.581     13.626 6.903
Phi for ID       0.807 0.147      0.437    0.845      0.983 0.947

Deviance Information Criterion (DIC) ...............: 536.51
Deviance Information Criterion (DIC, saturated) ....: 147.54
Effective number of parameters .....................: 45.08

Watanabe-Akaike information criterion (WAIC) ...: 533.26
Effective number of parameters .................: 32.68

Marginal log-Likelihood:  -233.31 
 is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
~~~




## Turing implementation in Julia

This is a simple Julia implementation of a CAR model. 

We use randomly generated data (through R, to show how to use RCall). 


```julia


using RData
using Turing

# this file was saved at the start of his file ( "toy_spatial_data()" )
o = load( joinpath( homedir(), "tmp", "toy_spatial_data.RData" ), convert=true)   # load data; alter file path as required   

dat = o["out"]["dat"]
W = o["out"]["W"]


source_directory = joinpath( homedir(), "bio", "carstm") # alter this to location of car_functions.jl

include( joinpath( source_directory, "startup.jl"  ))  # load support functions

include( joinpath( source_directory, "car_functions.jl"  ))  # load support functions

nb = adjacency_matrix_to_nb(W)
node1, node2, scaling_factor = nodes(nb) # pre-compute required vars from adjacency_matrix outside of modelling step

nAU = length(nb)
nY = size(dat)[1]
nX = 2  #    x1, x2
ysd = std(dat.y / dat.size)  # naive estimate of sd

y = dat.y
Ntrial = floor.(Int, dat.size)
dat.intercept = ones(100)
X = Matrix( dat[:,[ :x1, :x2]] )
 
m = turing_icar_direct_bym2_binomial(y, Ntrial, X, nX, nAU, node1, node2, scaling_factor) 

# ~ 800 sec

n_samples, n_adapts, n_chains = 1000, 1000, 4  
target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.001   

# turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth )
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # sample 

showall( summarize(o) )
```
~~~
  parameters      mean       std      mcse    ess_bulk    ess_tail      rhat   ess_per_sec 
      Symbol   Float64   Float64   Float64     Float64     Float64   Float64       Float64 

       beta0   -0.2411    0.0748    0.0028    718.1295   1185.5090    1.0024        0.7695
    betas[1]    1.0598    0.0597    0.0013   2050.3045   1704.8361    1.0043        2.1970
    betas[2]    1.0409    0.0552    0.0012   2061.9761   1672.1992    1.0025        2.2095
    phi[100]   -0.1663    0.5084    0.0096   2779.5159   1762.7478    1.0047        2.9784
     sum_phi   -0.0024    0.0988    0.0016   3892.9513   2115.3173    1.0011        4.1715
       sigma    0.5311    0.0749    0.0025    921.6516   1935.5351    1.0003        0.9876
         rho    0.9417    0.0828    0.0055    328.6297    343.1122    1.0104        0.3521
  
recovers sigma ~ 0.4
rho : nn when distance <= 1 .. exp(-0.1 *1) = 0.94


Comparing to INLA's solution:

Fixed effects:
             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
(Intercept) -0.452 0.040     -0.531   -0.452     -0.372 -0.452   0
x1           1.101 0.058      0.988    1.101      1.216  1.101   0
x2           1.041 0.053      0.938    1.040      1.145  1.040   0

Model hyperparameters:
                  mean    sd 0.025quant 0.5quant 0.975quant  mode
Precision for ID 7.029 2.012      3.937    6.742     11.790 6.193
Phi for ID       0.806 0.146      0.439    0.844      0.982 0.946


And brms's solutions (below, currently broken)

Note brms's solutions results are slightly different as another random seed is used

Family: binomial 
  Links: mu = logit 
Formula: y | trials(size) ~ x1 + x2 + car(W) 
   Data: dat (Number of observations: 100) 
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Correlation Structures:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
car       0.91      0.11     0.61     1.00 1.02      574     1132
sdcar     0.47      0.10     0.29     0.68 1.01      420      751

Population-Level Effects: 
          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept     0.92      0.13     0.65     1.17 1.04      227      268
x1            1.03      0.05     0.93     1.13 1.00     3272     3225
x2            1.02      0.06     0.91     1.13 1.00     2471     2760

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).

~~~


##  CAR analysis of Scottish Lip cancer, using Julia/Turing

### Small samples to check for model functionality. 

```julia
 
using Turing, PDMats, LinearAlgebra , StatsModels, DataFrames, SparseArrays, Graphs, LazyArrays

source_directory = joinpath( homedir(), "bio", "carstm") # alter this to location of car_functions.jl

include( joinpath( source_directory, "startup.jl"  ))  # load support functions

include( joinpath( source_directory, "car_functions.jl"  ))  # load support functions

# data source:  https://mc-stan.org/users/documentation/case-studies/icar_stan.html
D, W, X, log_offset, y, nX, nAU, node1, node2, scaling_factor = scottish_lip_cancer_data()  # data and pre-computed parameters 
  # y: the observed lip cancer case counts on a per-county basis
  # x: an area-specific continuous covariate that represents the proportion of the population employed in agriculture, fishing, or forestry (AFF)
  # E: the expected number of cases, used as an offset .. log_offset=log.(E),
  # adj: a list of region ids for adjacent regions
  # num: a list of the number of neighbors for each region
  # node1 node2: the nodes for the adjacency matrix
  # scaling factor: re-scaling variance to be equal to 1, using Reibler's solution
auid =1:length(y)

# define sampler params:
# use NUTS: see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
# Morris uses 0.97 for target_acceptance, stan default is 0.95; such high acceptance rate does not work well -- divergent chains
#= for larger runs  
  n_samples, n_adapts, n_chains = 5_000, 1_000, 4
  target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.1   
=#

n_samples, n_adapts, n_chains = 100, 500, 1
target_acceptance, max_depth, init_ϵ = 0.65, 7, 0.05

turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)


# NOTE: param estimates are for short samples ... they will be unstable

m = turing_car(D, W, X, log_offset, y, nX )        # 4 min  -- poorer mixing? (beta)
o = sample(m, turing_sampler, n_samples); showall( summarize(o) )
showall( summarize(o) )

   parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
      Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 

       alpha    0.9647    0.0376    0.0056    30.6683    65.7714    1.0340        0.0799
         tau    1.5766    0.5011    0.0896    28.8793    66.9638    0.9906        0.0752
     beta[1]    0.0492    0.4006    0.1338    10.4083    12.9860    0.9935        0.0271
     beta[2]    0.2394    0.0889    0.0138    44.4606    34.3628    1.0037        0.1158
      phi[1]    1.2524    0.5034    0.1301    16.1530    14.2115    1.0579        0.0421

m = turing_car_prec(D, W, X, log_offset, y, nX )   # 1 min 
o = sample(m, turing_sampler, n_samples); showall( summarize(o) )
showall( summarize(o) )
 
   parameters      mean       std      mcse   ess_bulk   ess_tail      rhat   ess_per_sec 
      Symbol   Float64   Float64   Float64    Float64    Float64   Float64       Float64 

       alpha    0.9550    0.0331    0.0034    60.9872    99.2247    0.9958        0.2404
         tau    1.4834    0.3652    0.1376     7.1995    76.6284    1.1528        0.0284
     beta[1]    0.0428    0.1977    0.1071     3.9234    16.4027    1.1826        0.0155
     beta[2]    0.2143    0.0646    0.0119    31.6640    49.8569    0.9953        0.1248
      phi[1]    1.2456    0.3423    0.0749    17.5361    58.4190    1.0310        0.0691
   
# Morris' "simple_iar" testing difference formulation .. 
# using same run specs: results are similar with much better ess than stan   
m = turing_icar_direct_test( node1, node2 )  # just testing to make sure soft sum contraint is working
o = sample(m, turing_sampler, n_samples); showall( summarize(o) )
showall( summarize(o) )

m = turing_icar_direct_bym(X, log_offset, y, nX, node1, node2 ) # 20 sec
o = sample(m, turing_sampler, n_samples); showall( summarize(o) )
showall( summarize(o) )

     beta[1]   -0.2172    0.6336    0.1480    19.9397    34.3628    0.9919        2.6909
     beta[2]    0.2975    0.1051    0.0107   108.0927    78.3393    1.0019       14.5874
     sum_phi   -0.0035    0.0576    0.0056   109.3275    73.5692    1.0321       14.7540
   tau_theta    4.4587    1.0101    0.1143    75.8465    76.8221    1.0039       10.2357
     tau_phi    2.6273    1.0678    0.1092    90.9362    96.8032    0.9976       12.2721

  
# var(theta) = 1/ 2.4183 = 0.4135136252739528
# var(phi) = 1/4.7161  = 0.21203960899896102
# "rho"(below) ~  sqrt( 0.212 / (0.212 + 0.4135 ) ) = 0.5821

# bym2 requires a "scaling factor"  and auid 
m = turing_icar_direct_bym2(X, log_offset, y, auid, nX, nAU, node1, node2, scaling_factor )   # 35 sec
o = sample(m, turing_sampler, n_samples); showall( summarize(o) )

     beta[1]    0.0797    0.1463     0.0065    0.0145    80.4474    0.9981        2.2399
     beta[2]    0.2970    0.1023     0.0046    0.0067   147.8925    1.0104        4.1177
     sum_phi   -0.0029    0.0542     0.0024    0.0026   481.6455    1.0053       13.4103
       sigma    0.7567    0.1185     0.0053    0.0089   154.5590    0.9984        4.3033
         rho    0.9478    0.0709     0.0032    0.0069    89.4564    0.9999        2.4907

 
```

ESS seem reasonably good even with such small samples. And inference is similar to other methods.

### Larger number of samples 

(requires 254.44 sec!)

```julia

n_samples, n_adapts, n_chains = 500, 1000, 1
target_acceptance, max_depth, init_ϵ = 0.65, 7, 0.05
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

# bym2 requires a "scaling factor"  and auid 
m = turing_icar_direct_bym2(X, log_offset, y, auid, nX, nAU, node1, node2, scaling_factor )   # 10 sec
o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress
showall( summarize(o) )

# if on windows and threads are still not working, use single processor mode:
# o = mapreduce(c -> sample(m, turing_sampler, n_samples), chainscat, 1:n_chains)
 
     beta[1]    0.1295    0.1409     0.0063    0.0135    84.2400    1.0116        3.1415
     beta[2]    0.2932    0.1011     0.0045    0.0098   121.6998    1.0006        4.5385
     sum_phi   -0.0051    0.0579     0.0026    0.0033   317.5338    0.9992       11.8416
       sigma    0.7366    0.0986     0.0044    0.0070   156.0503    1.0068        5.8195
         rho    0.9274    0.1256     0.0056    0.0156    36.0899    1.0488        1.3459

(looks like more samples are required)  

similar but some difference relative to below:
 
Inference for Stan model: 

                mean se_mean   sd  2.5%   50% 97.5% n_eff Rhat
beta0          -0.28    0.00 0.16 -0.61 -0.28  0.04  1943    1
beta1           0.42    0.00 0.16  0.09  0.41  0.74  1819    1
sigma           0.52    0.00 0.09  0.37  0.51  0.70   978 1.00
rho             0.88    0.01 0.13  0.52  0.93  1.00   605 1.01

vs INLA:
 
                  mean        sd 0.025quant   0.5quant 0.975quant       mode          kld
(Intercept) -0.2215948 0.1265029 -0.4711830 -0.2215091 0.02705429 -0.2214959 1.472228e-08
x            0.3706808 0.1320332  0.1054408  0.3725290 0.62566048  0.3762751 4.162445e-09


```

### Grouped areas 

Function not complete: turing_icar_direct_bym2_groups needs to use model matrix construction ...  

Continued in .. ![](notes_gaussian_process_comparisons.qmd)  

```julia
  # bym2 grouped model (multiple groups or disconnected groups): this is not finished 
  groups_unique = unique(sort(groups))
  gi = Vector{Vector{Int64}}()
  for g in groups_unique
      o =  findall(x -> x==g, groups) 
      push!(gi, o)
  end

  scaling_factor = scaling_factor_bym2(node1, node2, groups)  # same function (overloaded)
  m = turing_icar_direct_bym2_groups(X, log_offset, y, auid, nX, nAU, node1, node2, scaling_factor, groups)  # incomplete for now


```

 

## brms (in R)

Fit a CAR model with brms.

NOTE: example or brm seems broken 2023/06/22 .. (skip for now)

```r

library(brms)

fit = brm(y | trials(size) ~ x1 + x2 + car(W, gr=ID), 
           data = dat, data2 = list(W = W),
           family = binomial()) 
summary(fit)
```

~~~
 Family: binomial 
  Links: mu = logit 
Formula: y | trials(size) ~ x1 + x2 + car(W) 
   Data: dat (Number of observations: 100) 
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Correlation Structures:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
car       0.91      0.11     0.61     1.00 1.02      574     1132
sdcar     0.47      0.10     0.29     0.68 1.01      420      751

Population-Level Effects: 
          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept     0.92      0.13     0.65     1.17 1.04      227      268
x1            1.03      0.05     0.93     1.13 1.00     3272     3225
x2            1.02      0.06     0.91     1.13 1.00     2471     2760

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
~~~

```r
# specify autocor terms within the formula
y ~ x + arma(p = 1, q = 1) + car(M)
#> y ~ x + arma(p = 1, q = 1) + car(M)
#> <environment: 0x0000000048121c30>

# specify autocor terms in the 'autocor' argument
bf(y ~ x, autocor = ~ arma(p = 1, q = 1) + car(M))
#> y ~ x 
#> autocor ~ arma(p = 1, q = 1) + car(M)

# specify autocor terms via 'acformula'
bf(y ~ x) + acformula(~ arma(p = 1, q = 1) + car(M))
#> y ~ x 
#> autocor ~ arma(p = 1, q = 1) + car(M)

```

## Greta (in R)

Using Scotish lip data in R

Warning:  installation is a long process.

Documentation:
https://forum.greta-stats.org/t/simulations-with-greta-including-time-series/140
https://forum.greta-stats.org/t/spatial-model-in-greta/47

```r
install.packages("igraph")
install.packages("DiagrammeR")

remotes::install_github("greta-dev/greta")
remotes::install_github("greta-dev/greta.gp")
remotes::install_github("greta-dev/greta.dynamics")

install.packages('pop')
devtools::install_github('goldingn/pop')

require(greta)
require(greta.gp)
require(greta.dynamics)
require(pop)

lip = carstm::scottish_lip_cancer_data()
attach(lip)

X <- as_data(X)
y <- as_data(O)
W <- as_data(W)
log_offset <- as_data(log(E))
D <- as_data(D)

alpha <- uniform(0, 1)
beta <- normal(0, 1, dim = c(2, 1))
tau <- gamma(2, 2)

mu <- t(zeros(N))
sigma <- tau * (D - alpha * W)
phi <- t(multivariate_normal(mu, solve(sigma)))  # solve(sigma) is computing the inverse

distribution(y) <- poisson(exp(X %*% beta + phi + log_offset))    
model <- model(beta, phi, alpha, tau)
draws <- mcmc(model, n_samples = 1000, warmup = 1000, chains = 4, one_by_one = TRUE, n_cores=2 )
summary(draws)
```

~~~
greta solution:
             Mean     SD Naive SE Time-series SE
beta[1,1] -0.0424 0.2596  0.00411        0.02285
beta[2,1]  0.2769 0.0980  0.00155        0.00302
alpha      0.9301 0.0658  0.00104        0.00195
tau        1.6422 0.4986  0.00788        0.01339
 
stan solution:
              mean     se_mean         sd          2.5%          25%
 beta[1]  -0.01792750 0.013038796 0.28861331  -0.631593939  -0.17496830
 beta[2]   0.27432038 0.001415310 0.09442595   0.087170821   0.21124791
 alpha     0.93210676 0.001040790 0.06480737   0.759120479   0.91085712
 tau       1.63267268 0.006709894 0.49358056   0.849935240   1.27696174

Benchmark:
greta: Time difference of 2.937537 mins
Stan: Time difference of 4.895678 mins
 
~~~

# Greta BYM sparse version (in R)

```r

require("carstm")

nb = carstm::adjacency_matrix_to_nb( as.matrix(W))

n = nodes(nb)
node1 = n$node1
node2 = n$node2
scale_factor =  n$scale_factor


alpha <- uniform(0, 1)
beta <- normal(0, 1, dim = c(2, 1))
tau <- gamma(2, 2)
sigma <- tau * (D - alpha * W)

mu <- t(zeros(N))
phi <- t(multivariate_normal(mu, solve(sigma)))  # solve(sigma) is computing the inverse

nnb = length(node1)
zeros <- as_data(rep(0L, nnb))   ## Used to marginalize z out of likelihood
nlk = -log( 0.5 * (phi[node1] - phi[node2])^2 )
distribution(zeros) <- poisson(nlk)

distribution(y) <- poisson(exp(X %*% beta + phi + log_offset))
# distribution(y) <- poisson(exp(X %*% beta + phi[id_loc] + log_offset)) # multiple obs

model <- model(beta, phi, alpha, tau)

draws <- mcmc(model, n_samples = 1000, warmup = 1000, chains = 4, one_by_one = TRUE)
summary(draws)

              Mean     SD Naive SE Time-series SE
beta[1,1]  0.04733 0.0213 0.000336        0.00172
beta[2,1]  0.44173 0.0948 0.001500        0.01297
alpha      0.48532 0.0116 0.000184        0.00185
tau        1.29418 0.1670 0.002641        0.02017

 
 
```



## Greta sparse version :

Not tested

```r

library(greta)
library(R6)

# get greta internals
distribution_node <- .internals$nodes$node_classes$distribution_node
as.greta_array <- .internals$greta_arrays$as.greta_array
check_dims <- .internals$utils$checks$check_dims
distrib <- .internals$nodes$constructors$distrib
fl <- .internals$utils$misc$fl
tf_sum <- greta:::tf_sum  # oops, looks like this one didn't make it into .internals!


sparse_car_distribution <- R6Class(
  "sparse_car_distribution",
  inherit = distribution_node,
  public = list(
    initialize = function(tau, alpha, 
                          W_sparse, D_sparse,
                          lambda, dim) {
      
      tau <- as.greta_array(tau)
      alpha <- as.greta_array(alpha)
      W_sparse <- as.greta_array(W_sparse)
      D_sparse <- as.greta_array(D_sparse)
      lambda <- as.greta_array(lambda)
      
      dim <- check_dims(tau, alpha, target_dim = dim)
      super$initialize("sparse_car", dim)
      
      self$add_parameter(tau, "tau")
      self$add_parameter(alpha, "alpha")
      
      self$add_parameter(W_sparse, "W_s")
      self$add_parameter(D_sparse, "D_s")
      self$add_parameter(lambda, "lbd")

    },
    
    tf_distrib = function(parameters, dag) {
      
      alpha <- parameters$alpha
      tau <- parameters$tau
      
      W_sparse <- parameters$W_s
      D_sparse <- parameters$D_s
      lambda <- parameters$lbd
      
      log_prob <- function(x) {
        
        nloc <- length(lambda)
        npair <- nrow(W_sparse)
        
        phit_d <- x * D_sparse
        phit_w <- rep(0, 4)
        for(i in 1:3){
          phit_w[W_sparse[i, 1]] <- phit_w[W_sparse[i, 1]] + 
            x[W_sparse[i, 2]]
          phit_w[W_sparse[i, 2]] <- phit_w[W_sparse[i, 2]] + 
            x[W_sparse[i, 1]]
        }
        
        ldet_terms <- log(fl(1) - alpha * lambda)
        
        res <- fl(0.5) * (nloc * log(tau) + tf_sum(ldet_terms) -
                        tau * (phit_d * phi - alpha *
                                 (phit_w * phi)))
        
       
      }
      list(log_prob = log_prob, cdf = NULL, log_cdf = NULL)
    },
    
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL
  )
)

sparse_car <- function (tau, alpha, 
                        W_sparse, D_sparse,
                        lambda, dim = NULL) {
  distrib("sparse_car", tau, alpha, 
          W_sparse, D_sparse,
          lambda, dim)
}


# a simple example with 4 locations
loc <- data.frame(x=runif(4),
                  y=runif(4))
d_mat <- matrix(0, ncol = 4, nrow=4)
d_mat[lower.tri(d_mat)][dist(loc) < 0.4] <- 1
d_mat[upper.tri(d_mat)] <- t(d_mat)[upper.tri(d_mat)]
# make W_sparse, D_sparse and lambda
nloc <- nrow(loc)
npair <- sum(d_mat) / 2
counter <- 1
W_sparse <- matrix(0, nrow= npair, ncol = 2)
for (i in 1:(nloc - 1)) {
  for (j in (i + 1):nloc) {
    if (d_mat[i, j] == 1) {
      W_sparse[counter, 1] = i
      W_sparse[counter, 2] = j
      counter = counter + 1
    }
  }
}

# number of neighbour per site
D_sparse <- rowSums(d_mat)
# compute eigenvalues
invsqrtD <- base::diag(1 / sqrt(D_sparse), nrow = nloc, ncol = nloc)
lambda <- base::eigen(t(invsqrtD) %*% d_mat %*% invsqrtD)$values

# the parameter
alpha <- uniform(0, 1)
tau <- gamma(2, 2)

phi <- sparse_car(tau, alpha, W_sparse, D_sparse, lambda, dim = 4)

mm <- model(alpha, tau)
dd <- mcmc(mm)


```

## STAN version

REFERENCE is:
https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/geobugs12manual.pdf
https://www.paulamoraga.com/book-geospatial/sec-arealdataexamplespatial.html


```r

# library(SpatialEpi)
# data(scotland)
# names(scotland)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Define MCMC parameters 
niter <- 1E4   # definitely overkill, but good for comparison
nchains <- 4

# load datamessage("Scottish lip cancer data" )
message("source is https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html" )

lip = carstm::scottish_lip_cancer_data()
attach(lip)

W <- A # adjacency matrix
scaled_x <- c(scale(x))
X <- model.matrix(~scaled_x)
	
full_d <- list(
	n = nrow(X),         # number of observations
	p = ncol(X),         # number of coefficients
	X = X,               # design matrix
	y = O,               # observed number of cases
	log_offset = log(E), # log(expected) num. cases
	W = W)               # adjacency matrix

```

and the STAN code:

```stan
    data {
      int<lower = 1> n;
      int<lower = 1> p;
      matrix[n, p] X;
      int<lower = 0> y[n];
      vector[n] log_offset;
      matrix<lower = 0, upper = 1>[n, n] W;
    }
    transformed data{
      vector[n] zeros;
      matrix<lower = 0>[n, n] D;
      {
        vector[n] W_rowsums;
        for (i in 1:n) {
          W_rowsums[i] = sum(W[i, ]);
        }
        D = diag_matrix(W_rowsums);
      }
      zeros = rep_vector(0, n);
    }
    parameters {
      vector[p] beta;
      vector[n] phi;
      real<lower = 0> tau;
      real<lower = 0, upper = 1> alpha;
    }
    model {
      phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
      beta ~ normal(0, 1);
      tau ~ gamma(2, 2);
      y ~ poisson_log(X * beta + phi + log_offset);
    }

```

https://github.com/ConnorDonegan/Stan-IAR
https://github.com/ConnorDonegan/survey-HBM#car-models-in-stan
https://mc-stan.org/users/documentation/case-studies/icar_stan.html


## paper:
https://www.mdpi.com/1660-4601/18/13/6856


---

 
The ICAR prior specified with a binary connectivity matrix reduces to a function of the pairwise differences of neighboring values under the constraint that the parameter vector sums to zero. Morris et al. program the ICAR prior in Stan using the following function:

```stan
/**
  * intrinsic autoregressive prior
  * @return lpdf of IAR prior minus any constant terms
  */
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]) +
      normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
```

where phi is the N-length vector of parameters to which the ICAR prior is assigned. node1 and node2 contain the indices of each pair of connected nodes. It is simplified by keeping the ICAR prior on unit scale. So instead of assigning the prior directly to phi, the ICAR prior is assigned to a vector of standard normal deviates phi_tilde; then phi = phi_tilde * phi_scale gets passed into the linear predictor of the model.

This function contains two restrictions. The first is that the connectivity structure must consist of binary entries only (ones for neighboring observations, zero otherwise).

The second restriction is that the graph structure needs to be fully connected. For graph structures that are not fully connected, the sum to zero constraint needs to be applied to each connected region separately; as a result, it is best for each connected component to have its own intercept.

Finally, the ICAR prior is typically used in conjunction with a spatially unstructured term theta to capture variation around the local mean (the local mean being modeled by phi.) The BYM model consists of the combination of local and global partial-pooling (so-called random effects terms). This is refered to as the convolution term, convolution = phi + theta.

The following Stan function calculates the log probability of the ICAR prior, adjusting as needed for use with the BYM model. In short, observations with zero neighbors will be handled differently depending on the inclusion of theta; if the model has theta, then those phi values need to drop out. If the model does not have theta, this code assigns to the zero-neighbor observations an independent Gaussian prior with scale equal to phi_scale.

```stan
/**
 * Log probability of the intrinsic conditional autoregressive (ICAR) prior,
 * excluding additive constants.
 *
 * @param phi Vector of parameters for spatial smoothing (on unit scale)
 * @param spatial_scale Scale parameter for the ICAR model
 * @param node1
 * @param node2
 * @param k number of groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param has_theta If the model contains an independent partial pooling term, phi for singletons can be zeroed out; otherwise, they require a standard normal prior. Both BYM and BYM2 have theta.
 *
 * @return Log probability density of ICAR prior up to additive constant
 **/
real icar_normal_lpdf(vector phi, real spatial_scale,
              int[] node1, int[] node2,
              int k, int[] group_size, int[] group_idx,
              int has_theta) {
  real lp;
  int pos=1;
  lp = -0.5 * dot_self(phi[node1] - phi[node2]);
  if (has_theta) {
    for (j in 1:k) {
      /* sum to zero constraint for each connected group; singletons zero out */
      lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * group_size[j]);
      pos += group_size[j];
    }
  } else {
    /* does not have theta */
    for (j in 1:k) {
      if (group_size[j] > 1) {
    /* same as above for non-singletons: sum to zero constraint */
    lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * group_size[j]);
      } else {
    /* its a singleton: independent Gaussian prior on phi */
    lp += normal_lpdf(phi[ segment(group_idx, pos, group_size[j]) ] | 0, spatial_scale);
      }
      pos += group_size[j];
    }
  }
  return lp;
}
```

---

## Riebler parameterization: BYM convolution term

Riebler et al. (2016) proposed an adjustment to the ICAR model to enable more meaningful priors to be placed on phi_scale. The idea is to adjust the scale of phi for the additional variance present in the covariance matrix of the ICAR model, relative to a covariance matrix with zeroes on the off-diagonal elements. This is introduced through a scale_factor term, which we will code as inv_sqrt_scale_factor = sqrt(1/scale_factor) (to relate this directly to other implementations you may find).

The following function is used to combine phi_tilde with phi_scale as well as the scale_factor (which may be a vector ones, to be ignored).

```stan
/**
 * Create phi from phi_tilde, inv_sqrt_scale_factor, and spatial_scale.
 *
 * @param phi_tilde local component (spatially autocorrelated)
 * @param phi_scale scale parameter for phi
 * @param inv_sqrt_scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA);
 *                              transformed from 1/scale^2 --> scale. Or, a vector of ones.
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 *
 * @return phi vector of spatially autocorrelated coefficients
 */
vector make_phi(vector phi_tilde, real phi_scale,
              vector inv_sqrt_scale_factor,
              int n, int k,
              int[] group_size, int[] group_idx
              ) {
  vector[n] phi;
  int pos=1;
  for (j in 1:k) {
      phi[ segment(group_idx, pos, group_size[j]) ] = phi_scale * inv_sqrt_scale_factor[j] * phi_tilde[ segment(group_idx, pos, group_size[j]) ];
    pos += group_size[j];
  }
  return phi;
}
```

One way this model can be extended is by assigning a separate scale parameters for each connected component of the graph. Once you assign separate intercepts and scale parameters for each disconnected region, you have independent prior models for each region. Imposing the constraint that disconnected regions have the same scale parameter may seem unreasonable for some applications, and also may slow down sampling. For example, why would the spatial autocorrelation parameters for the counties of Hawaii have the same scale as those for the continental U.S.?

Implementing this extension requires declaring vector<lower=0>[1+m] phi_scale in the parameters block and then adjusting the make_phi function as follows:

```stan
phi[ segment(group_idx, pos, group_size[j]) ] = phi_scale[j] * inv_sqrt_scale_factor[j] * phi_tilde[ segment(group_idx, pos, group_size[j]) ];`.
```

The icar-functions.stan file contains a function called make_phi2 with that adjustment made.



The BYM model includes the parameter vector assigned the ICAR prior plus a vector theta assigned a normal prior with unknown scale: θ ∼ N(0, η), with η assigned some prior such as η ∼ N(0, 1). Again, in practice we assign theta_tilde a standard normal prior and then multiply it by its scale theta_scale. Then the convolution term is

convolution = phi + theta = phi_tilde * spatial_scale + theta_tilde * theta_scale


or optionally with the scaling factor:

convolution = phi_tilde * inv_sqrt_scale_factor * spatial_scale + theta_tilde * theta_scale

The following function combines terms to create the BYM convolution term, making adjustments as needed for disconnected graph structures and observations with zero neighbors. The input for phi should be the parameter vector returned by make_phi (as demonstrated below).

```stan
/**
 * Combine local and global partial-pooling components into the convolved BYM term.
 *
 * @param phi spatially autocorrelated component (not phi_tilde!)
 * @param theta global component (not theta_tilde!)
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 *
 * @return BYM convolution vector
 */
vector convolve_bym(vector phi, vector theta,
              int n, int k,
              int[] group_size, int[] group_idx
              ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:k) {
     if (group_size[j] == 1) {
        convolution[ segment(group_idx, pos, group_size[j]) ] = theta[ segment(group_idx, pos, group_size[j]) ];
    } else {
    convolution[ segment(group_idx, pos, group_size[j]) ] =
      phi[ segment(group_idx, pos, group_size[j]) ] + theta[ segment(group_idx, pos, group_size[j]) ];
  }
      pos += group_size[j];
  }
  return convolution;
}
```
 
Riebler et al. (2016) also proposed to combine theta with phi using a mixing parameter rho and a single scale spatial_scale, such that

convolution = spatial_scale * (sqrt(rho * scale_factor^-1) * phi + sqrt(1 - rho) * theta)

The following function creates the convolution term for the BYM2 model and makes adjustments for disconnected graph structures and zero-neighbor observations. If you use this function, do not also use the make_phi function (see BYM2.stan).

```stan
/**
 * Combine local and global partial-pooling components into the convolved BYM2 term.
 *
 * @param phi_tilde local (spatially autocorrelated) component
 * @param theta_tilde global component
 * @param spatial_scale scale parameter for the convolution term
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param inv_sqrt_scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA);
 *                              transformed from 1/scale^2 --> scale. Or, a vector of ones.
 * @param rho proportion of convolution that is spatially autocorrelated
 *
 * @return BYM2 convolution vector
 */
vector convolve_bym2(vector phi_tilde, vector theta_tilde,
          real spatial_scale,
              int n, int k,
              int[] group_size, int[] group_idx,
              real rho, vector inv_sqrt_scale_factor
              ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:k) {
    if (group_size[j] == 1) {
        convolution[ segment(group_idx, pos, group_size[j]) ] = spatial_scale * theta_tilde[ segment(group_idx, pos, group_size[j]) ];
    } else {
    convolution[ segment(group_idx, pos, group_size[j]) ] = spatial_scale * (
     sqrt(rho) * inv_sqrt_scale_factor[j] * phi_tilde[ segment(group_idx, pos, group_size[j]) ] +
     sqrt(1 - rho) * theta_tilde[ segment(group_idx, pos, group_size[j]) ]
      );
  }
  pos += group_size[j];
  }
  return convolution;
}
```
 
The ICAR model requires the following input as data:

    n number of observations
    k number of groups of connected observations (i.e. a connected graph with one disconnected island has k=1)
    group_size an integer array with the number of observations per group
    n_edges total number of edges (sum of all edge_size)
    node1, node2 these contain the indices for each pair of connected nodes
    group_idx integer array (size n), see below. Allows us to extract observations by their group membership
    inv_sqrt_scale_factor a k-length vector of scale factors, one per connected group. Or, a k-length vector of ones.
    m number of connected graph components requiring their own intercept
        for a single fully connected graph, this is zero; add one to m for each additional connected component with group_size > 1.
    A an n by m matrix of dummy variables indicating to which connected component an observation belongs (for any extra intercepts).

The demonstration shows how to use some R code to very easily obtain all these items at once. The Stan code appears more complicated than some applications will require, but it is designed to function under a variety of common circumstances.

The Stan code below also includes the following terms, just to provide a complete example:

    prior_only Binary indicator, skip the likelihood and just compute the prior or sample from posterior distribution given data
    y outcome data (integer array) (ignored if prior_only=1)
    offset offset term (ignored if prior_only=1).

This is for the BYM model:

```stan

// The BYM model //
functions {
#include icar-functions.stan
}
data {
  int n;    // no. observations
  int<lower=1> k; // no. of groups
  int group_size[k]; // observational units per group
  int group_idx[n]; // index of observations, ordered by group
  int<lower=0> m; // no of components requiring additional intercepts
  matrix[n, m] A; // dummy variables for any extra graph component intercepts
  int<lower=1> n_edges;
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges];
  int<lower=1, upper=k> comp_id[n];
  vector[k] inv_sqrt_scale_factor; // can be a vector of ones, as a placeholder
  int<lower=0, upper=1> prior_only;
  int y[n];
  vector[n] offset; // e.g., log of population at risk
}

transformed data {
  int<lower=0,upper=1> has_theta=1;
}

parameters {
  real alpha;
  vector[m] alpha_phi;
  vector[n] phi_tilde;
  real<lower=0> spatial_scale;
  vector[n] theta_tilde;
  real<lower=0> theta_scale;
}

transformed parameters {
  vector[n] phi = make_phi(phi_tilde, spatial_scale, inv_sqrt_scale_factor, n, k, group_size, group_idx);
  vector[n] theta = theta_tilde * theta_scale;
  vector[n] convolution = convolve_bym(phi, theta, n, k, group_size, group_idx);
  vector[n] eta = offset + alpha + convolution;
  if (m) eta += A * alpha_phi;
}

model {
   // keep the following lines as they are:
   phi_tilde ~ icar_normal(spatial_scale, node1, node2, k, group_size, group_idx, has_theta);
   theta_tilde ~ std_normal();
   // the rest of the priors may need to be adjusted for your own model.
   spatial_scale ~ std_normal();
   theta_scale ~ std_normal();
   alpha ~ normal(0, 10); // this is the prior for the mean log rate.
   if (m) alpha_phi ~ normal(0, 2);
   if (!prior_only) y ~ poisson_log(eta);
}
```
 

The following are R codes to run the BYM2 model is in the BYM2.stan file.


Scaling the ICAR prior

To follow Reibler et al.’s adjustment to the scale of the model, you can use the INLA R package and the following R code:

```r
icar.data <- prep_icar_data(C)

## calculate the scale factor for each of k connected group of nodes, using the scale_c function from M. Morris
k <- icar.data$k
scale_factor <- vector(mode = "numeric", length = k)
for (j in 1:k) {
  g.idx <- which(icar.data$comp_id == j)
  if (length(g.idx) == 1) {
    scale_factor[j] <- 1
    next
  }
  Cg <- C[g.idx, g.idx]
  scale_factor[j] <- scale_c(Cg)
}

## update the data list for Stan
icar.data$inv_sqrt_scale_factor <- 1 / sqrt( scale_factor )


```

