---
title: "Gaussian process models"
header: "CAR/GP methods"
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
	Compare Guassian Process / CAR methods using simple data.
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

Needs some clean up and aggregation.


## Gaussian Processes and Kernels in Julia

*Formal definition*: 

    - [Gaussian process](https://en.wikipedia.org/wiki/Gaussian_process) is a collection of random variables (i.e., a stochastic process; $X_t$) 

    - Where any *finite* subset of these random variables has a *joint* multivariate normal distribution (Grimmet & Stirzaker 2001). 

    - They can be represented as one of two equivalent ways (see Rassmussen and Williams, 2006):

        - Function space -- Distribution of functions: GP( f(x), g(x) ); where f and g are mean and covariance functions of x and GP is a MVN .. versitile and expandable

        - Weight-space -- (parameters of regression) projected into a higher dimensional *feature* space via some function (kernel)
 
Limitation: 

    - speed ... computation of a determinant and inverse of a matrix in $R^{n×n}$, generally taking $O(n^3)$ operations to complete, and the computation of the Cholesky decomposition is typically a prerequisite for sampling from a multivariate normal, which also takes $O(n^3)$ time to complete (Golub & Van Loan 1996). Many spatial models, including those which tackle nonstationarity as in Bornn et al. (2012), parametrize the covariance matrix Σ, and hence for Monte Carlo-based inference require repeated recalculations of $Σ^{−1}$ and $|Σ|$.

The *Function space* representations is versatile: 

    - A GP model is simply a set of variables that is multivariate normally distributed with mean $\mu$ and variance/covariance matrix $\Sigma$:

$Y \sim \mathrm{MVN} (\mu, \Sigma)$

    - Usually the mean vector is set to $0$, which means the Gaussian process is fully defined by its choice of variance/covariance matrix $\Sigma$.

    - The variance/covariance matrix is defined by a kernel function which defines the covariance between any two variables:

$\Sigma_{i, j} = K (i, j)$

Kernels are flexible ways of creating/defining the covariance functions that serve as the basis for GPs.




## Kernel methods

kernels: https://www.cs.toronto.edu/~duvenaud/cookbook/

- multiplying two kernels can be thought of as an AND operation. 

    - linear kernel times a periodic results in functions which 
    are periodic with increasing amplitude as we move away from the origin. 

    - A linear kernel times another linear kernel results in functions which are quadratic! 
    This trick can be taken to produce Bayesian polynomial regression of any degree. 

- adding two kernels can be thought of as an OR operation

    - linear kernel plus a periodic results in functions which 
    are periodic with increasing mean as we move away from the origin. 

- An ARD Kernel simply expresses each dimension as being independent from the others. 
    As such Σ is diagonal.




## Preliminaries: Matrix manipulations

```julia
# Aside: see following for matrix specifics wrt julia
# https://docs.juliahub.com/CalculusWithJulia/AZHbv/0.0.5/differentiable_vector_calculus/vectors.html

using Symbolics, LinearAlgebra

a = [10  11  12] # row vector
b = [10; 11; 12] # column vector 

[a; a.+2]   # vertically combine
[b b.+2]   # horizontally combine

M = [3 4 -5; 5 -5 7; -3 6 9] # matrix(3,3)
@variables x1 x2 x3
x = [x1, x2, x3]

M * x - b

@variables c d e f g h i j k
A = [ c d e;  f g h;  i j k ]
A * A'
A * x - b
dot(A, A)

a = [1, 2, 3]
b = [4, 2, 1]
cross(a, b)

# etc ... 
```


## Example problem: 1-D curve fitting

Create some fake data (truth) to model, a one-dimensional toy problem, a fourth-order polynomial:

$f(x)=(x+4)(x+1)(x−1)(x−3)f(x)=(x+4)(x+1)(x−1)(x−3)$


```julia

using KernelFunctions
using LinearAlgebra
using Distributions
using Turing 

using Plots;
default(; lw=2.0, legendfontsize=11.0, ylims=(-150, 500));

using Random: seed!
seed!(42);


# load support functions
# alter this to location of car_functions.jl
code_directory = "/home/jae/bio/carstm/"

include( joinpath( code_directory, "regression_functions.jl"  ))  

Xl, yl, Xo, yo, Xp = example_nonlinear_data() 

plot(Xl, yl; label=raw"$f(x)$ latent 'truth'", legend=:top)  # latent "truth"
scatter!(Xo, yo; seriescolor=:orange, label="observations with error")


```


## 2. Linear regression 

Regression: $y = X \beta + \epsilon$; 

where $y$ is a column vector of length n; $X$ is a matrix (n, k); where k = no covariates + 1 (for intercept if any); $\beta$ is a column vector of length k  (NOTE: it is often referred to as $\beta$, the coefficients of the covariates); and $\epsilon = y - X \beta$ is column vector of length n errors. 

Minimizing $\epsilon$ amounts to least-squares solution where MSE($\beta$) = $n^{-1} \epsilon^{T} \epsilon =  (y - X \beta)^{T}  (y - X \beta)$ which after algebra and minimizing MSE wrt $\beta$ gives:

$X^{T} X \beta - X^{T} y = 0$ 

$\beta = (X^{T} X)^{-1} X^{T} y$  

Predictions $\hat{y}_*$ at test inputs $x_{*}$ are obtained from:  $\hat{y}_* = x_{*}  \beta$

```julia
 
beta = (Xo' * Xo) \ (Xo' * yo)

yp = Xp * beta 
 
plot!(Xp, yp; marker=(:circle,3), label="linear fit")

```

## 3. Adding "features" (higher order polynomial terms)

Fit non-linear data sets by introducing a feature space

We can improve the fit by including additional features
generalizing to:

$\tilde{X} = (ϕ(x_n))_{n=1}^N $​, 

where $\phi(x)$  constructs a feature vector for each input $x$:

$\phi(x) = (1, x, x^2, \dots, x^d)$. 

In this case, it is simply a polynomial regression of degree $d$.

``` julia

function featurize_poly(Xin, degree=1)
    # add higher order polynomials
    return repeat(Xin, 1, degree + 1) .^ (0:degree)'
end

Xl, yl, Xo, yo, Xp = example_nonlinear_data() 

degree = 3
Xof = featurize_poly(Xo, degree)
Xpf = featurize_poly(Xp, degree)
yp = linear_regression(Xof, yo, Xpf ) 
plot!(Xp, yp; width=3, marker=(:star,5),label ="fit of order $degree", legend=:top )

degree = 4
Xof = featurize_poly(Xo, degree)
Xpf = featurize_poly(Xp, degree)
yp = linear_regression(Xof, yo, Xpf )
plot!(Xp, yp; width=3, marker=(:star,5),label ="fit of order $degree", legend=:top )


```

## 4. Ridge regression: regularization


In a higher-dimensional feature space, we can overfit the data; ridge regression introduces regularization to avoid this (L2​ regularization of the weights ($\beta$); Tikhonov regularization): 
 
$\beta=(X^{T} X)^{-1} X^{T} y$ 

we add a ridge parameter λ:

$\beta=(X^{T} X + \lambda 1)^{-1} X^{T} y$ 
 
```julia

function ridge_regression(X, y, Xstar, lambda)
    beta = (X' * X + lambda * I) \ (X' * y)
    return Xstar * beta
end

Xl, yl, Xo, yo, Xp = example_nonlinear_data() 

degree = 3
Xof = featurize_poly(Xo, degree)
Xpf = featurize_poly(Xp, degree)

lambda = 1e-3
yp = ridge_regression(Xof, yo, Xpf, lambda)  

plot!(Xp, yp; width=3, marker=(:triangle,3), label ="fit of order $degree and \$\\lambda=$lambda\$", legend=:top )


degree = 4
Xof = featurize_poly(Xo, degree)
Xpf = featurize_poly(Xp, degree)

lambda = 1e-1
yp = ridge_regression(Xof, yo, Xpf, lambda) 

plot!(Xp, yp; width=3, marker=(:triangle,3), label ="fit of order $degree and \$\\lambda=$lambda\$", legend=:top )

```


## 5. Kernel ridge regression

Instead of constructing the feature matrix explicitly, we can use kernels to replace inner products of feature vectors with a kernel evaluation ... this makes it easy to generalize and extend to/use other kernels ("feature" representations): 

$ \langle \phi (x), \phi (x^T) \rangle = k(x, x^T)$

or 

$\tilde{X} \tilde{X}^{T} = K$, where $K_{ij} = k(x_i,x_j)$.

To apply this "kernel trick" to ridge regression weights ($\beta$):
 
$\beta = ( X^{T} X +\lambda 1)^{-1} X^{T} y$

with the matrix inversion lemma (https://tlienart.github.io/posts/2018/12/13-matrix-inversion-lemmas/index.html)

$\beta = X^{T} ( X X^{T} X +\lambda 1)^{-1} y$

and replacing the inner product ($X X^{T}$) with the kernel matrix ($K$) to give:

$ \beta =X^{T} (K+\lambda 1)^{-1} y$

with predictions:

$\hat{y}_* = x_*  \beta = \langle x_* \beta \rangle = k_* (K + \lambda 1)^{-1} y$

where $(k_*)_n = k(x_*,x_n)$

```julia
 
function kernel_ridge_regression(X, y, Xstar, lambda, kern)
    K = kernelmatrix(kern, X)
    kstar = kernelmatrix(kern, Xstar, X)
    return kstar * ((K + lambda * I) \ y)
end;

Xl, yl, Xo, yo, Xp = example_nonlinear_data() 

# some possible kernels:
kernel = PolynomialKernel(; degree=1, c=1) #  polynomial feature expansion
kernel = PolynomialKernel(; degree=4, c=1)
kernel = SqExponentialKernel() #  infinite-dimensional feature expansion with another kernel
kernel = MaternKernel(ν=0.5)  # Matern 12
kernel = MaternKernel(ν=1.5)  # Matern 32
kernel = MaternKernel(ν=2.5)  # Matern 52

lambda = 1e-3
yp = kernel_ridge_regression(Xo, yo, Xp, lambda, kernel )
plot!(Xp, yp;  marker=(:square,3), label="fit of order $degree and \$\\lambda=$lambda\$")
 
```


## 6. Cholesky and covariance matrix (in kernel ridge regression)

Here are some options for working with matrix inversion (of the covariance matrix).

L = cholesky($\Sigma$);

that is: $\Sigma = L L^{T}$; where L = lower factor (or upper)

For : $\hat{y}_* = X_*  \beta$

and : $\beta = X^{T} \Sigma^{-1} y = X^{T} (K+\lambda 1)^{-1} y$


Usages to solve: 

$\begin{array}{rl}
u = & X^{T} \Sigma^{-1} X \\
= & X^{T} (LL^{T})^{-1} X \\
= & X^{T} L^{-T} L^{-1} X \\
= & (L^{-1} X)^{T} (L^{-1} X)  \\
= & z^{T} z
\end{array}$

$z=L^{-1} X$ needs to be computed just once and $Lz = X$

Or in this case:

$\begin{array}{rl}
\hat{y} = & X^{T} \Sigma^{-1} y \\
= & X^{T} (LL^{T})^{-1} y \\
= & X^{T} L^{-T} L^{-1} y \\
= & (L^{-1}X)^{T} (L^{-1}y)
\end{array}$


```julia

using KernelFunctions 
using LinearAlgebra
using Distributions
using Plots
using BenchmarkTools
 

Xl, yl, Xo, yo, Xp = example_nonlinear_data() 


plot(Xl, yl; label=raw"$f(x)$ latent 'truth'", legend=:top)  # latent "truth"
scatter!(Xo, yo; seriescolor=:orange, label="observations with error")

# here there are 4 possible parameters to tune: lambda, variance compoents (2) and scale 
# this example is about implementation and not tuning so, just picking a few values:

lambda = 1e-4
sigma = 500
scale = 0.1

k = (sigma * SqExponentialKernel() ) ∘ ScaleTransform(scale)
# k = with_lengthscale(sigma * SqExponentialKernel(), 1/scale )  # alternative .. careful with scaling
# k = (sigma * SqExponentialKernel() + sigma2 * Matern32Kernel()) ∘ ScaleTransform(scale)

Ko = kernelmatrix(k, Xo)
Kp = kernelmatrix(k, Xp, Xo)

#  options: (should try other decompositions, but Cholesky is known to be a good choice)

# 0: directly
yp = gp_kernel_ridge_regression_cholesky( sigma, scale, lambda, Xo, Xp, yo )
plot!(Xp, yp;  marker=(:square,3),  label="functional cholesky \$\\lambda=$lambda\$")

@benchmark let
    yp = gp_kernel_ridge_regression_cholesky( sigma, scale, lambda, Xo, Xp, yo )
    # (mean ± σ):   28.534 μs ± 51.886 μs  
end
 
# 1: let julia decide
yp = Kp * ( (Ko + lambda * I ) \ yo )  # note:: A\B == inv(A)*B
# yp = kernelmatrix(k, Xp, Xo) * ((kernelmatrix(k, Xo) + lambda * I) \ yo)
plot!(Xp, yp;  marker=(:triangle,3),   label="direct kernel matrix \$\\lambda=$lambda\$")

@benchmark let
    yp = Kp * ( inv(Ko + lambda * I ) * yo )
    # (mean ± σ):   7.475 μs ±  10.684 μs    # modal
end


# 2: direct inverse

    #= note:: 
        A = Ko + lambda * I 
        B = yo
        A\B == inv(A)*B  
    =#

yp = Kp * ( inv(Ko + lambda * I ) * yo ) 
plot!(Xp, yp;  marker=(:star,3),  label="direct inverse \$\\lambda=$lambda\$")

@benchmark let
    yp = Kp * ( (Ko + lambda * I ) \ yo )  
    #  (mean ± σ):   4.620 μs ±   7.988 μs # long tail
end


# 3: using Cholesky:
# note:: A\B == inv(A)*B

L = cholesky(Ko + lambda * I).L
yp = Kp * ( L' \ (L \ yo) )  
plot!(Xp, yp;  marker=(:star,3),  label="Cholesky -- fastest \$\\lambda=$lambda\$")

@benchmark let
    L = cholesky(Ko + lambda * I)
    yp = Kp * ( L.U \ (L.L \ yo) )  
    #  (mean ± σ):   3.810 μs ±   9.623 μs   <-- winner, but long tailed
end


KP = Kp'
L = cholesky(Ko + lambda * I).L
yp = transpose( L \ KP ) * ( L \ yo )  
plot!(Xp, yp;  marker=(:star,3),  label="Cholesky tranposed \$\\lambda=$lambda\$")

@benchmark let
    L = cholesky(Ko + lambda * I).L
    yp = transpose( L \ KP ) * ( L \ yo )  
    #  (mean ± σ):   16.492 μs ± 29.256 μs   # very slow
end


```


## 7. General MVN samples with Cholesky:

### https://math.stackexchange.com/questions/163470/generating-correlated-random-numbers-why-does-cholesky-decomposition-work
 

Solution: 

$s$ is a random sample

$s = L v + u$;  

with $LL^T = \Sigma = K+\lambda 1$ 

and $v \sim N(0,1)$ and 

$u$ is the mean

$u =  X  \beta$ ; 

$\beta = X^{T} \Sigma^{-1} y$ ; 

$\beta = X^{T} (K+\lambda 1)^{-1} y$ 

$u = X X^{T} (K+\lambda 1)^{-1} y$

so:

$s = L v + X X^{T} (K+\lambda 1)^{-1} y$

$s = L v + X X^{T} (LL^T )^{-1} y$

$s = L v + X X^{T} L^{-T} L^{-1} y$ ; ....  note: inverse distributes right/left

$s = L v + X (X L^{-1})^{T} L^{-1} y$ ; ....  note: inverse distributes right/left

because: 

$\Sigma = E[ (x-u)(x-u)^T ] = A \cdot E[v v^T] \cdot A^T \ = A A^T$


## Same data but with more noise

```julia 

using KernelFunctions
using LinearAlgebra
using Distributions

using Plots;
default(; lw=2.0, legendfontsize=11.0, ylims=(-150, 500));

using Random: seed!
seed!(42);


# same data 
Xl, yl, Xo, yo, Xp = example_nonlinear_data() 
yo = yo .+ rand(Normal(0, 20), length(yo))

plot(Xl, yl; label=raw"$f(x)$ latent 'truth'", legend=:top)  # latent "truth"
scatter!(Xo, yo; seriescolor=:orange, label="observations with more error")
   
lambda = 1e-4
k = (500 * SqExponentialKernel() ) ∘ ScaleTransform(0.1)

Ko  = kernelmatrix(k, Xo)
Kp = kernelmatrix(k, Xl, Xo)
L = cholesky(Ko + lambda * I)

# y_mean_process = Ko * ( inv(Ko + lambda * I ) * yo )
y_mean_process = Ko * ( L.U \ (L.L \ yo) )   # faster equivalent, # note:: A\B == inv(A)*B

plot!(Xp, yp;  marker=(:star,3),  label="Cholesky kernel matrix with more noise \$\\lambda=$lambda\$")
 
# random normal of correct length
v =  rand(Normal(0, 1), length(y_mean_process))
y_random_sample = L.L * v + y_mean_process

plot!(Xp, y_random_sample;  marker=(:star,3),  label="Random samples")
 
```

## Speed tests:


```julia 
using BenchmarkTools

lambda = 1e-4
k = (500 * SqExponentialKernel() ) ∘ ScaleTransform(0.1)

Ko  = kernelpdmat(k, Xo)

L = Ko.chol.L  # kernelpdmat precomputes cholesky
U = Ko.chol.U
# L = cholesky(Ko + lambda * I)

y_mean_process = Ko * ( U \ (L \ yo) )   # faster equivalent, # note:: A\B == inv(A)*B

ny = length(y_mean_process)


@benchmark let
    y_random_sample = L * rand(Normal(0, 1), ny) + y_mean_process
end
# 13.6  μs

@benchmark let
     y_random_sample = rand( MvNormal(y_mean_process, Ko) )
end
# 1.9  μs 

```
  
## ![Autoregressive (AR) processes](https://en.wikipedia.org/wiki/Autoregressive_model)

A set of random variables indexed by time: 

$Y = \{ Y_1, Y_2, \ldots, Y_n \}$ 

An AR model assumes that $Y$ is correlated over time. Defining $Y_t$ in terms of $Y_{t- 1}$ gives:

$Y_t = \rho Y_{t - 1} + \epsilon_t$ 

$\epsilon_t \sim N (0, \sigma_{\epsilon}^2)$ 

$\rho \in R$ is the relationship (correlation) between $Y_t$  and $Y_{t - 1}$

This is an AR process of order 1, AR(1), because $Y_t$ only depends on $Y_{t - 1}$.
  

The AR1 process as a Gaussian Process (Source: [Herb Susmann  (2019)](http://herbsusmann.com/2019/08/09/autoregressive-processes-are-gaussian-processes) )
 
 
Makes them easier to integrate into modeling (space)-time varying processes.

AR models are typically a set of conditional distributions but they can also be represented as a *joint distribution*. 

Rearranging to emphasize that this is the conditional distribution of $Y_t$ given $Y_{t - 1}$:

$\begin{array}{rl}
  Y_t \mid Y_{t - 1} \sim & N (\rho Y_{t-1}, \sigma_{\epsilon}^2)\\
  Y_t \sim & N \left( 0, \frac{\sigma_{\epsilon}^2}{1 -
  \rho^2} \right)
\end{array}$

Where the variance of $Y_t$ comes from the *unconditional* variance (derived below). 

The stationarity condition of an AR process is that
each $Y_t$ has the same distribution; that is,

$\mu = E (Y_i) = E(Y_j)$  and 

$\sigma^2 = \mathrm{Var} (Y_i) = \mathrm{Var} (Y_j)$ for all $i, j$.

The *unconditional mean and variance* of $Y_t$ is then:

$\begin{array}{rl}
  E (Y_t) & = E (\rho Y_{t - 1} + \epsilon_t)\\
  & = \rho E (Y_{t - 1})\\
  \mu & = \rho \mu \text{~(apply~stationarity~condition)}\\
  \mu & = 0\\
  \mathrm{Var} (Y_t) & = \mathrm{Var} (\rho Y_{t - 1} +
  \epsilon_t)\\
  & = \rho^2 \mathrm{Var} (Y_{t - 1}) + \mathrm{Var}
  (\epsilon_t)\\
  & = \rho^2 \mathrm{Var} (Y_{t - 1}) +
  \sigma_{\epsilon}^2\\
  \sigma^2 & = \rho^2 \sigma^2 + \sigma_{\epsilon}^2
  \text{~(apply~stationarity~condition)}\\
  \sigma^2 (1 - \rho^2) & = \sigma_{\epsilon}^2\\
  \sigma^2 & = \frac{\sigma_{\epsilon}^2}{1 - \rho^2}
\end{array}$


$Y$ is jointly normally distributed with some mean vector and variance/covariance matrix.

$E (Y_t) = 0$, so the mean vector of its joint distribution will be $0$.

Covariance between $Y_{t_1}$ and $Y_{t_2}$ (or in this simpler case, $Y_t$ and $Y_{t + 1} )$:

$\begin{array}{rl}
  \mathrm{cov} (Y_t, Y_{t + 1}) & = E [(Y_t - E [Y_t]) (Y_{t+ 1} - E [Y_{t + 1}])] \text{~(definition~of~covariance)~}\\
  & = E [Y_t Y_{t + 1}] \text{~(because~} E [Y_t] = E [Y_{t + 1}] = 0 \text{)}\\
  & = E [Y_t (\rho Y_t + \epsilon_{t + 1})]\\
  & = E [\rho Y_t^2 + Y_t \epsilon_{t + 1}]\\
  & = \rho E [Y_t^2]\\
  & = \rho (\mathrm{Var} (Y_t) + E [Y_t]^2)\\
  & = \rho \frac{\sigma_{\epsilon}^2}{1 - \rho^2}
\end{array}$

for $Y$'s separated by more than one time point, iterating the above
result yields the expression

$\begin{array}{r}
  \mathrm{cov} (Y_{t_1}, Y_{t_2}) = \rho^{\mid t_1 - t_2
  \mid} \frac{\sigma_{\epsilon}^2}{1 - \rho^2}
\end{array}$

The *joint distribution* of $Y$:

$Y \sim \mathrm{MVN} (0, \Sigma)$

where $\Sigma_{i, j} = \rho^{\mid i - j \mid}
\frac{\sigma_{\epsilon}^2}{1 - \rho^2}$. 

This is a Gaussian process.



## Combining kernel functions 

The nice thing about Gaussian processes is that we can combine multiple
kernel functions to model processes with dependence from different
sources. Two ways kernels can be combined are by multiplication and
addition. Multiplying two kernels is like an "AND" operation: the
correlation between points will be high if the correlation from both
kernels is high. Adding two kernels together is like an "OR" operation:
correlation is high if either kernel indicates high covariance.

If the spatial and temporal processes are independent, the "AND" (multiplication) is useful:

As an example, let's build a Gaussian process that combines an AR
process (for temporal correlation) and a spatial process (for spatial
correlation) by combining two kernel functions. 

First, we need to define
an outcome variable $Y$ that varies in time and space: let $Y_{c, t}$ be
a random variable indexed by spatial site $c$ at timepoint $t$. We take
the AR covariance as the first kernel function, to model temporal
correlation:

$K_1 (i, j) = \rho^{\mid t_i - t_j \mid}
\frac{\sigma_{\epsilon}^2}{1 - \rho^2}$

and a squared-exponential kernel function to model spatial dependence:

$K_2 (i, j) = \alpha^2 \exp ~ \left( - \frac{d (i,
j)}{2 \lambda^2} \right)$

where $d (i, j)$ is the spatial distance between sites $i$ and $j$,
$\lambda$ is a length-scale parameter, and $\alpha^2$ is a parameter
controlling the magnitude of the covariance.

Combine the two kernel functions so that two data points are correlated
if they are close together in time and space:

$\begin{array}{rl}
  K (i, j) & = K_1 (i, j) \times K_2 (i, j)\\
  & = \rho^{\mid t_i - t_j \mid}
  \frac{\sigma_{\epsilon}^2}{1 - \rho^2} \alpha^2 \exp
  ~ \left( - \frac{d (i, j)}{2 \lambda^2} \right)
\end{array}$

Note the parameters $\sigma_{\epsilon}^2$ and $\alpha^2$, which are
multipled together, would be unidentifiable in parameter estimation and
should be replaced by a single parameter that controls the magnitude of
the covariance.
 
EG:  Gaussian process using:  

temporal parameters $\rho = 0.9$, $\sigma_{\epsilon}^2 = 1$ and 

spatial parameters $\alpha = 1$ and $\lambda = 2$.


And the spatial distribution over time of $Y_{c, t}$ is shown below:

Visually we can see that the Gaussian process generates data that is
correlated in both time and space.
 

## Modeling using the mean and the covariance  

The spatio-temporal Gaussian process we defined in the previous section
does its modeling through the variance/covariance matrix, with its mean
function set to zero. An alternative way to think about a
spatio-temporal process is akin to the first AR representation we looked
at, and define $Y_t$ (the set of all $Y_{c, t}$ at time $t$t) relative
to $Y_{t - 1}$:

$\begin{array}{r}
  Y_t = \rho Y_{t - 1} + \epsilon_t
\end{array}$

where $\epsilon_t \sim \mathrm{MVN} (0, \Sigma_{\epsilon})$.

If we set $\Sigma_{\epsilon}$ to be the diagonals matrix
$\Sigma_{\epsilon} =
\sigma_{\epsilon}^2 I_n$ then we will have an independent AR(1)
independent process for each spatial site. It gets more interesting if
we define $\Sigma_{\epsilon}$ by a covariance function so we can include
dependence between sites, for example dependence based on the distance
between the sites. For now, let's use the squared exponential kernel and
define $\Sigma_{i, j} =
\alpha^2 \exp ~ \left( - \frac{d (i, j)}{2
\lambda^2} \right)$.

Is this process also equivalent to a mean zero Gaussian process with
some covariance kernel? We'll answer this question by deriving the
covariance between any two points.

The mean of $Y_t$ can be shown to be zero in the same way we showed a
univariate AR process has mean 0. We also need to know the overall
variance/covariance matrix of $Y_t$, which we'll call $\Phi$; the logic
is imilar to the univariate case, and I'll show it here for
completeness:

$\begin{array}{rl}
  \mathrm{Var} (Y_t) & = \mathrm{Var} (\rho^2 Y_{t - 1} + \epsilon_t)\\
  & = \rho^2 \mathrm{Var} (Y_{t - 1}) + \mathrm{Var}(\epsilon_t)\\
  \Phi & = \rho^2 \Phi + \Sigma_{\epsilon}\\
  \Phi - \rho^2 \Sigma & = \Sigma_{\epsilon}\\
  \Phi (I - \rho^2 I) & = \Sigma_{\epsilon}\\
  \Phi & = \Sigma_{\epsilon} (I - \rho^2
  I)^{- 1}
\end{array}$

If we pull out two sites at the same time point, their covariance is
$\mathrm{cov} (Y_{t, c_1}, Y_{t, c_2}) = \frac{\Sigma_{\epsilon, c_1, c_2}}{1- \rho^2}$, which looks very similar to the unidimensional AR(1) process
variance.

Now we derive the covariance between any two sites that are one time
point apart:

$\begin{array}{rl}
  \mathrm{cov} (y_{c_1, t}, y_{c_2, t + 1}) & = E [(y_{c_1, t} - E [y_{c_1, t}]) (y_{c_2, t} - E
[y_{c_2, t}])]\\
  & = E [y_{c_1, t} y_{c_2, t}]\\
  & = E [y_{c_1, t} [\rho y_{c_2, t} + \epsilon_{c_2, t + 1}]]\\
  & = \rho E [y_{c_1, t} y_{c_2, t}]\\
  & = \rho \mathrm{cov} (y_{c_1, t} y_{c_2, t})\\
  & = \rho \frac{\Sigma_{i, j}}{1 - \rho^2}\\
  & = \rho \frac{1}{1 - \rho^2} \Sigma_{i, j}\\
  & = \rho \frac{1}{1 - \rho^2} \alpha^2 \exp ~
  \left( - \frac{d (i, j)}{2 \lambda^2} \right)
\end{array}$

for sites more than one time point away from each other, we can iterate
the above result to get a general expression of the covariance between
any two points:

$ \mathrm{cov} (y_{c_1, t_1}, y_{c_2, t_2})  = \rho^{\mid t_1 - t_2 \mid} \frac{1}{1 - \rho^2} \alpha^2 \exp ~ \left( - \frac{d (i, j)}{2 \lambda^2} \right) $

if we reparameterize $\alpha$ to be the product of two parameters
$\alpha = \sigma_{\epsilon}^2 \alpha$, we get

$\begin{array}{rl}
  \mathrm{cov} (y_{c_1, t_1}, y_{c_2, t_2}) & = \rho^{\mid t_1 - t_2 \mid} \frac{\sigma_\epsilon^2} {1 - \rho^2} \alpha^2 \exp ~ \left( - \frac{d (i, j)}{2 \lambda^2} \right)\\
      & = K_1 (i, j) \times K_2 (i, j)
\end{array}$

which is the product of an AR(1) and squared exponential kernel function
as defined in the previous section. In practice we wouldn't want to
separate these parameters because both of them will not be identifiable
given observed data, but I separated them here to show how the
covariance structure is the product of two kernel functions.

Therefore, we can write this process in the form of a Gaussian process
with mean zero and covariance kernel given by the product of a temporal
and spatial kernel:

$\begin{array}{rl}
  Y \sim & \mathrm{MVN} (0, \Sigma)\\
  \Sigma_{i, j} = & K_1 (i, j) \times K_2 (i, j)
\end{array}$

The spatio-temporal processes defined as a set of conditional
distributions and as a Gaussian process are equivalent.

 

To summarize, AR processes can be written as a Gaussian process model,
which is useful because a temporal process can then be easily combined
with other sources of dependence. In general, we can build our models by
defining conditional distributions with a given mean and covariance, or
a joint distribution with mean zero where the model is fully defined by
a variance/covariance kernel function. In a future post I will look at
Bayesian parameter estimation in these models using Stan.


 
## Rasmussen and Williams (2006) treats AR(1) processes in their Appendix B. 

Have a slightly different solution ... need to verify previous method s

They start with:

$X_t = \Sigma^p_{k=1} a_k X_{t-k} + b_0 Z_t$ 

where $Z_t \sim N(0, 1)$ is IID with order-$p=1$; $a =$ rho; and $b_0 = $ sigma scales the $Z_t$ 

which can be expressed as a SDE:

$X^\prime(t) + a_0 X(t) = b_0Z(t),$ with $a_0 >0$ for stationarity

This gives a power spectrum for AR(1) (Rasmussen and Williams 2006, eq B.29):

$S(s) = \frac{b_0^2}{ (2 \pi s)^2+ a_0^2}$

or alternatively as (Rasmussen and Williams 2006: 214, eq. B.43):

$S(s) = \frac{b_0^2}{1-2a_1\cos(\omega_i s) + a_1^2}$


It has a Fourier transform: $k(t) = \frac{b^2_0}{2 a_0} e^{-a_0|t|} $



# D. Additional resources: 

https://luiarthur.github.io/TuringBnpBenchmarks/gp

https://github.com/JuliaGaussianProcesses/Stheno.jl

https://juliagaussianprocesses.github.io/KernelFunctions.jl/stable/examples/gaussian-process-priors/

https://storopoli.io/Bayesian-Julia/pages/11_multilevel_models/

https://gist.github.com/luiarthur/7a1dfa6a980d11a00862152a371a5cf3


# E. Examples:
  
```julia 

using Random
using Distributions
using Plots

n = 1000
x = rand(n)
xr = -6:0.1:6 

sigma = 2.0
rho = 0.75

using AbstractGPs
K1 = 1 * Matern32Kernel()
K = sigma^2 * Matern32Kernel()

# plot distribution
plot(  xr, [ K1(1,i) for i in xr] )
plot!( xr, [ K1(-1,i) for i in xr] )
plot!( xr, [ K(0,i) for i in xr] )

f = GP(0.0, K) # Define GP prior with zero mean and cov that is Matern-3/2 kernel
rand(f([1]))

fx = f(x)
mean(fx)
cov(fx)  # covariance  ==  kernelmatrix(sigma^2 * Matern32Kernel(), x)

ys = rand( fx )  # a sample from the GP 


## Visualize it

# Load required packages
using LinearAlgebra, KernelFunctions

using Plots, Plots.PlotMeasures
default(; lw=1.0, legendfontsize=8.0)

using Random: seed!
seed!(42); # reproducibility

num_inputs = 101
xlim = (-5, 5)
X = range(xlim...; length=num_inputs);

num_samples = 7
v = randn(num_inputs, num_samples);
 
function mvn_sample(K, v)
    L = cholesky(K + 1e-6 * I)
    f = L.L * v
    return f
end;

KM = kernelmatrix(K, X)  # 101 X 101 
f = mvn_sample(KM, v)  # 101 X 7 matrix  

plot(X, f; c="blue", title=raw"$f(x)$", ylim=(-3, 3), label=nothing)

COV = cov(KM) # cov matrix
f2 =  ( rand(MvNormal(COV + 1e-6 * I))) .* v  # 101 X 7 matrix


function visualize(k::Kernel, v)
    K = kernelmatrix(k, X)
    f = mvn_sample(K, v)

    p_kernel_2d = heatmap(
        X,
        X,
        K;
        yflip=true,
        colorbar=false,
        ylabel=string(nameof(typeof(k))),
        ylim=xlim,
        yticks=([xlim[1], 0, xlim[end]], ["\u22125", raw"$x'$", "5"]),
        vlim=(0, 1),
        title=raw"$k(x, x')$",
        aspect_ratio=:equal,
        left_margin=5mm,
    )

    p_kernel_cut = plot(
        X,
        k.(X, 0.0);
        title=string(raw"$k(x, x_\mathrm{ref})$"),
        label=raw"$x_\mathrm{ref}=0.0$",
        legend=:topleft,
        foreground_color_legend=nothing,
    )
    plot!(X, k.(X, 1.5); label=raw"$x_\mathrm{ref}=1.5$")

    p_samples = plot(X, f; c="blue", title=raw"$f(x)$", ylim=(-3, 3), label=nothing)

    return plot(
        p_kernel_2d,
        p_kernel_cut,
        p_samples;
        layout=(1, 3),
        xlabel=raw"$x$",
        xlim=xlim,
        xticks=collect(xlim),
    )
end;

plot(visualize(SqExponentialKernel(), v); size=(800, 210), bottommargin=5mm, topmargin=5mm)

kernels = [
    Matern12Kernel(),
    Matern32Kernel(),
    Matern52Kernel(),
    SqExponentialKernel(),
    WhiteKernel(),
    ConstantKernel(),
    LinearKernel(),
    compose(PeriodicKernel(), ScaleTransform(0.2)),
    NeuralNetworkKernel(),
    GibbsKernel(; lengthscale=x -> sum(exp ∘ sin, x)),
]

plot(
    [visualize(k, v) for k in kernels]...;
    layout=(length(kernels), 1),
    size=(800, 220 * length(kernels) + 100),
)


# AbstractGPs test:
# https://github.com/JuliaGaussianProcesses/AbstractGPs.jl/blob/master/test/ppl/turing.jl

using AbstractGPs
using Turing
import ForwardDiff

k = SqExponentialKernel()
y = randn(3)
X = randn(3, 1)
x = [rand(1) for _ in 1:3]

Turing.@model function GPRegression(y, X)
    # Priors.
    α ~ LogNormal(0.0, 0.1)
    ρ ~ LogNormal(0.0, 1.0)
    σ² ~ LogNormal(0.0, 1.0)

    # Realized covariance function
    kernel = α * (SqExponentialKernel() ∘ ScaleTransform(1 / ρ))   
    f = GP(kernel)

    # Sampling Distribution.
    return y ~ f(X, σ²)
end

# Test for matrices
m = GPRegression(y, RowVecs(X))
res = Turing.sample(m, Turing.NUTS(), 5) 


# Test for vectors of vector
m = GPRegression(y, x)



X = randn(3, 1)
x = [rand(1) for _ in 1:3]
y = rand.(Poisson.(exp.(randn(3))))

Turing.@model function latent_gp_regression(y, X)
    f = GP(Matern32Kernel())
    u ~ f(X)  
    λ = exp.(u)
    return y .~ Poisson.(λ)
end
m = latent_gp_regression(y, RowVecs(X))
res = Turing.sample(m, Turing.NUTS(), 5) 

# Test for vectors of vector
m = latent_gp_regression(y, x)
res = Turing.sample(m, Turing.NUTS(), 5)


```

AbstractGPs test:

copied from: https://github.com/JuliaGaussianProcesses/AbstractGPs.jl/blob/master/test/ppl/turing.jl



```julia

using AbstractGPs
using Turing
import ForwardDiff

k = SqExponentialKernel()
y = randn(3)
X = randn(3, 1)
x = [rand(1) for _ in 1:3]

Turing.@model function GPRegression(y, X)
    # Priors.
    α ~ LogNormal(0.0, 0.1)
    ρ ~ LogNormal(0.0, 1.0)
    σ² ~ LogNormal(0.0, 1.0)

    # Realized covariance function
    kernel = α * (SqExponentialKernel() ∘ ScaleTransform(1 / ρ))   
    f = GP(kernel)

    # Sampling Distribution.
    return y ~ f(X, σ²)
end

# Test for matrices
m = GPRegression(y, RowVecs(X))
res = Turing.sample(m, Turing.NUTS(), 5) 


# Test for vectors of vector
m = GPRegression(y, x)
 

X = randn(3, 1)
x = [rand(1) for _ in 1:3]
y = rand.(Poisson.(exp.(randn(3))))

Turing.@model function latent_gp_regression(y, X)
    f = GP(Matern32Kernel())
    u ~ f(X)
    λ = exp.(u)
    return y .~ Poisson.(λ)
end
m = latent_gp_regression(y, RowVecs(X))
res = Turing.sample(m, Turing.NUTS(), 5) 

# Test for vectors of vector
m = latent_gp_regression(y, x)
res = Turing.sample(m, Turing.NUTS(), 5)

```

 

```julia 

using Distributions, Plots, Random

function ar1( ; x_init, n, rho, sigma )
  x = zeros(n)
  x[1] = x_init
  for i in 2:n 
    x[i] = rho * x[i-1] + rand(Normal(0, sigma))
  end
  return x
end

n = 1000
x = rand(n)
xr = -6:0.1:6 

sigma = 2.0
rho = 0.75

u = ar1( x_init=0, n=n, rho=rho, sigma=sigma )
plot(1:n, u )
mean(u) # 0.043842140590446714
std(u) # 1.1339969369576197

histogram(u)
 
# increments
du = u[1:(n-1)] - u[2:n]
mean(du) # -0.0002578562327382s6734
std(du) # 1.1339969369576197

histogram(du)
 

```


```julia  

using AbstractGPs

cov_ar1 = (rho) /(1-rho^2)
Kar1 = (1/sigma^2)  * ( ConstantKernel(c=cov_ar1) + WhiteKernel() )  
f = GP(0.0, Kar1) # Define GP prior with zero mean and cov that is Matern-3/2 kernel
 
fx= f([1,2])

mean(fx)
cov(fx)  # covariance  ==  kernelmatrix(Matern32Kernel(), x)

rw =  vec( rand( fx, 1000 ))  # a sample from the GP 
plot(rw)
histogram(rw)
mean(rw) #  
std(rw)  # 

# increments
du = rw[1:(n-1)] - rw[2:n]
mean(du) # -0.00025785623273826734
std(du) # 1.1339969369576197

histogram(du)
 
```




```julia 

using Distributions, Plots, Random
using AbstractGPs

n = 5000
x = rand(n)
xr = -6:0.1:6 

sigma = 2.0
rho = 0.75

K1 = 1 * Matern32Kernel()
K = sigma^2 * Matern32Kernel()

# plot distribution
plot(  xr, [ K1(1,i) for i in xr] )
plot!( xr, [ K1(-1,i) for i in xr] )
plot!( xr, [ K(0,i) for i in xr] )

f = GP(0.0, K) # Define GP prior with zero mean and cov that is Matern-3/2 kernel
rand(f([1]))

fx = f(x)
mean(fx)
cov(fx)  # covariance  ==  kernelmatrix(sigma^2 * Matern32Kernel(), x)

ys = rand( fx )  # a sample from the GP 

## simple AR1 process
function ar1( ; x_init, n, rho, sigma )
  x = zeros(n)
  r = rand(Normal(0, sigma), n-1)
  x[1] = x_init
  for i in 2:n 
    x[i] = rho * x[i-1] + r[i-1]
  end
  return x
end

u = ar1( x_init=0, n=n, rho=rho, sigma=sigma )
plot(1:n, u )
mean(u) # 0
std(u) # 3

histogram(u)
  
   
# GP   
cov_ar1 =  ((sigma^2 * rho) /(1-rho^2))
Kar1 = (1-rho^2)/rho * EyeKernel() * ConstantKernel(c=cov_ar1)  
f = GP(0.0, Kar1) # Define GP prior with zero mean and cov that is Matern-3/2 kernel
 
fx= f([1,2])

mean(fx)
cov(fx)  

rw =  vec( rand( fx, n ))  # a sample from the GP 
plot(rw)
histogram(rw)
mean(rw) #  
std(rw)  # 2.6142129622102956
 



using AbstractGPs, Plots

# Generate toy synthetic data.
x = rand(10);  y = sin.(x)

# Define GP prior with Matern-3/2 kernel
f = GP(Matern32Kernel())

# Finite projection of `f` at inputs `x`.
# Added Gaussian noise with variance 0.001.
fx = f(x, 0.001)

# Log marginal probability of `y` under `f` at `x`.
# Quantity typically maximised to train hyperparameters.
logpdf(fx, y)

# Exact posterior given `y`. This is another GP.
p_fx = posterior(fx, y)

# Log marginal posterior predictive probability.
logpdf(p_fx(x), y)

# Plot posterior.
scatter(x, y; label="Data")
plot!(-0.5:0.001:1.5, p_fx; label="Posterior")

x2 = rand(10).*2
ys = similar(x2)

rand!(p_fx(x2), ys)  # sample and store into ys
plot!( x, ys, seriestype=:scatter, legend=:none)

x = rand(10)
scatter(x, y; label="Data")
plot!(-0.5:0.001:1.5, p_fx; label="Posterior")
ys = rand(p_fx(x))  # sample and store into ys
plot!( x, ys, seriestype=:scatter, legend=:none)


```

  

  

# MISC

AbstractGPs: example 
https://github.com/JuliaGaussianProcesses/AbstractGPs.jl/blob/master/examples/0-intro-1d/script.jl


```julia

using AbstractGPs
using Distributions
using FillArrays
using StatsFuns

using Plots
default(; legend=:outertopright, size=(700, 400))

using Random
Random.seed!(42)  # setting the seed for reproducibility of this notebook
#md nothing #hide

# Load toy regression
# [dataset](https://github.com/GPflow/GPflow/blob/7705cee6723f78066981f27954130daaede55dfc/doc/sphinx/notebooks/basics/data/regression_1D.csv)
# taken from GPflow examples.

x = [
    0.8658165855998895,
    0.6661700880180962,
    0.8049218148148531,
    0.7714303440386239,
    0.14790478354654835,
    0.8666105548197428,
    0.007044577166530286,
    0.026331737288148638,
    0.17188596617099916,
    0.8897812990554013,
    0.24323574561119998,
    0.028590102134105955,
]
y = [
    1.5255314337144372,
    3.6434202968230003,
    3.010885733911661,
    3.774442382979625,
    3.3687639483798324,
    1.5506452040608503,
    3.790447985799683,
    3.8689707574953,
    3.4933565751758713,
    1.4284538820635841,
    3.8715350915692364,
    3.7045949061144983,
]
scatter(x, y; xlabel="x", ylabel="y", legend=false)

# We split the observations into train and test data.

x_train = x[1:8]
y_train = y[1:8]
x_test = x[9:end]
y_test = y[9:end]
#md nothing #hide

# We instantiate a Gaussian process with a Matern kernel. The kernel has
# fixed variance and length scale parameters of default value 1.

f = GP(Matern52Kernel())
#md nothing #hide

# We create a finite dimensional projection at the inputs of the training dataset
# observed under Gaussian noise with variance $\sigma^2 = 0.1$, and compute the
# log-likelihood of the outputs of the training dataset.

fx = f(x_train, 0.1)
logpdf(fx, y_train)

# We compute the posterior Gaussian process given the training data, and calculate the
# log-likelihood of the test dataset.

p_fx = posterior(fx, y_train)
logpdf(p_fx(x_test), y_test)

# We plot the posterior Gaussian process (its mean and a ribbon of 2 standard deviations
# around it) on a grid along with the observations.

scatter(
    x_train,
    y_train;
    xlim=(0, 1),
    xlabel="x",
    ylabel="y",
    title="posterior (default parameters)",
    label="Train Data",
)
scatter!(x_test, y_test; label="Test Data")
plot!(0:0.001:1, p_fx; label=false, ribbon_scale=2)



# ## Markov Chain Monte Carlo
#
# Previously we computed the log likelihood of the untuned kernel parameters of the GP.
# We now also perform approximate inference over said kernel parameters using different
# Markov chain Monte Carlo (MCMC) methods. I.e., we approximate the posterior distribution
# of the kernel parameters with samples from a Markov chain.
#
# We define a function which returns the log-likelihood of the data for different variance
# and inverse lengthscale parameters of the Matern kernel. We ensure that these parameters are
# positive with the softplus function
# ```math
# f(x) = \log (1 + \exp x).
# ```

function gp_loglikelihood(x, y)
    function loglikelihood(params)
        kernel =
            softplus(params[1]) * (Matern52Kernel() ∘ ScaleTransform(softplus(params[2])))
        f = GP(kernel)
        fx = f(x, 0.1)
        return logpdf(fx, y)
    end
    return loglikelihood
end

const loglik_train = gp_loglikelihood(x_train, y_train)
#md nothing #hide

# We define a Gaussian prior for the joint distribution of the two transformed kernel
# parameters. We assume that both parameters are independent with mean 0 and variance 1.

logprior(params) = logpdf(MvNormal(Eye(2)), params)
#md nothing #hide

# ### Hamiltonian Monte Carlo
#
# We start with a Hamiltonian Monte Carlo (HMC) sampler. More precisely, we use the
# [No-U-Turn sampler (NUTS)](http://www.jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf),
# which is provided by the Julia packages
# [AdvancedHMC.jl](https://github.com/TuringLang/AdvancedHMC.jl/) and
# [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl/).
#
# #### AdvancedHMC
#
# We start with performing inference with AdvancedHMC.

using AdvancedHMC
using ForwardDiff

# Set the number of samples to draw and warmup iterations.

n_samples = 2_000
n_adapts = 1_000
#md nothing #hide

# AdvancedHMC and DynamicHMC below require us to implement the LogDensityProblems interface for the log joint probability.

using LogDensityProblems

struct LogJointTrain end

## Log joint density
LogDensityProblems.logdensity(::LogJointTrain, θ) = loglik_train(θ) + logprior(θ)

## The parameter space is two-dimensional
LogDensityProblems.dimension(::LogJointTrain) = 2

## `LogJointTrain` does not allow to evaluate derivatives of the log density function
function LogDensityProblems.capabilities(::Type{LogJointTrain})
    return LogDensityProblems.LogDensityOrder{0}()
end
#md nothing #hide

# We use [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to compute the derivatives of the log joint density with automatic differentiation.

using LogDensityProblemsAD

const logjoint_train = ADgradient(Val(:ForwardDiff), LogJointTrain())
#md nothing #hide

# We define an Hamiltonian system of the log joint probability.

metric = DiagEuclideanMetric(2)
hamiltonian = Hamiltonian(metric, logjoint_train)
#md nothing #hide

# Define a leapfrog solver, with initial step size chosen heuristically.

initial_params = rand(2)
initial_ϵ = find_good_stepsize(hamiltonian, initial_params)
integrator = Leapfrog(initial_ϵ)
#md nothing #hide

# Define an HMC sampler, with the following components:
# - multinomial sampling scheme,
# - generalised No-U-Turn criteria, and
# - windowed adaption for step-size and diagonal mass matrix

proposal = NUTS{MultinomialTS,GeneralisedNoUTurn}(integrator)
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
#md nothing #hide

# We draw samples from the posterior distribution of kernel parameters. These samples
# are in the unconstrained space $\mathbb{R}^2$.

samples, _ = sample(
    hamiltonian, proposal, initial_params, n_samples, adaptor, n_adapts; progress=false
)
#md nothing #hide

# We transform the samples back to the constrained space and compute the mean of both
# parameters:

samples_constrained = [map(softplus, p) for p in samples]
mean_samples = mean(samples_constrained)

# We plot a histogram of the samples for the two parameters.
# The vertical line in each graph indicates the mean of the samples.

histogram(
    reduce(hcat, samples_constrained)';
    xlabel="sample",
    ylabel="counts",
    layout=2,
    title=["variance" "inverse length scale"],
    legend=false,
)
vline!(mean_samples'; linewidth=2)

# We approximate the log-likelihood of the test data using the posterior Gaussian processes
# for kernels with the sampled kernel parameters. We can observe that there is a significant
# improvement over the log-likelihood of the test data with respect to the posterior
# Gaussian process with default kernel parameters of value 1.

function gp_posterior(x, y, p)
    kernel = softplus(p[1]) * (Matern52Kernel() ∘ ScaleTransform(softplus(p[2])))
    f = GP(kernel)
    return posterior(f(x, 0.1), y)
end

mean(logpdf(gp_posterior(x_train, y_train, p)(x_test), y_test) for p in samples)

# We sample 5 functions from each posterior GP given by the final 100 samples of kernel
# parameters.

plt = plot(; xlim=(0, 1), xlabel="x", ylabel="y", title="posterior (AdvancedHMC)")
for (i, p) in enumerate(samples[(end - 100):end])
    sampleplot!(
        plt,
        0:0.02:1,
        gp_posterior(x_train, y_train, p);
        samples=5,
        seriescolor="red",
        label=(i == 1 ? "samples" : nothing),
    )
end
scatter!(plt, x_train, y_train; label="Train Data", markercolor=1)
scatter!(plt, x_test, y_test; label="Test Data", markercolor=2)
plt

# #### DynamicHMC
#
# We repeat the inference with DynamicHMC.

using DynamicHMC

samples =
    mcmc_with_warmup(
        Random.GLOBAL_RNG, logjoint_train, n_samples; reporter=NoProgressReport()
    ).posterior_matrix
#md nothing #hide

# We transform the samples back to the constrained space and compute the mean of both
# parameters:

samples_constrained = map(softplus, samples)
mean_samples = vec(mean(samples_constrained; dims=2))

# We plot a histogram of the samples for the two parameters.
# The vertical line in each graph indicates the mean of the samples.

histogram(
    samples_constrained';
    xlabel="sample",
    ylabel="counts",
    layout=2,
    title=["variance" "inverse length scale"],
    legend=false,
)
vline!(mean_samples'; linewidth=2)

# Again we can observe that there is a significant improvement over the log-likelihood
# of the test data with respect to the posterior Gaussian process with default kernel
# parameters.

mean(logpdf(gp_posterior(x_train, y_train, p)(x_test), y_test) for p in eachcol(samples))

# We sample a function from the posterior GP for the final 100 samples of kernel
# parameters.

plt = plot(; xlim=(0, 1), xlabel="x", ylabel="y", title="posterior (DynamicHMC)")
scatter!(plt, x_train, y_train; label="Train Data")
scatter!(plt, x_test, y_test; label="Test Data")
for i in (n_samples - 100):n_samples
    p = @view samples[:, i]
    sampleplot!(plt, 0:0.02:1, gp_posterior(x_train, y_train, p); seriescolor="red")
end
plt

# ### Elliptical slice sampling
#
# Instead of HMC, we use
# [elliptical slice sampling](http://proceedings.mlr.press/v9/murray10a/murray10a.pdf)
# which is provided by the Julia package
# [EllipticalSliceSampling.jl](https://github.com/TuringLang/EllipticalSliceSampling.jl/).

using EllipticalSliceSampling

# We draw 2000 samples from the posterior distribution of kernel parameters.

samples = sample(ESSModel(
    MvNormal(Eye(2)), # Gaussian prior
    loglik_train,
), ESS(), n_samples; progress=false)
#md nothing #hide

# We transform the samples back to the constrained space and compute the mean of both
# parameters:

samples_constrained = [map(softplus, p) for p in samples]
mean_samples = mean(samples_constrained)

# We plot a histogram of the samples for the two parameters.
# The vertical line in each graph indicates the mean of the samples.

histogram(
    reduce(hcat, samples_constrained)';
    xlabel="sample",
    ylabel="counts",
    layout=2,
    title=["variance" "inverse length scale"],
)
vline!(mean_samples'; layout=2, labels="mean")

# Again we can observe that there is a significant improvement over the log-likelihood
# of the test data with respect to the posterior Gaussian process with default kernel
# parameters.

mean(logpdf(gp_posterior(x_train, y_train, p)(x_test), y_test) for p in samples)

# We sample a function from the posterior GP for the final 100 samples of kernel
# parameters.

plt = plot(;
    xlim=(0, 1), xlabel="x", ylabel="y", title="posterior (EllipticalSliceSampling)"
)
scatter!(plt, x_train, y_train; label="Train Data")
scatter!(plt, x_test, y_test; label="Test Data")
for p in samples[(end - 100):end]
    sampleplot!(plt, 0:0.02:1, gp_posterior(x_train, y_train, p); seriescolor="red")
end
plt

# ## Variational Inference
#
# Sanity check for the Evidence Lower BOund (ELBO) implemented according to
# M. K. Titsias's _Variational learning of inducing variables in sparse Gaussian processes_.

elbo(VFE(f(rand(5))), fx, y_train)

# We use the LBFGS algorithm to maximize the given ELBO. It is provided by the Julia
# package [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl).

using Optim

# We define a function which returns the negative ELBO for different variance and inverse
# lengthscale parameters of the Matern kernel and different pseudo-points. We ensure that
# the kernel parameters are positive with the softplus function
# ```math
# f(x) = \log (1 + \exp x),
# ```
# and that the pseudo-points are in the unit interval $[0,1]$ with the logistic function
# ```math
# f(x) = \frac{1}{1 + \exp{(-x)}}.
# ```

jitter = 1e-6  # "observing" the latent process with some (small) amount of jitter improves numerical stability

function objective_function(x, y)
    function negative_elbo(params)
        kernel =
            softplus(params[1]) * (Matern52Kernel() ∘ ScaleTransform(softplus(params[2])))
        f = GP(kernel)
        fx = f(x, 0.1)
        z = logistic.(params[3:end])
        approx = VFE(f(z, jitter))
        return -elbo(approx, fx, y)
    end
    return negative_elbo
end
#md nothing #hide

# We randomly initialize the kernel parameters and 5 pseudo points, and minimize the
# negative ELBO with the LBFGS algorithm and obtain the following optimal parameters:

x0 = rand(7)
opt = optimize(objective_function(x_train, y_train), x0, LBFGS())

#-

opt.minimizer

# The optimized value of the variance is

softplus(opt.minimizer[1])

# and of the inverse lengthscale is

softplus(opt.minimizer[2])

# We compute the log-likelihood of the test data for the resulting approximate
# posterior. We can observe that there is a significant improvement over the
# log-likelihood with the default kernel parameters of value 1.

opt_kernel =
    softplus(opt.minimizer[1]) *
    (Matern52Kernel() ∘ ScaleTransform(softplus(opt.minimizer[2])))
opt_f = GP(opt_kernel)
opt_fx = opt_f(x_train, 0.1)
ap = posterior(VFE(opt_f(logistic.(opt.minimizer[3:end]), jitter)), opt_fx, y_train)
logpdf(ap(x_test), y_test)

# We visualize the approximate posterior with optimized parameters.

scatter(
    x_train,
    y_train;
    xlim=(0, 1),
    xlabel="x",
    ylabel="y",
    title="posterior (VI with sparse grid)",
    label="Train Data",
)
scatter!(x_test, y_test; label="Test Data")
plot!(0:0.001:1, ap; label=false, ribbon_scale=2)
vline!(logistic.(opt.minimizer[3:end]); label="Pseudo-points")

# ## Exact Gaussian Process Inference
#
# Here we use Type-II MLE to train the hyperparameters of the Gaussian process.
# This means that our loss function is the negative log marginal likelihood.

# We re-calculate the log-likelihood of the test dataset with the
# default kernel parameters of value 1 for the sake of comparison.

logpdf(p_fx(x_test), y_test)

# We define a function which returns the negative log marginal
# likelihood for different variance and inverse lengthscale parameters
# of the Matern kernel and different pseudo-points. We ensure that the
# kernel parameters are positive with the softplus function
# ``f(x) = \log (1 + \exp x)``.

function loss_function(x, y)
    function negativelogmarginallikelihood(params)
        kernel =
            softplus(params[1]) * (Matern52Kernel() ∘ ScaleTransform(softplus(params[2])))
        f = GP(kernel)
        fx = f(x, 0.1)
        return -logpdf(fx, y)
    end
    return negativelogmarginallikelihood
end

#md nothing #hide

# We randomly initialize the kernel parameters, and minimize the
# negative log marginal likelihood with the LBFGS algorithm
# and obtain the following optimal parameters:

θ = randn(2)
opt = Optim.optimize(loss_function(x_train, y_train), θ, LBFGS())

#-

opt.minimizer

# The optimized value of the variance is

softplus(opt.minimizer[1])

# and of the inverse lengthscale is

softplus(opt.minimizer[2])

# We compute the log-likelihood of the test data for the resulting optimized
# posterior. We can observe that there is a significant improvement over the
# log-likelihood with the default kernel parameters of value 1.

opt_kernel =
    softplus(opt.minimizer[1]) *
    (Matern52Kernel() ∘ ScaleTransform(softplus(opt.minimizer[2])))

opt_f = GP(opt_kernel)
opt_fx = opt_f(x_train, 0.1)
opt_p_fx = posterior(opt_fx, y_train)
logpdf(opt_p_fx(x_test), y_test)

# We visualize the posterior with optimized parameters.

scatter(
    x_train,
    y_train;
    xlim=(0, 1),
    xlabel="x",
    ylabel="y",
    title="posterior (optimized parameters)",
    label="Train Data",
)
scatter!(x_test, y_test; label="Test Data")
plot!(0:0.001:1, opt_p_fx; label=false, ribbon_scale=2)





# Import Libraries
using Turing
using Turing: Variational
using Distributions
using Distances
using PyPlot
using StatsFuns
import Random
using Flux
import LinearAlgebra

# Squared-exponential covariance function
sqexp_cov_fn(D, phi, eps=1e-3) = exp.(-D^2 / phi) + LinearAlgebra.I * eps

# Exponential covariance function
exp_cov_fn(D, phi) = exp.(-D / phi)

@model function GP(y, X, m=0, s=1, cov_fn=exp_cov_fn)
    # Dimensions of predictors .
    N, P = size(X)
    
    # Distance matrix.
    D = pairwise(Distances.Euclidean(), X, dims=1)
    
    # Priors.
    mu ~ Normal(m, s)
    sig2 ~ LogNormal(0, 1)
    phi ~ LogNormal(0, 1)
    
    # Realized covariance function
    K = cov_fn(D, phi)
    
    # Sampling Distribution.
    # The latent variables have been marginalized out here,
    # so there's only one level.
    y ~ MvNormal(mu * ones(N), K + sig2 * LinearAlgebra.I(N))

    # # Prior for linear model coefficients
    #  beta ~ Normal(zeros(J), s_beta)  # you may need to choose s_beta carefully.
#   # Sampling Distribution.
#   y ~ MvNormal(mu * ones(N) + Z * beta, K + sig2 * LinearAlgebra.I(N))


end

# Generate some non-linear data
N = 150
X = randn(N, 1) * 3
y = sin.(vec(X)) + randn(N) * 0.1
scatter(vec(X), y)

# Fit via ADVI. You can also use HMC.
Random.seed!(0)

m = GP(y, X)


q0 = Variational.meanfield(m)  # initialize variational distribution (optional)
advi = ADVI(1, 2000)  # num_elbo_samples, max_iters
@time q = vi(m, advi, q0, optimizer=Flux.ADAM(1e-1));


nsamples = 1000  # number of MCMC samples
nadapt = 1000  # number of iterations to adapt tuning parameters in NUTS
iterations = nsamples + nadapt
target_accept_ratio = 0.8

chain = sample(m, NUTS(nadapt, target_accept_ratio, max_depth=10), iterations);




---

using CSV
using DataFrames
using Turing
using Distributions
using LinearAlgebra
using KernelFunctions

df = DataFrame(
    subject = repeat(1:5, inner=3),
    obese = repeat(rand(Bool, 5), inner=3),
    timepoint = [1,2,3,1,3,4,1,2,5,1,4,5,1,3,5],
    bug = rand(Beta(0.9, 5), 15),
    nutrient = rand(Beta(0.9,5), 15)
)


# Define a parameterized kernel.
sqexpkernel(alpha::Real, rho::Real) = 
    alpha^2 * SqExponentialKernel() ∘ ScaleTransform(1/rho*sqrt(2))


@model function myGP(y, Z, X, jitter=1e-6)
    # Dimensions of GP predictors. Note that if X is a single column 
    # it needs to have dimension Nx1 (a matrix), not an array of 
    # length N, so that it can be processed properly by the distance function.
    N, P = size(X)
    
    # Dimensions of linear model predictors
    J = size(Z, 2)  # Z should be N x J
    
    # Priors.
    mu ~ Normal(0, 1)
    sig2 ~ LogNormal(0, 1)
    alpha ~ LogNormal(0, 0.1)
    rho ~ LogNormal(0, 1)

    # Prior for linear model coefficients
    beta ~ filldist(Normal(0, 1), J)
    
    # GP Covariance matrix
    kernel = sqexpkernel(alpha, rho)  # covariance function
    K = kernelmatrix(kernel, X, obsdim=1)  # cov matrix
    K += LinearAlgebra.I * (sig2 + jitter)

    # Sampling Distribution.
    y ~ MvNormal(mu .+ Z * beta, K)
end

y = df.bug
Z = Matrix(df[!,[:subject, :obese, :nutrient]])
X = Matrix{Float64}(df[!,[:timepoint]])

mygp = myGP(y, Z, X)

chain = sample(mygp, HMC(0.01, 100), 200)


 
---


# k = SEKernel()
# jitter = 1e-6
# # adds jitter to diagonal of kernel matrix 
# k +=  jitter * WhiteKernel()  
# # `EyeKernel` is an alias for `WhiteKernel`
# # So, `k += jitter * EyeKernel()` also works.

using AbstractGPs, KernelFunctions, Distributions


# Define model.
@model GPRegression(y, X, jitter=1e-6) = begin
    # Priors.
    alpha ~ LogNormal(0.0, 0.1)
    rho ~ LogNormal(0.0, 1.0)
    sigma ~ LogNormal(0.0, 1.0)
     
    # GP (implicit zero-mean, cov=sqexp).
    gpreg = GP(sqexpkernel(alpha, rho))
    
    # Sampling Distribution
    y ~ gpreg(X, sigma^2 + jitter)  # add jitter for numerical stability.
end;

y = df.bug
Z = Matrix(df[!,[:subject, :obese, :nutrient]])
X = Matrix{Float64}(df[!,[:timepoint]])

mygpreg = GPRegression(y, vec(X))

chain = sample(mygpreg, HMC(0.01, 100), 200)
 


```


## SPDE via RandomFields (in R)

```r
install.packages("RandomFieldsUtils")

library(rstan)
library(INLA)
library(RandomFields)

dat <- data.frame(x = runif(100),

 y = runif(100))

spat_field <- raster::raster(RFsimulate(RMmatern(nu=1, var = 1, scale = 0.1),
        x = seq(0, 1, length.out = 100),
        y = seq(0, 1, length.out = 100),
        spConform = FALSE))
dat$resp <- rnorm(100, mean = 1 + raster::extract(spat_field, dat), sd = 1)
bnd <- inla.mesh.segment(matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE))

# a coarse mesh
mesh <- inla.mesh.2d(max.edge = 0.2, offset = 0, boundary = bnd)

# derive the FEM matrices
fem <- inla.mesh.fem(mesh)
# put the matrices in one object
G <- array(c(as(fem$c1, "matrix"),
             as(fem$g1, "matrix"),
             as(fem$g2, "matrix")),
          dim = c(mesh$n, mesh$n, 3))
G <- aperm(G, c(3, 1, 2))

# derive the projection matrix
A <- inla.spde.make.A(mesh, loc = as.matrix(dat[,c("x", "y")]))
A <- as(A, "matrix")

y <- dat$resp

# the model parameters
mu <- normal(0, 1)
kappa <- lognormal(0, 1) # variance of spatial effect
tau <- gamma(0.1, 0.1) # scale of spatial effect
sigma <- lognormal(0, 1)
z <- normal(0, 1, dim = mesh$n) # std-normal variable for non-centered parametrization

# the precision matrix of the spatial effect
S <- (tau ** 2) * (G[1,,] * kappa ** 4 + G[2,,] * 2 * kappa ** 2 + G[3,,])
# drawing the spatial effect coefficient
beta <- backsolve(chol(S), z)

# the linear predictor
linpred <- mu + A %*% beta
distribution(y) <- normal(linpred, sigma)

# fitting the model
m_g <- model(mu, sigma, tau, kappa)
m_draw <- mcmc(m_g, one_by_one = TRUE)
bayesplot::mcmc_trace(m_draw)
coda::gelman.diag(m_draw)

# NOTE: lost model code?

```
