
function example_data(; N=1000)
    # make random data for analysis
    # NOTE: utility in terms of creating model matrix using schemas, etc
    
    xvar = vec(randn(N)*3.0  )
    
    df = DataFrame(
        xvar = xvar,
        covar1 = string.(rand( [1, 2, 3, 4], N)  ),  # factorial
        covar2 = vec(randn(N))
    )
    
    # group means
    cov2 = replace(df.covar1, "1"=>1.0, "2" =>1.25, "3"=>2, "4"=>1.5) 
    yvar = sin.(vec(xvar)) + df.covar2 .* 0.1  .+ 0.5 .* rand.(Normal.(cov2, 0.1))
    schm = StatsModels.schema(df, Dict(:covar1 => EffectsCoding()) )
    dm = StatsModels.apply_schema(StatsModels.term.((:xvar, :covar1, :covar2)), schm ) 
    dmat = StatsModels.modelcols(StatsModels.MatrixTerm( dm ), df)
    dmatnm = StatsModels.coefnames( dm )   # coef names of design matrix 
    tnm = StatsModels.termnames( dm)  # coef names of model data
    
    # alternative access:
    # fm = StatsModels.@formula( yvar ~ 1 + xvar + covar1 + covar2)
    # resp = response(fm, df)  # response
    # cols = modelcols(z, df)
    # o = reduce(hcat, cols)
    
    return yvar, df, dmat, dmatnm, tnm
end

function example_nonlinear_data(Xl = -7:0.1:7, ) 
    
    # function that describes latent process
    Y(x) = (x + 4) * (x + 1) * (x - 1) * (x - 3)  

    # Latent data ("truth")
    # Xl = -7:0.1:7
    yl = Y.(Xl)

    # "Observations" with noise
    Xo = -5:0.5:5   
    yo = Y.(X) + rand(Uniform(-20, 20), length(Xo))
    
    return Xl, yl, Xp, Xo, yo
end


function linear_regression(X, y, Xstar)
    beta = (X' * X) \ (X' * y)
    return Xstar * beta
end;

function featurize_poly(Xin, degree=1)
    # add higher order polynomials
    return repeat(Xin, 1, degree + 1) .^ (0:degree)'
end
 

function ridge_regression(X, y, Xstar, lambda)
    beta = (X' * X + lambda * I) \ (X' * y)
    return Xstar * beta
end


function kernel_ridge_regression(k, X, y, Xstar, lambda)
    K = kernelmatrix(k, X)
    kstar = kernelmatrix(k, Xstar, X)
    return kstar * ((K + lambda * I) \ y)
end;
 

function gp_kernel_ridge_regression_cholesky( v, s, r, Xo, Xp, Yo )
    # Gaussian Process Kernel Ridge Regression
    k = ( v * SqExponentialKernel() ) ∘ ScaleTransform(s)
    ko = kernelmatrix(k, Xo)
    kp = kernelmatrix(k, Xp, Xo)
    L = cholesky(ko + r * I)
    Yp = kp * ( L.U \ (L.L \ Yo) )  # mean process  
    return Yp 
end


function gp_kernel_ridge_regression( v, s, r, Xo, Xp, Yo )
    # need to check algebra
    # Gaussian Process Kernel Ridge Regression
    # k = ( v * SqExponentialKernel() ) ∘ ScaleTransform(s)
    # ko = kernelmatrix(k, Xo) + r * I  
    # kp = kernelmatrix(k, Xp, Xo) 
    # Yp = rand( MvNormal(  ko\Yo, ko ) ) # faster
    return Yp
end


function quantiles(X, q; dims, drop=false)
    Q = mapslices(x -> quantile(x, q), X, dims=dims)
    out = drop ? dropdims(Q, dims=dims) : Q
    return out
  end
  
       
  # Squared-exponential covariance function
  sqexp_cov_fn(D, phi, eps=1e-3) = exp.(-D^2 / phi) + LinearAlgebra.I * eps
  
  # Exponential covariance function
  exp_cov_fn(D, phi) = exp.(-D / phi)
  
  @model function GP(Y, X, m=0, s=1, cov_fn=exp_cov_fn)
  
      # Dimensions of predictors 
      N, P = size(X)
      
      # Distance matrix
      D = pairwise(Distances.Euclidean(), X, dims=1)
      
      # Priors
      mu ~ Normal(m, s)
      sig2 ~ LogNormal(0, 1)
      phi ~ LogNormal(0, 1)
      
      # Realized covariance function
      K = cov_fn(D, phi)
      
      # Sampling Distribution
      # The latent variables have been marginalized out here,
      # so there's only one level.
      Y ~ MvNormal(mu * ones(N), K + sig2 * LinearAlgebra.I(N))
  end



# This funciton returns a function for predicting at new points given parameter values.
function make_gp_predict_fn(Xnew, Y, X, cov_fn)
    N = size(X, 1)
    M = size(Xnew, 1)
    Q = N + M
    Z = [Xnew; X]
    D = pairwise(Euclidean(), Z, dims=1)
    
    return (mu, sig2, phi) -> let
        K = cov_fn(D, phi)
        Koo_inv = inv(K[(M+1):end, (M+1):end])
        Knn = K[1:M, 1:M]
        Kno = K[1:M, (M+1):end]
        C = Kno * Koo_inv
        m = C * (Y .- mu) .+ mu
        S = Matrix(LinearAlgebra.Hermitian(Knn - C * Kno'))
        mvn = MvNormal(m, S + sig2 * LinearAlgebra.I)
        rand(mvn)
    end
end


function make_extractor_avdi(m, q, nsamples=1000)
    # To extract parameters from ADVI model.
    qsamples = rand(q, nsamples)
    _, sym2range = Bijectors.bijector(m, Val(true));
    return sym -> qsamples[collect(sym2range[sym][1]), :]
end
