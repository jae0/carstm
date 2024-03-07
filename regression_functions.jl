
function example_data(; N=1000,  cov2lev = ("1"=>1, "2"=>1.25, "3"=>2, "4"=>1.5), alpha=0.1 )
    # make random data for analysis
    # NOTE: utility in terms of creating model matrix using schemas, etc
    
    xvar = vec(randn(N)*3.0  )
    
    df = DataFrame(
        xvar = xvar,
        covar1 = string.(rand( [1, 2, 3, 4], N)  ),  # factorial
        covar2 = vec(randn(N)),
        covar3 = vec(trunc.( Int, randn(N)*3 ) )
    )

    cov2 = replace(df.covar1, cov2lev[1], cov2lev[2], cov2lev[3], cov2lev[4] ) 
    df.yvar = sin.(vec(xvar)) + df.covar2 .* alpha .+ rand.(Normal.(cov2, 0.1))
    schm = StatsModels.schema(df, Dict(:covar1 => EffectsCoding()) )
    dm = StatsModels.apply_schema(StatsModels.term.((:xvar, :covar1, :covar2)), schm ) 
    modelcols = StatsModels.modelcols(StatsModels.MatrixTerm( dm ), df)
    coefnames = StatsModels.coefnames( dm )   # coef names of design matrix 
    termnames = StatsModels.termnames( dm)  # coef names of model data
    
    # alternative access:
    # fm = StatsModels.@formula( yvar ~ 1 + xvar + covar1 + covar2)
    # resp = response(fm, df)  # response
    # cols = modelcols(z, df)
    # o = reduce(hcat, cols)
    
    return df, modelcols, coefnames, termnames, cov2lev
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


function summarize_samples(S)
  smean = mean(S, dims=1)
  slb = quantiles(S, 0.025, dims=1)
  sub = quantiles(S, 0.975, dims=1)
  return smean, slb, sub
end
    
# Squared-exponential covariance function
sqexp_cov_fn(D, phi, eps=1e-3) = exp.(-D^2 / phi) + LinearAlgebra.I * eps

# Exponential covariance function
exp_cov_fn(D, phi) = exp.(-D / phi)


Turing.@model function gaussian_process_basic(; Y, D, cov_fn=exp_cov_fn, nY=length(Y) )
    mu ~ Normal(0.0, 1.0); # mean process
    sig2 ~ LogNormal(0, 1) # "nugget" variance
    phi ~ LogNormal(0, 1) # phi is ~ lengthscale along Xstar (range parameter)
    # sigma = cov_fn(D, phi) + sig2 * LinearAlgebra.I(nY) # Realized covariance function + nugget variance

    vcov = cov_fn(D, phi) + sig2 .* LinearAlgebra.I(nY )
    Y ~ MvNormal(mu * ones(nY), Symmetric(vcov) )     # likelihood

end


Turing.@model function gaussian_process_covars(; Y, X, D, cov_fn=exp_cov_fn, nY=length(Y), nF=size(X,2) )
    # model matrix for fixed effects (X)
    beta ~ filldist( Normal(0.0, 1.0), nF); 
    sig2 ~ LogNormal(0, 1) # "nugget" variance
    phi ~ LogNormal(0, 1) # phi is ~ lengthscale along Xstar (range parameter)
    # sigma = cov_fn(D, phi) + sig2 * LinearAlgebra.I(nY) # Realized covariance function + nugget variance
    
    mu = X * beta # mean process
    vcov = cov_fn(D, phi) + sig2 .* LinearAlgebra.I(nY )
    Y ~ MvNormal(mu, Symmetric(vcov) )     # likelihood
end



Turing.@model function gaussian_process_ar1(; Y, X, D, ar1, cov_fn=exp_cov_fn, 
    nY=length(Y), nF=size(X,2), nT=maximum(ar1)-minimum(ar1)+1 )

    alpha_ar1 ~ Normal(0,1)
    beta_ar1 ~ Normal(0,1)
    sigma_ar1 ~ LogNormal(0, 1) # "nugget" variance

    yr ~ filldist( Normal(0.0, 1.0), nT);  # -- means by time 
    # -- collects like years and 
    for n in 2:nT
        yr[n] ~ Normal(alpha_ar1 + beta_ar1 * yr[n-1], sigma_ar1 );
    end

    # # mean process model matrix  
    beta ~ filldist( Normal(0.0, 1.0), nF); 
    mu = X * beta .+ yr[ar1] 

    sig2 ~ LogNormal(0, 1) # "nugget" variance
    phi ~ LogNormal(0, 1) # phi is ~ lengthscale along Xstar (range parameter)
    # sigma = cov_fn(D, phi) + sig2 * LinearAlgebra.I(nY) # Realized covariance function + nugget variance
    # vcov = cov_fn(D, phi) + sig2 .* LinearAlgebra.I(nY )
    
    Y ~ MvNormal(mu, Symmetric(cov_fn(D, phi) .+ sig2 * eps() * I(nY) ) )     # likelihood
end






function gp_predictions(; Y, D, mu, sig2, phi, cov_fn=exp_cov_fn, nN=length(Xnew), nP=size(chain, 1) ) 
    ynew = Vector{Float64}()
    for i in sample(1:size(chain,1), nP, replace=true)
        K = cov_fn(D, phi[i])
        Koo_inv = inv(K[(nN+1):end, (nN+1):end])
        Knn = K[1:nN, 1:nN]
        Kno = K[1:nN, (nN+1):end]
        C = Kno * Koo_inv
        mvn = MvNormal( 
            C * (Y .- mu[i]) .+ mu[i], 
            Matrix(LinearAlgebra.Hermitian(Knn - C * Kno')) + sig2[i] * LinearAlgebra.I 
        ) 
        ynew = vcat(ynew, [rand(mvn) ] )
    end
    ynew = stack(ynew, dims=1)  # rehape to matrix   
    return ynew
end

