

function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end

function scottish_lip_cancer_data()
 
    # data source:  https://mc-stan.org/users/documentation/case-studies/icar_stan.html

    # y: the observed lip cancer case counts on a per-county basis
    # x: an area-specific continuous covariate that represents the proportion of the population employed in agriculture, fishing, or forestry (AFF)
    # E: the expected number of cases, used as an offset,
    # adj: a list of region ids for adjacent regions
    # num: a list of the number of neighbors for each region

    N = 56

    y   = [ 9, 39, 11, 9, 15, 8, 26, 7, 6, 20, 13, 5, 3, 8, 17, 9, 2, 7, 9, 7,
    16, 31, 11, 7, 19, 15, 7, 10, 16, 11, 5, 3, 7, 8, 11, 9, 11, 8, 6, 4,
    10, 8, 2, 6, 19, 3, 2, 3, 28, 6, 1, 1, 1, 1, 0, 0]
    
    E = [1.4, 8.7, 3.0, 2.5, 4.3, 2.4, 8.1, 2.3, 2.0, 6.6, 4.4, 1.8, 1.1, 3.3, 7.8, 4.6,
    1.1, 4.2, 5.5, 4.4, 10.5,22.7, 8.8, 5.6,15.5,12.5, 6.0, 9.0,14.4,10.2, 4.8, 2.9, 7.0,
    8.5, 12.3, 10.1, 12.7, 9.4, 7.2, 5.3,  18.8,15.8, 4.3,14.6,50.7, 8.2, 5.6, 9.3, 88.7, 
    19.6, 3.4, 3.6, 5.7, 7.0, 4.2, 1.8]
    
    x = [16,16,10,24,10,24,10, 7, 7,16, 7,16,10,24, 7,16,10, 7, 7,10,
    7,16,10, 7, 1, 1, 7, 7,10,10, 7,24,10, 7, 7, 0,10, 1,16, 0, 
    1,16,16, 0, 1, 7, 1, 1, 0, 1, 1, 0, 1, 1,16,10]
    
    # fake groups
    groups = [1, 1, 1,1 ,1,1,1, 1, 1,1,1,1,1,2,2,2,2, 2, 2,1,
    3,3,3, 1, 1, 1, 1, 1,1,1, 1,1,1, 1, 1, 1,1, 1,1, 1,
    1,1,2, 2, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1,3,2]
    
    
    adj = [ 5, 9,11,19, 7,10, 6,12, 18,20,28, 1,11,12,13,19,
    3, 8, 2,10,13,16,17, 6, 1,11,17,19,23,29, 2, 7,16,22, 1, 5, 9,12,
    3, 5,11, 5, 7,17,19, 31,32,35, 25,29,50, 7,10,17,21,22,29,
    7, 9,13,16,19,29, 4,20,28,33,55,56, 1, 5, 9,13,17, 4,18,55,
    16,29,50, 10,16, 9,29,34,36,37,39, 27,30,31,44,47,48,55,56,
    15,26,29, 25,29,42,43, 24,31,32,55, 4,18,33,45, 9,15,16,17,21,23,25,
    26,34,43,50, 24,38,42,44,45,56, 14,24,27,32,35,46,47, 14,27,31,35,
    18,28,45,56, 23,29,39,40,42,43,51,52,54, 14,31,32,37,46,
    23,37,39,41, 23,35,36,41,46, 30,42,44,49,51,54, 23,34,36,40,41,
    34,39,41,49,52, 36,37,39,40,46,49,53, 26,30,34,38,43,51, 26,29,34,42,
    24,30,38,48,49, 28,30,33,56, 31,35,37,41,47,53, 24,31,46,48,49,53,
    24,44,47,49, 38,40,41,44,47,48,52,53,54, 15,21,29, 34,38,42,54,
    34,40,49,54, 41,46,47,49, 34,38,49,51,52, 18,20,24,27,56,
    18,24,30,33,45,55]
            
    num = [4, 2, 2, 3, 5, 2, 5, 1,  6,  4, 4, 3, 4, 3, 3, 6, 6, 6 ,5, 
    3, 3, 2, 6, 8, 3, 4, 4, 4,11,  6, 7, 4, 4, 9, 5, 4, 5, 6, 5, 
    5, 7, 6, 4, 5, 4, 6, 6, 4, 9, 3, 4, 4, 4, 5, 5, 6]
    
    # areal unit id of data (y) 
    auid = 1:N  # simple here as 1:1 au:data correspondence
    nAU = N
    
    N_edges = Integer( length(adj) / 2 );
    node1 =  fill(0, N_edges); 
    node2 =  fill(0, N_edges); 
    i_adjacency = 0;
    i_edge = 0;
    for i in 1:N
    for j in 1:num[i]
        i_adjacency = i_adjacency + 1;
        if i < adj[i_adjacency]
            i_edge = i_edge + 1;
            node1[i_edge] = i;
            node2[i_edge] = adj[i_adjacency];
        end
    end
    end
    
    e = Edge.(node1, node2)
    g = Graph(e)
    W = adjacency_matrix(g)
    
    scaling_factor = scaling_factor_bym2(W)  # 0.22585017121540601 
    D = diagm(vec( sum(W, dims=2) ))
    
    # x = collect(x)
    x_scaled = (x .- mean(x)) ./ std(x)
    X = DataFrame( Intercept=ones( N ), x=x_scaled )
    X = Matrix(X)
    
    nX = size(X)[2]
    log_offset = log.(E)

    return D, W, X, log_offset, y, nX, nAU, node1, node2, scaling_factor 

end


function adjacency_matrix_to_nb( W )
    nau = size(W)[1]
    # W = LowerTriangular(W)  # using LinearAlgebra
    nb = [Int[] for _ in 1:nau]
    for i in 1:nau
        nb[i] = findall( isone, W[i,:] )
    end
    return nb
end


function nb_to_adjacency_matrix( nb )
    nau = Integer( length( unique( reduce(vcat, nb) )) )
    W = zeros( Int8, nau, nau )
    for i in 1:nau
        for j in 1:length( nb[i] )
            k = nb[i][j]
            W[i, k] = 1
        end
    end
    return(W)
end


function nodes( adj )
    nau = length(adj)
    N_edges = Integer( length( reduce(vcat, adj) )/2 )
    node1 =  fill(0, N_edges); 
    node2 =  fill(0, N_edges); 
    i_edge = 0;
    for i in 1:nau
        u = adj[i]
        num = length(u)
        for j in 1:num
            k = u[j]
            if i < k
                i_edge = i_edge + 1;
                node1[i_edge] = i;
                node2[i_edge] = k;
            end
        end
    end

    e = Edge.(node1, node2)
    g = Graph(e)
    W = Graphs.adjacency_matrix(g)
    
    # D = diagm(vec( sum(W, dims=2) ))
    scalefactor = scaling_factor_bym2(W)

    return node1, node2, scalefactor
end




function scaling_factor_bym2( adjacency_mat )
    # re-scaling variance using Reibler's solution and 
    # Buerkner's implementation: https://codesti.com/issue/paul-buerkner/brms/1241)  
    # Compute the diagonal elements of the covariance matrix subject to the 
    # constraint that the entries of the ICAR sum to zero.
    # See the inla.qinv function help for further details.
    # Q_inv = inla.qinv(Q, constr=list(A = matrix(1,1,nbs$N),e=0))  # sum to zero constraint
    # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    # scaling_factor = exp(mean(log(diag(Q_inv))))
    N = size(adjacency_mat)[1]
    asum = vec( sum(adjacency_mat, dims=2)) 
    asum = float(asum) + N .* max.(asum) .* sqrt(1e-15)  # small perturbation
    Q = Diagonal(asum) - adjacency_mat
    A = ones(N)   # constraint (sum to zero)
    S = Q \ Diagonal(A)  # == inv(Q)
    V = S * A
    S = S - V * inv(A' * V) * V'
    # equivalent form as inv is scalar
    # S = S - V / (A' * V) * V'
    scale_factor = exp(mean(log.(diag(S))))
    return scale_factor
 
end



function scaling_factor_bym2(node1, node2, groups=ones(length(node1))) 
    ## calculate the scale factor for each of k connected group of nodes, 
    ## copied from the scale_c function from M. Morris
    gr = unique( groups )
    n_groups = length(gr)
    scale_factor = ones(n_groups)
    for j in 1:n_groups 
      k = findall( x -> x==j, groups)
      if length(k) > 1 
        e = Edge.(node1[k], node2[k])
        g = Graph(e)
        adjacency_mat = adjacency_matrix(g)
        scale_factor[j] = scaling_factor_bym2( adjacency_mat )
      end
    end
    return scale_factor
end
  




Turing.@model function turing_car(D, W, X, log_offset, y, nX )
    # base model .. slow, MVN of var-covariance matrix
    alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
    tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
    beta ~ filldist( Normal(0.0, 1.0), nX );
    prec = tau .* (D - alpha .* W)
    if !isposdef(prec)
        # check postive definiteness
        # phi ~ MvNormal( zeros(N) ); 
        Turing.@addlogprob! -Inf
        return nothing
    end
    sigma = inv( Symmetric(prec) )
    # sigma = Symmetric( prec) \ Diagonal(ones(nX))  # alternatively
    phi ~ MvNormal( sigma );  # mean zero
    lambda = exp.( X * beta .+ phi .+ log_offset )
    @. y ~ Poisson( lambda );
end

   
Turing.@model function turing_car_prec(D, W, X, log_offset, y, nX )
    # MVN of precision matrix
    alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
    tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
    beta ~ filldist( Normal(0.0, 1.0), nX);
    prec = tau .* (D - alpha .* W)
    if !isposdef(prec)
        # check postive definiteness
        # phi ~ MvNormal( zeros(N) ); 
        Turing.@addlogprob! -Inf
        return nothing
    end
    phi ~ MvNormalCanon( Symmetric(prec) );  # mean zero .. no inverse
    lambda = exp.( X * beta .+ phi .+ log_offset )
    @. y ~ Poisson( lambda );

end

 

Turing.@model function turing_icar_direct_test( node1, node2; ysd=std(skipmissing(y)), nY=size(y,1)  )
    # equivalent to Morris' "simple_iar' .. testing pairwise difference formulation
    # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)

    phi ~ filldist( Normal(0.0, ysd), nY)   # 10 is std from data: std(y)=7.9 stan goes from U(-Inf,Inf) .. not sure why 
    dphi = phi[node1] - phi[node2]
    lp_phi =  -0.5 * dot( dphi, dphi )
    Turing.@addlogprob! lp_phi
    
    # soft sum-to-zero constraint on phi)
    # equivalent to mean(phi) ~ normal(0,0.001)
    sum_phi = sum(phi)
    sum_phi ~ Normal(0, 0.001 * nY);  
  
    # no data likelihood -- just prior sampling  -- 
end

  
Turing.@model function turing_icar_direct_bym( X, log_offset, y, nX, node1, node2; ysd=std(skipmissing(y)),  nY=size(X,1) )
    # BYM
    # alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
     # tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
     beta ~ filldist( Normal(0.0, 5.0), nX);
     theta ~ filldist( Normal(0.0, 1.0), nY) # unstructured (heterogeneous effect)
     # phi ~ filldist( Laplace(0.0, ysd), nY) # spatial effects: stan goes from -Inf to Inf .. 
     phi ~ filldist( Normal(0.0, ysd), nY) # spatial effects: stan goes from -Inf to Inf .. 
 
     # pairwise difference formulation ::  prior on phi on the unit scale with sd = 1
     # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
     dphi = phi[node1] - phi[node2]
     lp_phi =  -0.5 * dot( dphi, dphi )
     Turing.@addlogprob! lp_phi
     
     # soft sum-to-zero constraint on phi)
     # equivalent to mean(phi) ~ normal(0, 0.001)
     sum_phi = sum(phi)
     sum_phi ~ Normal(0, 0.001 * nY);  

     tau_theta ~ Gamma(3.2761, 1.0/1.81);  # Carlin WinBUGS priors
     tau_phi ~ Gamma(1.0, 1.0);            # Carlin WinBUGS priors

     sigma_theta = inv(sqrt(tau_theta));  # convert precision to sigma
     sigma_phi = inv(sqrt(tau_phi));      # convert precision to sigma

     lambda = exp.( X * beta .+ phi .* sigma_phi .+ theta .* sigma_theta .+ log_offset )
  
     @. y ~ Poisson( lambda );
end
 

Turing.@model function turing_icar_direct_bym2( X, log_offset, y, auid, nX, nAU, node1, node2, scaling_factor )
    beta ~ filldist( Normal(0.0, 1.0), nX);
    theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
    phi ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 
    # pairwise difference formulation ::  prior on phi on the unit scale with sd = 1
    # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
    dphi = phi[node1] - phi[node2]
    Turing.@addlogprob! -0.5 * dot( dphi, dphi )
    # soft sum-to-zero constraint on phi)
    sum_phi = sum(phi)
    sum_phi ~ Normal(0, 0.001 * nAU);  
    sigma ~ truncated( Normal(0, 1.0), 0, Inf) ; 
    rho ~ Beta(0.5, 0.5);
    # variance of each component should be approximately equal to 1
    convolved_re =  sigma .*  ( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi );
    lambda = exp.( X * beta +  convolved_re[auid] + log_offset )
    @. y ~ Poisson( lambda );
end
 

Turing.@model function turing_icar_direct_bym2_binomial(y, Ntrial, X, nX, nAU, node1, node2, scaling_factor)
    # poor form to not pass args directly but simpler to use global vars as this is a one-off:
    beta0 ~ Normal(0.0, 1.0)
    betas ~ filldist( Normal(0.0, 1.0), nX); #coeff
    theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
    phi ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 
    dphi = phi[node1] - phi[node2]
    Turing.@addlogprob! (-0.5 * dot( dphi, dphi )) # directly add to logprob
    sum_phi = sum(phi)
    sum_phi ~ Normal(0.0, 0.001 * nAU);      # soft sum-to-zero constraint on phi), equivalent to 
    sigma ~ Gamma(1.0, 1.0)
    rho ~ Beta(0.5, 0.5);
    # variance of each component should be approximately equal to 1
    convolved_re =  sqrt(1 - rho) .* theta .+ sqrt(rho / scaling_factor) .* phi ;
    v = beta0 .+ X * betas .+ sigma .* convolved_re 
    y ~ arraydist(LazyArray(@~ BinomialLogit.(Ntrial, v)))  # 100 sec
end



Turing.@model function turing_icar_direct_bym2_groups( X, log_offset, y, auid, nX, nAU, node1, node2, scaling_factor, groups, gi; ysd=std(skipmissing(y)) )
    # BYM2
    # alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
     # tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
     beta ~ filldist( Normal(0.0, 5.0), nX);
     theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
     phi ~ filldist(Normal(0.0, ysd), nAU) # spatial effects: stan goes from -Inf to Inf .. 
        
     # pairwise difference formulation ::  prior on phi on the unit scale with sd = 1
     # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
     dphi = phi[node1] - phi[node2]
     lp_phi =  -0.5 * dot( dphi, dphi )
     Turing.@addlogprob! lp_phi
     
     sigma ~ truncated( Normal(0, 1.0), 0, Inf) ; 
     rho ~ Beta(0.5, 0.5);

     convolved_re = zeros(nAU)

     for j in 1:length(gi)
         ic = gi[j] 
        
         # soft sum-to-zero constraint on phi)
         # equivalent to mean(phi) ~ normal(0, 0.001)
         sum_phi = sum(phi[ic])
         sum_phi ~ Normal(0, 0.001 * nAU);  

         if  length(ic) == 1 
             convolved_re[ ic ] = sigma .* theta[ ic ];
         else  
             convolved_re[ ic ] = sigma .* ( sqrt.(1 .- rho) .* theta[ ic ]  +  sqrt(rho ./ scaling_factor[j] )  .* phi[ ic ] ) ;
         end 
     end
  
     # convolved_re =  sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi;
   
     lambda = exp.( X * beta .+  convolved_re[auid] .+ log_offset )
   
     @. y ~ Poisson( lambda );
  
    # to compute from posteriors
    #  real logit_rho = log(rho / (1.0 - rho));
    #  vector[N] eta = log_E + beta0 + x * betas + convolved_re * sigma; // co-variates
    #  vector[N] lambda = exp(eta);
end 


Turing.@model function turing_icar_latent_bym2( X, log_offset, y, auid, nX, nAU, node1, node2, 
    scaling_factor )

    beta ~ filldist( Normal(0.0, 1.0), nX); 
    theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
    phi ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 
    dphi = phi[node1] - phi[node2]
    Turing.@addlogprob! (-0.5 * dot( dphi, dphi ))
    sum_phi = sum(phi) 
    sum_phi ~ Normal(0, 0.001 * nAU);      # soft sum-to-zero constraint on phi)
    sigma ~ truncated( Normal(0.0, 1.0), 0, Inf) ; 
    rho ~ Beta(0.5, 0.5);

    # spatial effects:  nAU
    convolved_re = sigma .*( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi )

    lambda =  X * beta +  convolved_re[auid] + log_offset 
    # if all( lambda .< 42.0) & all( lambda .> -42.0 ) # on log scale .. forces results to also to be finite)
        # @show summarystats(lambda)
        # equivalent ways of expressing likelihood:
        @. y ~ LogPoisson( lambda);
        # y ~ arraydist([LogPoisson( lambda[i] ) for i in 1:nY ])
        # y ~ arraydist(LazyArray(Base.broadcasted((l) -> LogPoisson(l), lambda)))
        # y ~ arraydist(LazyArray( @~ LogPoisson.(lambda) ) )
    #else 
    #    Turing.@addlogprob! -Inf
    #end
    return nothing
end




Turing.@model function turing_icar_latent_bym2_gp1( X, G, log_offset, y, z, auid, nY, nX, nG, nAU, node1, node2, scaling_factor )
    # one kernel for all GP
    beta ~ filldist( Normal(0.0, 1.0), nX);
    theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
    phi ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 
    dphi = phi[node1] - phi[node2]
    Turing.@addlogprob! (-0.5 * dot( dphi, dphi ))
    sum_phi = sum(phi) 
    sum_phi ~ Normal(0, 0.001 * nAU);      # soft sum-to-zero constraint on phi)
    sigma ~ truncated( Normal(0.0, 1.0), 0, Inf) ; 
    rho ~ Beta(0.5, 0.5);
    # variance of each component should be approximately equal to 1
    convolved_re = sigma .*( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi )
     
    lambda0 =  X * beta + convolved_re[auid] + log_offset  #non GP components

    kernel_var ~ Gamma(2.0, 0.5) 
    kernel_scale ~ Gamma(2.0, 0.1) # even more left shifted with mode ~0.1 .. ~ magnitudes of ~ 0.1 seems to be stable
    l2reg ~ Gamma(2.0, 0.001)   # = 1e-4  # L2 regularization factor for ridge regression
    # plot(x->pdf(Beta(2.0, 100.), x), xlim=(0, 1))  # small value >0 ~ 1e-4
    # plot(x->pdf(Gamma(2, 0.5), x), xlim=(0, 10))   # centered over 1 variance prior; kernel functions models scaled and centered data
 
    k = ( kernel_var * SqExponentialKernel() ) ∘ ScaleTransform(kernel_scale) 

    kmat = Symmetric( kernelmatrix( k, G, obsdim=1 ) + l2reg*I ) # makes PD and cholesky, add regularization
    chkm = cholesky(kmat)
 
    # "latent" GP of G on unit scales
    mean_process = kmat * ( chkm.U \ (chkm.L \ (z-lambda0)) )
    
    # gps =  chkm.L * rand(Normal(0.0, 1.0), nY) + mean_process
    gps = rand( MvNormal( mean_process, kmat ) ) # faster

    lambda = lambda0 + sum(gps, dims=2)
    @. y ~ LogPoisson( lambda );

    return nothing
end



Turing.@model function turing_icar_latent_bym2_gp2( X, G, log_offset, y, z, auid, nY, nX, nG, nAU, node1, node2, scaling_factor )
    # product of separate kernels for all GP
 
    beta ~ filldist( Normal(0.0, 1.0), nX);
    theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
    phi ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 
    dphi = phi[node1] - phi[node2]
    Turing.@addlogprob! (-0.5 * dot( dphi, dphi ))
    sum_phi = sum(phi) 
    sum_phi ~ Normal(0, 0.001 * nAU);      # soft sum-to-zero constraint on phi)
    sigma ~ truncated( Normal(0.0, 1.0), 0, Inf) ; 
    rho ~ Beta(0.5, 0.5);
    # variance of each component should be approximately equal to 1
    convolved_re = sigma .*( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi )
    lambda0 =  X * beta + convolved_re[auid] + log_offset  # non GP component

    # centered over 1 variance prior; kernel functions models scaled and centered data
    # plot(x->pdf(Gamma(2, 0.5), x), xlim=(0, 10))   
    kernel_var ~ filldist( Gamma(2.0, 0.5), nG)
    
    # mode ~0.1 .. ~ magnitudes of ~ 0.1 seems to be stable
    kernel_scale ~ filldist( Gamma(2.0, 0.1), nG)
    
    # L2 regularization factor for ridge regression
    l2reg ~  Gamma(2.0, 0.001)
  
    k = prod( [ ( kernel_var[i] * SqExponentialKernel() ) ∘ ScaleTransform(kernel_scale[i]) for i in 1:nG ] )
  
    kmat = kernelmatrix( k, G, obsdim=1 ) + l2reg*I # add regularization
    chkm = cholesky(kmat) 
    # "latent" GP of G on unit scale 
    
    mean_process = kmat * ( chkm.U \ (chkm.L \ (z-lambda0) ) )

    # gps = chkm.L * rand(Normal(0.0, 1.0), size(G)) .+ mean_process
    gps = rand( MvNormal( mean_process, kmat ) ) # faster
    
    lambda = lambda0 + sum(gps, dims=2) 
    
    # equivalent ways of expressing likelihood:
    @. y ~ LogPoisson( lambda );

    return nothing
end



Turing.@model function turing_icar_latent_bym2_gp_hurdle(
    X, G, log_offset, y, pa, wt, z, auid, 
    nY, nX, nG, nAU, 
    node1, node2, scaling_factor, good )

    # one kernel for all GP
    beta ~ filldist( Normal(0.0, 1.0), nX);
    beta_hab ~ filldist( Normal(0.0, 1.0), nX);

    theta ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)
    theta_hab ~ filldist( Normal(0.0, 1.0), nAU)  # unstructured (heterogeneous effect)

    phi ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 
    phi_hab ~ filldist( Normal(0.0, 1.0), nAU) # spatial effects: stan goes from -Inf to Inf .. 

    dphi = phi[node1] - phi[node2]
    dphi_hab = phi[node1] - phi[node2]

    Turing.@addlogprob! (-0.5 * dot( dphi, dphi ))
    Turing.@addlogprob! (-0.5 * dot( dphi_hab, dphi_hab ))

    sum_phi = sum(phi) 
    sum_phi_hab = sum(phi_hab) 

    sum_phi ~ Normal(0, 0.001 * nAU);      # soft sum-to-zero constraint on phi)
    sum_phi_hab ~ Normal(0, 0.001 * nAU);      # soft sum-to-zero constraint on phi)

    sigma ~ truncated( Normal(0.0, 1.0), 0, Inf) ; 
    sigma_hab ~ truncated( Normal(0.0, 1.0), 0, Inf) ; 

    rho ~ Beta(0.5, 0.5);
    rho_hab ~ Beta(0.5, 0.5);

    # variance of each component should be approximately equal to 1
    convolved_re = sigma .*( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi )
    convolved_re_hab = sigma_hab .*( sqrt.(1 .- rho_hab) .* theta_hab .+ sqrt.(rho_hab ./ scaling_factor) .* phi )
     
    lambda0 =  X * beta + convolved_re[auid] + log_offset  #non GP components
    pr_habitat0 = X * beta_hab + convolved_re_hab[auid]  

    kernel_var ~ Gamma(2.0, 0.5) 
    kernel_var_hab ~ Gamma(2.0, 0.5) 

    kernel_scale ~ Gamma(2.0, 0.1) # even more left shifted with mode ~0.1 .. ~ magnitudes of ~ 0.1 seems to be stable
    kernel_scale_hab ~ Gamma(2.0, 0.1) # even more left shifted with mode ~0.1 .. ~ magnitudes of ~ 0.1 seems to be stable

    l2reg ~ Gamma(2.0, 0.001)   # = 1e-4  # L2 regularization factor for ridge regression
    l2reg_hab ~ Gamma(2.0, 0.001)   # = 1e-4  # L2 regularization factor for ridge regression

    # plot(x->pdf(Beta(2.0, 100.), x), xlim=(0, 1))  # small value >0 ~ 1e-4
    # plot(x->pdf(Gamma(2, 0.5), x), xlim=(0, 10))   # centered over 1 variance prior; kernel functions models scaled and centered data
 
    k = ( kernel_var * SqExponentialKernel() ) ∘ ScaleTransform(kernel_scale) 
    k_hab = ( kernel_var_hab * SqExponentialKernel() ) ∘ ScaleTransform(kernel_scale_hab) 

    kmat = kernelmatrix( k, G, obsdim=1 ) + l2reg*I # makes PD and cholesky, add regularization
    kmat_hab = kernelmatrix( k_hab, G, obsdim=1 ) + l2reg_hab*I # makes PD and cholesky, add regularization

    chkm = cholesky(kmat)
    chkm_hab = cholesky(kmat_hab)
 
    # "latent" GP of G on unit scales
    mean_process = kmat * ( chkm.U \ (chkm.L \ (z-lambda0)) )
    mean_process_hab = kmat_hab * ( chkm_hab.U \ (chkm_hab.L \ (pa-pr_habitat0)) )
    
    # gps =  chkm.L * rand(Normal(0.0, 1.0), nY) + mean_process
    gps = rand( MvNormal( mean_process, kmat ) ) # faster
    gps_hab = rand( MvNormal( mean_process_hab, kmat_hab ) ) # faster

    lambda = lambda0 + sum(gps, dims=2)
    pr_habitat = pr_habitat0 + sum(gps_hab, dims=2)

    # Hurdle process
    @. pa ~ Bernoulli( logistic(pr_habitat[auid]) )   
    @. y[good] ~ truncated( Poisson(  exp(lambda[good]) ),  2, nothing  ) ; # 1 or less is considered poor habitat
 #=
    # dynamics
    K ~ filldist( TruncatedNormal( PM.K[1], PM.K[2], PM.K[3], PM.K[4]), nAU)

    r ~  TruncatedNormal( PM.r[1], PM.r[2], PM.r[3], PM.r[4])   # (mu, sd)
    
    bpsd ~  TruncatedNormal( PM.bpsd[1], PM.bpsd[2], PM.bpsd[3], PM.bpsd[4] )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( PM.bosd[1], PM.bosd[2], PM.bosd[3], PM.bosd[4] )  ;  # slightly informative .. center of mass between (0,1)
    q ~ TruncatedNormal( PM.q[1], PM.q[2], PM.q[3], PM.q[4] )    
    qc ~ TruncatedNormal(PM.qc[1], PM.qc[2], PM.qc[3], PM.qc[4]  ) 

    m = tzeros( PM.nM, nAU )

    @.  m[1,:] ~ TruncatedNormal( PM.m0[1], PM.m0[2], PM.m0[3], PM.m0[4] )  ; # starting b prior to first catch event

    for i in 2: PM.nT
        for j in 1:nAU
            m[i,j] ~ TruncatedNormal( m[i-1,j] + r * m[i-1,j] * ( 1.0 - m[i-1,j] ) -  PM.removed[i-1,j]/K[j], bpsd, PM.mlim[1], PM.mlim[2])  ;
        end
    end
 
    if any( x -> x < 0.0 || x >1.0, m)
        Turing.@addlogprob! -Inf
        return nothing
    end
    
    # likelihood
    # observation model: Y = q X + qc ; X = (Y - qc) / q
    # m = abundance (prefishery) 
    PM.S[i] ~ Normal( q * ( m[i] - PM.removed[i]/K )+ qc, bosd )  ; # fall survey
=#

    return nothing
end

#=
m = turing_icar_latent_bym2_gp(
  X = X[io,:], 
  log_offset=log.(M.data_offset[io]), 
  y=M.totno[io], 
  auid=auid[io], 
  node1=node1, 
  node2=node2, 
  scaling_factor=scaling_factor 
) 

o = sample(m, Turing.NUTS( 5, 0.65 ), 5) 
showall( summarize(o) )
=#


function turing_icar_latent_bym2_predict( res, X; log_offset=0, scaling_factor=1.0, n_sample=-1, nAU, auid  )

    nchains = size(res)[3]
    nsims = size(res)[1]
    
    if n_sample == -1
        n_sample = nchains * nsims         # do all
    end

    # md = zeros(nM, nS, n_sample, 2)  # number normalized

    NX = size(X)[2]
    
    ntries_mult=2
    ntries = 0
    z = 0

    while z <= n_sample 
        ntries += 1
        ntries > ntries_mult * n_sample && break 
        z >= n_sample && break

        j = rand(1:nsims)  # nsims
        l = rand(1:nchains) # nchains

        sigma = res[j, Symbol("sigma"), l] 
        rho   = res[j, Symbol("rho"), l] 
        theta = [ res[j, Symbol("theta[$k]"), l] for k in 1:nAU]
        phi   = [ res[j, Symbol("phi[$k]"), l]   for k in 1:nAU]
        beta  = [ res[j, Symbol("beta[$k]"), l]  for k in 1:NX]
    
        convolved_re =  sigma .*  ( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi );
        lambda = exp.( X * beta +  convolved_re[auid] .+ log_offset )
        # @show summarystats(lambda)
        y = rand.(Poisson.(lambda));

        z += 1
        
    end  # while

    if z < n_sample 
    @warn  "Insufficient number of solutions" 
    end

    return md, mn, mb, trace, trace_bio, trace_time 

end
