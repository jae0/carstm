
scaling_factor_bym2 = function(adjacency_mat) {
    # re-scaling variance using Reibler's solution and 
    # Buerkner's implementation: https://codesti.com/issue/paul-buerkner/brms/1241)  
    # Compute the diagonal elements of the covariance matrix subject to the 
    # constraint that the entries of the ICAR sum to zero.
    # See the inla.qinv function help for further details.
    # Q_inv = inla.qinv(Q, constr=list(A = matrix(1,1,nbs$N),e=0))  # sum to zero constraint
    # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    # scaling_factor = exp(mean(log(diag(Q_inv))))

    N = nrow(adjacency_mat) 
    asum = rowSums(adjacency_mat)  
    asum = asum + N * max(asum) * sqrt(1e-15)   # small perturbation to make PD
    Q = diag(asum) - adjacency_mat
    S = solve(Q)
    A = matrix(1,N,1)   # constraint (sum to zero)
    V = S %*% A
    S = S - V %*% solve(t(A) %*% V) %*% t(V)
    scale_factor = exp(mean(log(diag(S))))
    return (scale_factor)
    
    if (0) {
        # test with STAN's version

        # data:  
        require("SpatialEpi")
        data(scotland)
        nb <- poly2nb(scotland$spatial.polygon, snap = 1)
        C <- nb2mat(nb, style = "B", zero.policy=T)


        # above function for group 1:
        group1 = setdiff( 1:N, c(6,8,11)) # to match test subset
        scaling_factor_bym2(C[group1, group1])  # 0.5578

        # using funcs from STAN
        scale_c <- function(C) {
            #' compute geometric mean of a vector
            geometric_mean <- function(x) exp(mean(log(x))) 
            
            N = dim(C)[1]
            
            # Create ICAR precision matrix  (diag - C): this is singular
            # function Diagonal creates a square matrix with given diagonal
            Q =  Diagonal(N, rowSums(C)) - C
            
            # Add a small jitter to the diagonal for numerical stability (optional but recommended)
            Q_pert = Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
            
            # Function inla.qinv provides efficient way to calculate the elements of the
            # the inverse corresponding to the non-zero elements of Q
            Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0))
            
            # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
            scaling_factor <- geometric_mean(Matrix::diag(Q_inv)) 
            return(scaling_factor) 
            }
            prep_icar_data <- function (C, inv_sqrt_scale_factor = NULL) {
            n <- nrow(C)
            E <- edges(C)
            G <- list(np = nrow(C), from = E$node1, to = E$node2, nedges = nrow(E))
            class(G) <- "Graph"
            nb2 <- spdep::n.comp.nb(spdep::graph2nb(G))
            k = nb2$nc
            if (inherits(inv_sqrt_scale_factor, "NULL")) inv_sqrt_scale_factor <- array(rep(1, k), dim = k)
            group_idx = NULL
            for (j in 1:k) group_idx <- c(group_idx, which(nb2$comp.id == j))
            group_size <- NULL
            for (j in 1:k) group_size <- c(group_size, sum(nb2$comp.id == j))
            # intercept per connected component of size > 1, if multiple.
            m <- sum(group_size > 1) - 1
            if (m) {
                GS <- group_size
                ID <- nb2$comp.id
                change.to.one <- which(GS == 1)
                ID[which(ID == change.to.one)] <- 1
                A = model.matrix(~ factor(ID))
                A <- as.matrix(A[,-1])
            } else {
                A <- model.matrix(~ 0, data.frame(C))
            }
            l <- list(k = k, 
                        group_size = array(group_size, dim = k), 
                        n_edges = nrow(E), 
                        node1 = E$node1, 
                        node2 = E$node2, 
                        group_idx = array(group_idx, dim = n), 
                        m = m,
                        A = A,
                        inv_sqrt_scale_factor = inv_sqrt_scale_factor, 
                        comp_id = nb2$comp.id)
            return(l)
        }

        edges <- function (w) {
            lw <- apply(w, 1, function(r) {
                which(r != 0)
            })
            all.edges <- lapply(1:length(lw), function(i) {
                nbs <- lw[[i]]
                if (length(nbs)) 
                data.frame(node1 = i, node2 = nbs, weight = w[i, nbs])
            })
            all.edges <- do.call("rbind", all.edges)
            edges <- all.edges[which(all.edges$node1 < all.edges$node2), ]
            return(edges)
        }


        icar.data <- prep_icar_data(C)
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
        scale_factor  # [1] 0.5578 1.0000 1.0000 1.0000

    }
 }
