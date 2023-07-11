
nodes = function( nb ){
    # unravel as nodes to speed up computations
    nau = length(nb)
    N_edges = sum( unlist( lapply(nb, length)) )/ 2

    node1 =  rep(0, N_edges); 
    node2 =  rep(0, N_edges); 
    i_edge = 0;

    for (i in 1:nau) {
        u = nb[[i]]
        num = length(u)
        for (j in 1:num) {
            k = u[j]
            if (i < k) {
                i_edge = i_edge + 1;
                node1[i_edge] = i;
                node2[i_edge] = k;
            }
        }
    }

    adj = nb_to_adjacency_matrix(nb)
    
    scale_factor = scaling_factor_bym2(adj)

    return (list(node1=node1, node2=node2, scale_factor=scale_factor))
}
