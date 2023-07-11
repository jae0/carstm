nb_to_adjacency_matrix = function( nb ){
    nau = length(nb)
    W = matrix(0, nau, nau )
    for (i in 1:nau) {
        for (j in 1:length( nb[[i]] )) {
            k = nb[[i]][[j]]
            W[i, k] = 1
        }
    }
    return(W)
}