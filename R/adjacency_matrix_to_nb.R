
adjacency_matrix_to_nb = function( W ){
    nau = nrow(W)
    nb = list()  
    for (i in 1:nau) nb[[i]] = which( W[i,] == 1 )
    return (nb)
}



