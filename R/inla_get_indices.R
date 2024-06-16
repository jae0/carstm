
inla_get_indices = function( X, model="", tag, start, len ) {

    ip = list()

    if (model=="direct_match")  {
        jp = which( tag %in% X )
        ip = list()
        for (i in 1:(length(jp))) {
            jj = jp[[i]]
            ip[[i]] = start[jj] + 0:(len[jj]-1)
        }
        return(ip) 
    }

    jp = list()
    for (m in 1:length(X)) {
        jp[[m]] = grep(X[[m]], tag 
        )
    }

    if (length(jp)>0) {
        ip = list()
        for (i in 1:(length(jp))) {
            jj = jp[[i]]
            ip[[i]] = start[jj] + 0:(len[jj]-1)
        }
    }

    if (model=="")  return(ip)

    if (model %in% c("bym2", "bym") ) {
        # bym2: Z = expand.grid( space=modelinfo$space, type = c("re_total", "re_neighbourhood"), time=modelinfo$time, stringsAsFactors =FALSE ) 
        ns = length(ip[[1]])/2
        ibym2 = list(
            re_total = ip[[1]][1:ns],   # random effects = neighbourhood effects + unstructured effects 
            re_neighbourhood = ip[[1]][(ns+1):(ns*2)]  ## neighbourhood effects 
        )
        return(ibym2)
    }

}
