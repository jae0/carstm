carstm_results_unpack = function( res, vn ) {
  # to facillitate data extraction depending upon vn data hierarchy 

  vn_len = length(vn)

  if (vn_len == 0) return( res)
  if (vn_len == 1) {
    if (is.null(names(res))) return( res)
  }

  switch( vn_len,
    res[[ vn[1] ]],
    res[[ vn[1] ]] [[vn[2] ]],
    res[[ vn[1] ]] [[vn[2] ]] [[ vn[3] ]],
    res[[ vn[1] ]] [[vn[2] ]] [[ vn[3] ]] [[vn[4] ]],
    res[[ vn[1] ]] [[vn[2] ]] [[ vn[3] ]] [[vn[4] ]] [[vn[5] ]]
  ) 
}
