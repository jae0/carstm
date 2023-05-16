carstm_results_unpack = function( res, vn ) {
  # to facillitate data extraction depending upon vn data hierarchy 
  return( switch( length(vn),
    res[[ vn[1] ]],
    res[[ vn[1] ]] [[vn[2] ]],
    res[[ vn[1] ]] [[vn[2] ]] [[ vn[3] ]],
    res[[ vn[1] ]] [[vn[2] ]] [[ vn[3] ]] [[vn[4] ]],
    res[[ vn[1] ]] [[vn[2] ]] [[ vn[3] ]] [[vn[4] ]] [[vn[5] ]]
  )) 
}
