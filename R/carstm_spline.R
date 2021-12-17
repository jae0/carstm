
carstm_spline = function( res, vn, statvar="mean", method="monoH.FC") {
  dta = carstm_results_unpack( res, vn ) 
  out = splinefun( dta$ID,  dta[,statvar], method=method)
  return(out)
}
