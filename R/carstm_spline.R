
carstm_spline = function( res, vn, statvar="mean", method="monoH.FC") {
  dta = carstm_results_unpack( res, vn ) 
  if (!exists("ID", dta)) {
    dta$ID = row.names(dta)
  }
  out = splinefun( dta$ID,  dta[,statvar], method=method)
  return(out)
}
