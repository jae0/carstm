
carstm_parse_formula = function( fm ) {

  tfm = terms.formula(fm, specials = c("f"), data = NULL)
  vars =  setdiff( as.character(attr(tfm, "variables")), "list" ) 
  terms = attr(tfm, "term.labels")
  nt = length(terms)

  dependent_variable = NULL
  if (attr(tfm, "response") > 0)  dependent_variable =  gsub( "[,()=\"\']|[[:space:]]", "", vars[1] )

  offset_variable = NULL
  if (length(attr(tfm, "offset")) > 0) {
    # remove "offset()"
    offset_variable = vars[grep("^offset[[:space:]]*[(]{1}", vars)]
    offset_variable = gsub( "offset[[:space:]]*[(]{1}|[)]{1}$", "", offset_variable )  
    offset_variable = gsub( "[,()=\"\']|[[:space:]]", "",  offset_variable )  
  }

  # random effects ("f()")
  rt = attr(tfm, "specials")$f  # random components
  n_random = length(rt)

  random_effects = NULL
  if (n_random > 0) {
    for (i in 1:length(rt)) {
      rfac = strsplit( vars[rt[i]], "model[[:space:]]*=" )[[1]]
      vnr = gsub( "f[[:space:]]*[(]{1}|[,()=\"\']|[[:space:]]", "",  rfac[1] )
      vnm = gsub( "[,()=\"\']|[[:space:]]", "", strsplit( rfac[2], ",")[[1]][1])
      random_effects  = rbind( random_effects, cbind( vn=vnr, model=vnm ) ) 
    }
  }

  ft = setdiff(1:length(vars), rt)
  ft = setdiff( ft, c(1, which(grepl("^offset[[:space:]]*[(]{1}", vars)) ) )
  n_fixed = length(ft)

  fixed_effects = NULL
  if (n_fixed > 0) {
    for (i in 1:length(ft)) {
      fixed_effects = rbind( fixed_effects, cbind( vn=vars[ft[i]], model="fixed" ) ) 
    }
  }
      
  if (attributes(tfm)$intercept) {
      fixed_effects = rbind( fixed_effects, cbind( vn="Intercept", model="fixed" ) )
  }
      
  return(
    list(  
      formula=fm, 
      vars= vars,
      dependent_variable = dependent_variable,
      offset_variable = offset_variable,
      fixed_effects = as.data.frame(fixed_effects),
      random_effects = as.data.frame(random_effects)
    )
  )
}
