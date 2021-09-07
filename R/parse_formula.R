
parse_formula = function( fm ) {

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
      vnr = gsub( "inla.group.*[[:space:]]*[(]{1}", " ", rfac[1] )
      vnr = strsplit( vnr, ",") [[1]][1]
      vnr = gsub( "f[[:space:]]*[(]{1}|[,()=\"\']|[[:space:]]", "",  vnr )
      vnm = gsub( "[,()=\"\']|[[:space:]]", "", strsplit( rfac[2], ",")[[1]][1])
      random_effects  = rbind( random_effects, cbind( vn=vnr, model=vnm ) ) 
      grps = grepl("[[:space:]]+group[[:space:]]*[=]", vars[rt[i]] )
      if (grps) {
        vnrg = strsplit( vars[rt[i]], "[[:space:]]+group[[:space:]]*[=][[:space:]]*" )[[1]][2]
        vnrg = strsplit( vnrg, "[[:space:]]*,")[[1]][1]
        vnrgmod = strsplit( vars[rt[i]], vnrg )[[1]][2]
        vnrgac = strsplit( vnrgmod, "model[[:space:]]*=" ) [[1]][2]
        vnrgac = strsplit( vnrgac, "," ) [[1]][1]
        vnrgac = gsub( "f[[:space:]]*[(]{1}|[,()=\"\']|[[:space:]]", "",  vnrgac  )
        random_effects = rbind( random_effects, cbind( vn=vnrg, model=vnrgac, submodel="group" ) )
      }
      reps = grepl("[[:space:]]+replicate[[:space:]]*[=]", vars[rt[i]] )
      if (reps) {
        vnrr = strsplit( vars[rt[i]], "[[:space:]]+replicate[[:space:]]*[=][[:space:]]*" )[[1]][2]
        vnrr = strsplit( vnrr, "[[:space:]]*,")[[1]][1]
        vnrr = gsub( "[,()=\"\']|[[:space:]]", "",  vnrr )  
        random_effects = rbind( random_effects, cbind( vn=vnrr, model="iid", submodel="replicate" ) )
      }
    }
  }
 
  ft = setdiff(1:length(vars), rt)
  ft = setdiff( ft, c(1, which(grepl("^offset[[:space:]]*[(]{1}", vars)) ) )
  n_fixed = length(ft)

  fixed_effects = NULL
  if (n_fixed > 0) {
    for (i in 1:length(ft)) {
      fixed_effects = rbind( fixed_effects, cbind( vn=vars[ft[i]], model="fixed", submodel="main" ) ) 
    }
  }
      
  if (attributes(tfm)$intercept) {
      fixed_effects = rbind( cbind( vn="Intercept", model="fixed", submodel="NA" ), fixed_effects )
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
