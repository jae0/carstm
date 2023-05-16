
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
     
      grps = grepl("[[:space:]]+group[[:space:]]*[=]", vars[rt[i]] )

      dimensionality=NA

      if (!grps) {
        if (grepl("^space", vnr)) dimensionality = "s"
        if (grepl("^time", vnr))  dimensionality = "t"
        if (grepl("^cyclic", vnr))  dimensionality = "c"
      } 

      if (grps) {
        grpvn = strsplit( vars[rt[i]], "[[:space:]]+group[[:space:]]*[=]" )[[1]][[2]]
        grpvn = strsplit( grpvn, ",")[[1]][[1]]
        grpvn = gsub( " ", "", grpvn)

        if (grepl("^space", vnr) & grepl("^time", grpvn) ) dimensionality = "st"
        if (grepl("^time", vnr) & grepl("^space", grpvn) )  dimensionality = "ts"
        if (grepl("^space", vnr) & grepl("^cyclic", grpvn) ) dimensionality = "sc"
        if (grepl("^time", vnr) & grepl("^cyclic", grpvn) )  dimensionality = "tc"
      }

      random_effects  = rbind( random_effects, cbind( vn=vnr, model=vnm, level="main", dimensionality=dimensionality ) ) 
      
      if (grps) {

        vnrg = strsplit( vars[rt[i]], "[[:space:]]+group[[:space:]]*[=][[:space:]]*" )[[1]][2]
        vnrg = strsplit( vnrg, "[[:space:]]*,")[[1]][1]
        vnrgmod = strsplit( vars[rt[i]], vnrg )[[1]][2]
        vnrgac = strsplit( vnrgmod, "model[[:space:]]*=" ) [[1]][2]
        vnrgac = strsplit( vnrgac, "," ) [[1]][1]
        vnrgac = gsub( "f[[:space:]]*[(]{1}|[,()=\"\']|[[:space:]]", "",  vnrgac  )

        if (grepl("^space", vnrg)) dimensionality = "s"
        if (grepl("^time", vnrg))  dimensionality = "t"

        random_effects = rbind( random_effects, cbind( vn=vnrg, model=vnrgac, level="group", dimensionality=dimensionality ) )

      }
      
      reps = grepl("[[:space:]]+replicate[[:space:]]*[=]", vars[rt[i]] )
      if (reps) {
        vnrr = strsplit( vars[rt[i]], "[[:space:]]+replicate[[:space:]]*[=][[:space:]]*" )[[1]][2]
        vnrr = strsplit( vnrr, "[[:space:]]*,")[[1]][1]
        vnrr = gsub( "[,()=\"\']|[[:space:]]", "",  vnrr )  

        if (grepl("^space", vnrr)) dimensionality = "s"
        if (grepl("^time", vnrr))  dimensionality = "t"

        random_effects = rbind( random_effects, cbind( vn=vnrr, model="iid", level="replicate", dimensionality="st" ) )
      }
    }
  }
  
  ft = setdiff(1:length(vars), rt)
  ft = setdiff( ft, c(1, which(grepl("^offset[[:space:]]*[(]{1}", vars)) ) )
  n_fixed = length(ft)

  fixed_effects = NULL
  if (n_fixed > 0) {
    for (i in 1:length(ft)) {

      dimensionality = NA
      if (grepl("^space", vars[ft[i]])) dimensionality = "s"
      if (grepl("^time", vars[ft[i]]))  dimensionality = "t"

      fixed_effects = rbind( fixed_effects, cbind( vn=vars[ft[i]], model="fixed", level="main", dimensionality=dimensionality ) ) 
    }
  }
      
  if (attributes(tfm)$intercept) {
      fixed_effects = rbind( cbind( vn="Intercept", model="fixed", level="main", dimensionality="i" ), fixed_effects )
  }
      
  fixed_effects = as.data.frame(fixed_effects)

  random_effects = as.data.frame(random_effects)

  # ID  random effects associated with space/time/cycle
  # main effects

  # space
  vnS = NULL
  js = which(random_effects$dimensionality=="s" & random_effects$level=="main")
  if (length(js)==1) vnS = random_effects$vn[js]

  # time
  vnT = NULL
  jt = which(random_effects$dimensionality=="t" & random_effects$level=="main")
  if (length(jt)==1) vnT = random_effects$vn[jt]

  # cyclic
  vnU = NULL
  ju = which(random_effects$dimensionality=="c" & random_effects$level=="main")
  if (length(ju)==1) vnU = random_effects$vn[ju]

  # space copy
  vnS2 = NULL
  js2 = which(random_effects$dimensionality=="st" & random_effects$level=="main")
  if (length(js2)==1) vnS2 = random_effects$vn[js2]

  # time copy
  vnT2 = NULL
  jt2 = which(random_effects$dimensionality=="ts" & random_effects$level=="main")
  if (length(jt2)==1) vnT2 = random_effects$vn[jt2]

  vn = list(
    S = vnS,
    T = vnT,
    U = vnU,
    S2 = vnS2,
    T2 = vnT2
  )


  return(
    list(  
      formula=fm, 
      vars= vars,
      dependent_variable = dependent_variable,
      offset_variable = offset_variable,
      fixed_effects = fixed_effects,
      random_effects = random_effects,
      vn=vn
    )
  )


}
