

  substrate_carstm = function( p=NULL, DS="parameters", redo=FALSE, ... ) {

    require( carstm)

    if ( is.null(p)) {
      p = aegis.substrate::substrate_parameters(...)
    } else {
      p = aegis.substrate::substrate_parameters(p=p, ...)
    }

  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.polygons", "aegis.coastline", "aegis.substrate", "carstm"   ) )

  # ------------------


    if (DS=="parameters_override") {
    # translate param values from one project to a unified representation
    # must be first to catch p
    pc = substrate_carstm(
      DS = "parameters",
      project_name = "substrate",
      variabletomodel = "substrate.grainsize",
      spatial_domain = p$spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_overlay = p$areal_units_overlay, # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = p$areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = p$areal_units_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      modeldir = p$modeldir,  # outputs all go the the main project's model output directory
      auid = p$auid
    )
    return(pc)
  }


  # ---------------------



  if (DS=="parameters") {

    p$libs = unique( c( p$libs, project.library ( "carstm" ) ) )


    if ( !exists("areal_units_strata_type", p)) p$areal_units_strata_type = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( p$spatial_domain == "SSE" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "groundfish_strata" #.. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( p$spatial_domain == "snowcrab" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla.default"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ("INLA" ) ) )
        p$carstm_model_label = "production"
        p$carstm_modelcall = paste('
          inla(
            formula =', p$variabletomodel, ' ~ 1
              + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
              + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "lognormal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # num.threads=4,
            # blas.num.threads=8,
            verbose=TRUE
          ) ' )
      }
      if ( grepl("glm", p$carstm_modelengine) ) {
      p$carstm_model_label = "default_glm"
       p$carstm_modelcall = paste('glm( formula =', p$variabletomodel, '~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  )  ' )  # for modelengine='glm'
      }
      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv" ) ) )
        p$carstm_model_label = "default_gam"
        p$carstm_modelcall = paste('gam( formula =', p$variabletomodel, '~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) ' ) # for modelengine='gam'
      }
    }

    p = carstm_parameters( p=p, DS="basic" )

    return(p)
  }



  # ------------------




  if ( DS=="carstm_inputs") {

    aggregate_data = FALSE
    if (exists("carstm_inputs_aggregated", p)) {
      if (p$carstm_inputs_aggregated)  aggregate_data = TRUE
    }


    if (aggregate_data) {
      fn = file.path( p$modeldir, paste( "substrate", "carstm_inputs", p$auid,
        p$inputdata_spatial_discretization_planar_km,
        "rdata", sep=".") )
    } else {
      fn = file.path( p$modeldir, paste( "substrate", "carstm_inputs", p$auid,
        "rawdata",
        "rdata", sep=".") )
    }


    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))


    # do this immediately to reduce storage for sppoly (before adding other variables)

    if (aggregate_data) {
      M = substrate.db ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!
      names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel
    } else {
      M = substrate.db( p=p, DS="lonlat.highres" )
      names(M)[which(names(M)=="grainsize" )] = p$variabletomodel
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

      # p$quantile_bounds_data = c(0.0005, 0.9995)
      if (exists("quantile_bounds_data", p)) {
        TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds_data, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
        keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
        # this was -1.7, 21.8 in 2015
      }

    }

    M = M[ which(is.finite(M[, p$variabletomodel] )), ]

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(substrate.grainsize.mean~plon+plat, data=M, aspect="iso")
    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
    M = M[ which(is.finite(M$StrataID)),]

    pB = bathymetry_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry

    M[, pB$variabletomodel] = lookup_bathymetry_from_surveys( p=pB, locs=M[, c("lon", "lat")] )

    # if any still missing then use a randomly chosen depth by StrataID
    kk =  which( !is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0) {
      AD = bathymetry.db ( p=pB, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$StrataID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      AD = AD[ which(is.finite(AD$StrataID)),]
      oo = tapply( AD[, paste(pB$variabletomodel, "mean", sep="." )], AD$StrataID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$StrataID[kk]), as.character( names(oo )) )
      M[kk, pB$variabletomodel] = oo[jj ]
    }


    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M$tag = "observations"

    sppoly_df = as.data.frame(sppoly)

    pB$carstm_model_label = "production"
    BM = carstm_model ( p=pB, DS="carstm_modelled" )  # modeled!
    kk = match( as.character(  sppoly_df$StrataID), as.character( BM$StrataID ) )
    sppoly_df[, pB$variabletomodel] = BM[[ paste(pB$variabletomodel, "predicted", sep=".") ]] [kk]

    sppoly_df[,  p$variabletomodel] = NA
    BM = NULL
    sppoly_df$StrataID = as.character( sppoly_df$StrataID )
    sppoly_df$tag ="predictions"

    vn = c( p$variabletomodel, pB$variabletomodel, "tag", "StrataID")

    M = rbind( M[, vn], sppoly_df[, vn] )
    sppoly_df = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    M$strata  = as.numeric( M$StrataID)
    M$zi = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$iid_error = 1:nrow(M) # for inla indexing for set level variation

    save( M, file=fn, compress=TRUE )
    return( M )
  }

}
