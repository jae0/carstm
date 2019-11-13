
temperature_carstm = function ( p=NULL, DS="parameters", redo=FALSE, ... ) {

  # over-ride default dependent variable name if it exists
   require( carstm)

    if ( is.null(p)) {
      p = aegis.temperature::temperature_parameters(...)
    } else {
      p = aegis.temperature::temperature_parameters(p=p, ...)
    }

  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.coastline","aegis.polygons", "aegis.temperature" , "carstm" ) )


  # OSD data series variables of interest

  # ----------------------

  if (DS=="parameters_override") {
    # translate param values from one project to a unified representation
    # must be first to catch p
    pc = temperature_carstm(
      DS = "parameters",
      project_name = "temperature",
      variabletomodel = "t",
      yrs = p$yrs,
      spatial_domain = p$spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_overlay = p$areal_units_overlay, # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = p$areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = p$areal_units_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = p$inputdata_temporal_discretization_yr,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      modeldir = p$modeldir,  # outputs all go the the main project's model output directory
      auid = p$auid
    )
    return(pc)
  }

  # ---------------------


  if (DS=="parameters") {

    p$libs = unique( c( p$libs, project.library ( "carstm" ) ) )

    if ( !exists("project_name", p)) p$project_name = "temperature"

    if ( !exists("areal_units_strata_type", p)) p$areal_units_strata_type = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( p$spatial_domain == "SSE" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "groundfish_strata" #.. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      p$inputdata_temporal_discretization_yr = 1/12  #  24/365,  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( p$spatial_domain == "snowcrab" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      p$inputdata_temporal_discretization_yr = 1/52  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling }
    }

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla.default"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )

        inla_nthreads = ifelse( exists("inla_nthreads", p ), p$inla_nthreads, 1 )
        inla_nthreads_blas = ifelse ( exists("inla_nthreads_blas", p ), p$inla_nthreads_blas, 1 )

        p$carstm_model_label = "production"
        p$carstm_modelcall = paste('
          inla(
            formula = ', p$variabletomodel, ' ~ 1
              + f(year_factor, model="ar1", hyper=H$ar1 )
              + f(dyri, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2 )
              + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
              + f(strata, model="bym2", graph=sppoly@nb, group= year_factor,  scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "normal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list(strategy="gaussian", int.strategy="eb") ,# to get empirical Bayes results much faster.
            # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
            control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            num.threads=', inla_nthreads, ' ,
            blas.num.threads=', inla_nthreads_blas, ' ,
            verbose=TRUE
          ) ' )
      }
         #     + f(tiyr2, model="seasonal", season.length=10 )
        # + f(dyear, model="ar1", hyper=H$ar1 )

      if ( grepl("glm", p$carstm_modelengine) ) {
        p$carstm_modelcall = paste( 'glm( formula = ', p$variabletomodel, ' ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="identity")  ) ' )  # for modelengine='glm'
      }
      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv"  ) ) )
        p$carstm_modelcall = paste( 'gam( formula = ', p$variabletomodel, ' ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="identity")  ) ' ) # for modelengine='gam'
      }
    }

    p = carstm_parameters( p=p, DS="basic" )  # fill in anything missing and some checks

    return(p)
  }


  # ---------------------------------


  if ( DS=="carstm_inputs") {


    aggregate_data = FALSE
    if (exists("carstm_inputs_aggregated", p)) {
      if (p$carstm_inputs_aggregated)  aggregate_data = TRUE
    }

    if (aggregate_data) {
      fn = file.path( p$modeldir, paste( "temperature", "carstm_inputs", p$auid,
        p$inputdata_spatial_discretization_planar_km,
        round(p$inputdata_temporal_discretization_yr, 6),
        "rdata", sep=".") )

    } else {
      fn = file.path( p$modeldir, paste( "temperature", "carstm_inputs", p$auid,
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

      M = temperature.db( p=p, DS="aggregated_data"  )
      names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel

    } else {
      M = temperature.db( p=p, DS="bottom.all"  )
      names(M)[which(names(M)=="t")] = p$variabletomodel
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

      names(M)[which(names(M)=="yr") ] = "year"
      M = M[ which(M$year %in% p$yrs), ]
      M$tiyr = lubridate::decimal_date ( M$date )

      # globally remove all unrealistic data
      keep = which( M[,p$variabletomodel] >= -3 & M[,p$variabletomodel] <= 25 ) # hard limits
      if (length(keep) > 0 ) M = M[ keep, ]

      # p$quantile_bounds_data = c(0.0005, 0.9995)
      if (exists("quantile_bounds_data", p)) {
        TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds_data, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
        keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
        # this was -1.7, 21.8 in 2015
      }

      keep = which( M$z >=  2 ) # ignore very shallow areas ..
      if (length(keep) > 0 ) M = M[ keep, ]

      # in case plon/plats are from an alternate projection  .. as there are multiple data sources
      M = lonlat2planar( M, p$aegis_proj4string_planar_km)

      M$dyear = M$tiyr - M$year

    }

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area

    pB = bathymetry_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry

    # pS = substrate_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry
    # M[, pS$variabletomodel] = lookup_substrate_from_surveys(  p=pS, locs=M[, c("lon", "lat")] )

    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area

    M$tiyr = M$year + M$dyear
    M$tag = "observations"

    # already has depth .. but in case some are missing data
    pB = bathymetry_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry

    if (!(exists(pB$variabletomodel, M ))) M[,pB$variabletomodel] = NA

    kk = which(!is.finite( M[, pB$variabletomodel] ))
    if (length(kk > 0)) M[kk, pB$variabletomodel] = lookup_bathymetry_from_surveys( p=pB, locs=M[kk, c("lon", "lat")] )

    # if any still missing then use a mean depth by StrataID
    kk =  which( !is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0) {
      AD = bathymetry.db ( p=pB, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$StrataID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      oo = tapply( AD[, paste(pB$variabletomodel, "mean", sep="." )], AD$StrataID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$StrataID[kk]), as.character( names(oo )) )
      M[kk, pB$variabletomodel] = oo[jj ]
    }

    APS = as.data.frame(sppoly)
    APS$StrataID = as.character( APS$StrataID )
    APS$tag ="predictions"
    APS[, p$variabletomodel] = NA

    pB$carstm_model_label = "production"
    BM = carstm_model ( p=pB, DS="carstm_modelled" )
    jj = match( as.character( APS$StrataID), as.character( BM$StrataID) )
    APS[, pB$variabletomodel] = BM[[ paste(pB$variabletomodel, "predicted", sep=".") ]] [jj]
    jj =NULL
    BM = NULL

    vn = c( p$variabletomodel, pB$variabletomodel, "tag", "StrataID"  )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")

    M = rbind( M[, names(APS)], APS )
    APS = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    M$strata  = as.numeric( M$StrataID)

    M$year = trunc( M$tiyr)
    M$year_factor = as.numeric( factor( M$year, levels=p$yrs))
    M$dyear =  M$tiyr - M$year  # reset in case it has been discretized
    # M$tiyri  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla


    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)

    M$zi = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )

    save( M, file=fn, compress=TRUE )
    return( M )
  }


  # -----------------------


}
