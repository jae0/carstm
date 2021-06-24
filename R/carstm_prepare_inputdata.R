carstm_prepare_inputdata = function( p, M, sppoly,
    lookup = c("bathymetry", "substrate", "temperature", "speciescomposition"),
    APS_data_offset=NULL
) {
 

  M$tag = "observations"

  if (grepl("space", p$aegis_dimensionality)) {

    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
 
    # observations
    M$AUID = st_points_in_polygons(
      pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname = "AUID"
    )
    M = M[!is.na(M$AUID),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    for (lu in lookup)  M = lookup_point_data( p=p, M=M, tolookup=lu )

    M$plon = M$plat = M$lon = M$lat = NULL

    if (any( grepl("offset", as.character(p$carstm_model_formula)))){
      if (!exists("data_offset", M)) {
        if (exists("data_offset", sppoly)) {
          message("data_offset not defined, using data_offset from sppoly")
          M$data_offset = sppoly$data_offset[ match(  M$AUID,  sppoly$AUID ) ]
        }
      }
    }
  }

  # in case it is purely spatial
  if ( grepl( "year", p$aegis_dimensionality ) ) {
    if (!exists("year", M)) {
      if (exists("yr", M)) names(M)[which(names(M)=="yr") ] = "year"
    }
  }

  if ( grepl( "season", p$aegis_dimensionality ) ) {
    if (!exists("dyear", M)) {
      if (exists("tiyr", M )) M$dyear = M$tiyr - M$year 
    }
  }


  # end observations
  # ----------


  # ----------
  # predicted locations (APS)

  if (grepl("space", p$aegis_dimensionality)) {

    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)

    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"

    if (any( grepl("offset", as.character(p$carstm_model_formula)))) {
      if (!is.null(APS_data_offset)) {
        APS$data_offset =  APS_data_offset   
      } else {
        message( "APS_data_offset is required as there is an offset in the formula ... using 1 for now")
        APS$data_offset = 1
      }
      APS = APS[ , c( "AUID", "tag", "data_offset" ) ]
    } else {
      APS = APS[ , c( "AUID", "tag" ) ]
    }

    APS[,p$variabletomodel] = NA
  }

  if ( "bathymetry" %in% lookup ) {
    require(aegis.bathymetry)
    pB = bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
    APS[, pB$variabletomodel] = bathymetry_lookup(  LOCS=sppoly,
      lookup_from = p$carstm_inputdata_model_source$bathymetry,
      lookup_to = "areal_units",
      vnames=pB$variabletomodel
    )
  }

  if ( "substrate" %in% lookup ) {
    require(aegis.substrate)
    pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
    APS[, pS$variabletomodel] = substrate_lookup(  LOCS=sppoly,
      lookup_from = p$carstm_inputdata_model_source$substrate,
      lookup_to = "areal_units",
      vnames=pS$variabletomodel
    )
  }


  # prediction surface in time
  # to this point APS is static, now add time dynamics (teperature)

  # in case it is purely spatial
  if ( grepl( "year", p$aegis_dimensionality ) | (grepl( "season", p$aegis_dimensionality )  ) ) {
    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS)[grep("rep.int", names(APS))] = "tiyr"

    APS$timestamp = lubridate::date_decimal( APS$tiyr, tz=p$timezone )

    if ( grepl( "year", p$aegis_dimensionality ) ) {
      APS$year = aegis_floor( APS$tiyr)
    }
    if ( grepl( "season", p$aegis_dimensionality ) ) {
      APS$dyear = APS$tiyr - APS$year
    }
  }

  vn  = names(APS) #update


  # ---------------------
  if ( "temperature" %in% lookup ) {
    require("aegis.temperature")
    pT = temperature_parameters( p=parameters_reset(p), project_class="carstm", year.assessment=p$year.assessment  )
    APS[, pT$variabletomodel] = temperature_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly,
      lookup_from = p$carstm_inputdata_model_source$temperature,
      lookup_to = "areal_units",
      vnames_from= paste( pT$variabletomodel, "predicted", sep="."),
      vnames=pT$variabletomodel ,
      year.assessment=p$year.assessment
    )
  }


  if ( "speciescomposition" %in% lookup ) {
    require("aegis.speciescomposition")

    pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" , year.assessment=p$year.assessment)

    APS[, pPC1$variabletomodel] = speciescomposition_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly,
      lookup_from = p$carstm_inputdata_model_source$speciescomposition,
      lookup_to = "areal_units",
      vnames_from= paste( pPC1$variabletomodel, "predicted", sep="."),
      vnames=pPC1$variabletomodel ,
      year.assessment=p$year.assessment
    )


    pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2", year.assessment=p$year.assessment )
    APS[, pPC2$variabletomodel] = speciescomposition_lookup(  LOCS=APS[ , c("AUID", "timestamp")], AU_target=sppoly,
      lookup_from = p$carstm_inputdata_model_source$speciescomposition,
      lookup_to = "areal_units",
      vnames_from= paste( pPC2$variabletomodel, "predicted", sep="."),
      vnames=pPC2$variabletomodel,
      year.assessment=p$year.assessment
    )
  }

  APS$timestamp = NULL  # time-based matching finished

  M = rbind( M[, names(APS)], APS )

  APS = NULL; gc()

  M$uid = 1:nrow(M)  # seems to require an iid model for each obs for stability .. use this for iid
  M$AUID  = as.character(M$AUID)  # revert to factors -- should always be a character
  M$space = as.character( M$AUID)

  if (exists("tiyr", M)) M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
  if (exists("tiyr", M)) M$yr = aegis_floor( M$tiyr)
  
  if (exists("yr", M))  M$time = as.character( M$yr )  # copy for INLA
  if (exists("yr", M))  M$yr_factor = factor(M$yr)
  
  if (exists("dyear", M)) M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )
  if (exists("dyri", M))  M$season = as.character( M$dyri )  # copy for INLA
   
  return(M)
}