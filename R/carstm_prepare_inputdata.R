carstm_prepare_inputdata = function( p, M, sppoly,
    lookup = c("bathymetry", "substrate", "temperature", "speciescomposition"),
    varstoretain = c("data_offset"),
    APS_data_offset=1
) {
 
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
 
    # observations
    M$AUID = st_points_in_polygons(
      pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname = "AUID"
    )
    M = M[!is.na(M$AUID),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    # in case it is purely spatial
    if (exists("yr", M)) {
      names(M)[which(names(M)=="yr") ] = "year"
      # M = M[ which(M$year %in% p$yrs), ]
      # M$tiyr = lubridate::decimal_date ( M$timestamp )
      # M$dyear = M$tiyr - M$year
    }


    for (lu in lookup)  M = lookup_point_data( p=p, M=M, tolookup=lu )

    M$plon = M$plat = M$lon = M$lat = NULL
    M = M[ which(!is.na(M$AUID)),]
    M$AUID = as.character( M$AUID )  # match each datum to an area
    M$tag = "observations"

    # end observations
    # ----------


    # ----------
    # predicted locations (APS)
    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)

    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA
    vn = p$variabletomodel


    if ( "bathymetry" %in% lookup ) {
      pB = bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
      APS[, pB$variabletomodel] = bathymetry_lookup(  LOCS=sppoly,
        lookup_from = p$carstm_inputdata_model_source$bathymetry,
        lookup_to = "areal_units",
        vnames="z"
      )
      vn = c(vn, pB$variabletomodel)
    }

    if ( "substrate" %in% lookup ) {
      pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
      APS[, pS$variabletomodel] = substrate_lookup(  LOCS=sppoly,
        lookup_from = p$carstm_inputdata_model_source$substrate,
        lookup_to = "areal_units",
        vnames="substrate.grainsize"
      )
      vn = c(vn, pS$variabletomodel)
    }

    vn = c( vn, "tag", "AUID" )

    # ---------------------
    if ( "temperature" %in% lookup | "speciescomposition" %in% lookup) {
      # to this point APS is static, now add time dynamics (teperature)
      APS = APS[, vn]

      # expand APS to all time slices
      n_aps = nrow(APS)
      APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
      names(APS) = c(vn, "tiyr")
      APS$year = aegis_floor( APS$tiyr)
      APS$dyear = APS$tiyr - APS$year
      APS$timestamp = lubridate::date_decimal( APS$tiyr, tz=p$timezone )

      if ( "temperature" %in% lookup ) {
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
    }

    # useful vars to have for analyses outside of carstm_summary
    for (vn in varstoretain) if (!exists( vn, APS)) APS[,vn] = NA
    APS$data_offset = APS_data_offset  # force to solve for unit area (1 km^2)

    M = rbind( M[, names(APS)], APS )

    APS = NULL; gc()

    M$AUID  = as.character(M$AUID)  # revert to factors -- should always be a character

    if (exists("tiyr", M)) {
      M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
      M$yr = aegis_floor( M$tiyr)
      M$dyear =  M$tiyr - M$yr   # revert dyear to non-discretized form
      M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )
      M$time = as.character( M$yr )  # copy for INLA
      M$season = as.character( M$dyri )  # copy for INLA
      M$yr_factor = factor(M$yr)
    }

    #required for carstm formulae
    M$space = as.character( M$AUID)
    M$uid = 1:nrow(M)  # seems to require an iid model for each obs for stability .. use this for iid

  return(M)
}