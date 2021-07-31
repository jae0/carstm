
carstm_prepare_inputdata = function( p, M, sppoly,
    lookup = c("bathymetry", "substrate", "temperature", "speciescomposition"),
    APS_data_offset=NULL, NA_remove=TRUE
) {
 
  setDT(M)

  M$tag = "observations"
 
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))

  # observations
  M$AUID = st_points_in_polygons(
    pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
    polys = sppoly[, "AUID"],
    varname = "AUID"
  )
  M = M[!is.na(M$AUID),]
  M$AUID = as.character( M$AUID )  # match each datum to an area

      
  if ("bathymetry" %in% lookup) {
    message( "lookup: bathymetry")
    require(aegis.bathymetry)
    pB = bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
    if ( !(exists( pB$variabletomodel, M )))  M[[pB$variabletomodel]] = NA
    iM = which(!is.finite( M[[pB$variabletomodel]] ))
    if (length(iM > 0)) {
      M[[pB$variabletomodel]][iM] = aegis_lookup(  
        data_class="bathymetry", 
        LOCS=M[ iM, c("lon", "lat")],  
        project_class="core", 
        output_format="points" , 
        DS="aggregated_data", 
        variable_name="z.mean", 
        returntype="vector" 
      ) 
    }
 
    if ( exists("spatial_domain", p)) {
        # need to be careful with extrapolation ...  filter depths
        if (NA_remove)  M = M[ is.finite(M[[pB$variabletomodel]] ) , ]
        ii = geo_subset( spatial_domain=p$spatial_domain, Z=M )
        if (length(ii)> 0 ) M = M[ ii , ] 
    }

    if ( p$carstm_inputdata_model_source$bathymetry %in% c("stmv", "hybrid") ) {
      pBD = bathymetry_parameters(  spatial_domain=p$spatial_domain, project_class=p$carstm_inputdata_model_source$bathymetry )  # full default
      LU = bathymetry_db( p=pBD, DS="baseline", varnames="all" )
      iML = match( 
        array_map( "xy->1", M[, c("plon","plat")], gridparams=p$gridparams ), 
        array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams ) 
      )
      vns = intersect(  c( "z", "dZ", "ddZ", "b.sdSpatial", "b.sdObs", "b.phi", "b.nu", "b.localrange" ), names(LU) )
      for (vn in setdiff( vns, "z") ) M[[ vn]] = LU[ iML, vn ]
      LU =  iML = vns = NULL
      if (NA_remove) {
        ii = which( is.finite( rowSums(M[, vns, with=FALSE] )  ))
        if (length(ii) > 0 ) M = M[ ii, ]
      }
    }

  }


    # --------------------------


    if ("substrate" %in% lookup) {
      message( "lookup: substrate")
      require(aegis.substrate)
      pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
      if (!(exists(pS$variabletomodel, M ))) M[[pS$variabletomodel]] = NA
      iM = which(!is.finite( M[[pS$variabletomodel]] ))
      if (length(iM > 0)) {
        M[[pS$variabletomodel]] [iM]  = aegis_lookup( 
          data_class="substrate", 
          LOCS=M[iM, c("lon", "lat")], 
          project_class="core", 
          output_format="points", 
          DS="aggregated_data", 
          variable_name="substrate.grainsize.mean" 
        )  
      }   
      

      if ( p$carstm_inputdata_model_source$substrate %in% c("stmv", "hybrid") ) {
        pBD = bathymetry_parameters(  spatial_domain=p$spatial_domain, project_class=p$carstm_inputdata_model_source$substrate  )  # yes substrate source coordinate system.. to match substrate source for the data
        LUB = bathymetry_db( p=pBD, DS="baseline", varnames="all" )
        iML = match( 
          array_map( "xy->1", M[, c("plon","plat")], gridparams=p$gridparams ), 
          array_map( "xy->1", LUB[,c("plon","plat")], gridparams=p$gridparams ) 
        )

        pSD = substrate_parameters(  spatial_domain=p$spatial_domain, project_class=p$carstm_inputdata_model_source$substrate )  # full default
        LU = substrate_db( p=pSD, DS="complete"  )

        vns = intersect(  c( 
            "substrate.grainsize", "substrate.grainsize.lb", "substrate.grainsize.ub", 
            "s.sdTotal", "s.sdSpatial", "s.sdObs", "s.phi", "s.nu", "s.localrange" 
          ), names(LU) )

        for (vn in vns  ) M[[ vn]] = LU[ iML, vn ]
        
        if (NA_remove) {
          ii = which( is.finite( rowSums(M[, vns, with=FALSE] )  ))
          if (length(ii) > 0 ) M = M[ ii, ]
        }
        LU =  iML = vns = NULL
      
      }
 
    }


    # --------------------------


    if ("temperature" %in% lookup) {
      require(aegis.temperature)
      pT = temperature_parameters( p=parameters_reset(p), project_class="carstm", year.assessment=p$year.assessment  )
      if (!(exists(pT$variabletomodel, M ))) M[[pT$variabletomodel]] = NA
      iM = which(!is.finite( M[[pT$variabletomodel]] ))
      if (length(iM > 0)) {
        M[[pT$variabletomodel]] [iM] = aegis_lookup(  
          data_class="temperature", 
          LOCS=M[ iM, c("lon", "lat", "timestamp")],
          project_class="core", 
          DS="aggregated_data", 
          output_format="points", 
          variable_name="t.mean", 
          tz="America/Halifax",
          year.assessment=p$year.assessment
        )
      }
      if (NA_remove) M = M[ is.finite(M[[ pT$variabletomodel]]  ) , ]
      M = M[ which( M[[ pT$variabletomodel]]  < 14 ) , ]  #

      # to to:  add st,v/hybrid 
    }


    # --------------------------


    if ("speciescomposition" %in% lookup) {
      require(aegis.temperature)

      pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" , year.assessment=p$year.assessment)
      if (!(exists(pPC1$variabletomodel, M ))) M[[pPC1$variabletomodel]] = NA
      iM = which(!is.finite( M[[pPC1$variabletomodel]] ))
      if (length(iM > 0)) {
        M[[pPC1$variabletomodel]][iM] = aegis_lookup(  
          data_class="speciescomposition", 
          LOCS=M[ iM, c("lon", "lat", "timestamp")],
          project_class="core", 
          DS="speciescomposition", 
          output_format="points", 
          variable_name="pca1", 
          tz="America/Halifax" ,
          year.assessment=p$year.assessment
        )
      }
      if (NA_remove) M = M[ which(is.finite(M[[pPC1$variabletomodel]] )), ]

      pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2", year.assessment=p$year.assessment )
      if (!(exists(pPC2$variabletomodel, M ))) M[,pPC2$variabletomodel] = NA
      iM = which(!is.finite( M[[pPC2$variabletomodel]] ))
      if (length(iM > 0)) {
         M[[pPC2$variabletomodel]][iM]  = aegis_lookup(  
           data_class="speciescomposition", 
           LOCS=M[ iM, c("lon", "lat", "timestamp")], 
           project_class="core", 
           DS="speciescomposition", 
           output_format="points", 
           variable_name="pca2", 
           tz="America/Halifax" ,
           year.assessment=p$year.assessment
         )
      }
      if (NA_remove) M = M[ which(is.finite(M[[pPC2$variabletomodel]] )),]


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
      

    if ( grepl( "year", p$aegis_dimensionality ) ) {
      if (!exists("year", M)) {
        if (exists("yr", M)) names(M)[which(names(M)=="yr") ] = "year"
      }
      # though dyear is required for seasonal models, it can also be used in annual as well 
      # that is, if not predicting the full seasonal array but only A GIVEN time slice 
      if (!exists("dyear", M)) {
        if (exists("tiyr", M )) M$dyear = M$tiyr - M$year 
      }
    }

    # to to:  add st,v/hybrid 


  # end observations
  # ----------


  # ----------
  # generate prediction surface locations (APS) .. use carstm predictions (project_class)
  if (grepl("space", p$aegis_dimensionality)) {

    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)
    setDT(APS)

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
    APS[[pB$variabletomodel]] = aegis_lookup( 
      data_class="bathymetry", 
      LOCS=sppoly$AUID,
      LOCS_AU=sppoly,
      project_class = "carstm", # lookup from modelled predictions from carstm
      output_format = "areal_units",
      variable_name= list("predictions"),
      statvars=c("mean"),
      raster_resolution=min(p$gridparams$res) /2,
      returntype = "vector"
    ) 

    iM = which(!is.finite( APS[[pB$variabletomodel]] )) 
    if (length(iM) > 0 ) {
      # depth is very important
      APS[[pB$variabletomodel]][iM] = aegis_lookup(  
        data_class="bathymetry", 
        LOCS=APS$AUID,
        LOCS_AU=sppoly,
        project_class = "stmv", # lookup from modelled predictions from carstm
        output_format = "areal_units",
        variable_name="z", 
        raster_resolution=min(p$gridparams$res) /2,
        returntype = "vector"
      ) 
 
    }
    iM = which(!is.finite( APS[[pB$variabletomodel]] )) 
    if (length(iM) > 0 ) {
      # depth is very important
      APS[[pB$variabletomodel]][iM] = aegis_lookup(  
        data_class="bathymetry", 
        LOCS=APS$AUID,
        LOCS_AU=sppoly,
        project_class = "core", # lookup from modelled predictions from carstm
        output_format = "areal_units",
        DS = "aggregated_data",  # needed for core 
        variable_name = "z.mean", 
        raster_resolution=min(p$gridparams$res) /2,
        returntype = "vector"
      ) 
 
    }

  }

  if ( "substrate" %in% lookup ) {
    require(aegis.substrate)
    pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
    APS[[pS$variabletomodel]]  = aegis_lookup( 
      data_class="substrate", 
      LOCS=sppoly$AUID,
      LOCS_AU=sppoly,
      project_class = "carstm", # lookup from modelled predictions from carstm
      output_format = "areal_units",
      variable_name = list("predictions"),
      statvars = c("mean"),
      raster_resolution = min(p$gridparams$res) /2,
      returntype = "vector"
    )  

    
    iM = which(!is.finite( APS[[pS$variabletomodel]] )) 
    if (length(iM) > 0 ) {

      APS[[pS$variabletomodel]][iM] = aegis_lookup(  
        data_class="substrate", 
        LOCS=APS$AUID,
        LOCS_AU=sppoly,
        project_class = "stmv", # lookup from modelled predictions from carstm
        output_format = "areal_units",
        variable_name = "substrate.grainsize", 
        raster_resolution=min(p$gridparams$res) /2,
        returntype = "vector"
      ) 

    }
  }
  
  # prediction surface in time
  # to this point APS is static, now add time dynamics (teperature),  expand APS to all time slices
  if ( grepl( "year", p$aegis_dimensionality ) | (grepl( "season", p$aegis_dimensionality )  ) ) {
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS)[ncol(APS)] = "tiyr"
    APS$timestamp = lubridate::date_decimal( APS$tiyr, tz=p$timezone )
    APS$year = aegis_floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$year
  }
  vn  = names(APS) #update


  # ---------------------
  if ( "temperature" %in% lookup ) {
    require("aegis.temperature")
    pT = temperature_parameters( p=parameters_reset(p), project_class="carstm", year.assessment=p$year.assessment  )
    APS[[ pT$variabletomodel ]] = aegis_lookup( 
      data_class="temperature", 
      LOCS=APS[ , c("AUID", "timestamp")], 
      LOCS_AU=sppoly,
      project_class = "carstm", # lookup from modelled predictions from carstm
      output_format = "areal_units",
      variable_name=list("predictions"),
      statvars=c("mean"),
      raster_resolution=min(p$gridparams$res) /2,
      year.assessment=p$year.assessment,
      returntype = "vector"
    )
  }


  if ( "speciescomposition" %in% lookup ) {
    require("aegis.speciescomposition")
    pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" , year.assessment=p$year.assessment)
    APS[[ pPC1$variabletomodel ]] = aegis_lookup( 
      data_class="speciescomposition", 
      LOCS=APS[ , c("AUID", "timestamp")], 
      LOCS_AU=sppoly, 
      project_class = "carstm", # lookup from modelled predictions from carstm 
      output_format = "areal_units",
      variable_name=list("predictions"),
      variabletomodel="pca1" ,
      statvars=c("mean"),
      raster_resolution=min(p$gridparams$res) /2,
      year.assessment=p$year.assessment,
      returntype = "vector"
    ) 


    pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2", year.assessment=p$year.assessment )
    APS[[ pPC2$variabletomodel ]] = aegis_lookup( 
      data_class="speciescomposition", 
      LOCS=APS[ , c("AUID", "timestamp")], 
      LOCS_AU=sppoly,
      project_class = "carstm", # lookup from modelled predictions from carstm
      output_format = "areal_units",
      variable_name=list("predictions"),
      variabletomodel="pca2" ,
      statvars=c("mean"),
      raster_resolution=min(p$gridparams$res) /2,
      year.assessment=p$year.assessment,
      returntype = "vector"
    ) 
 }

  if (exists("timestamp", APS)) APS$timestamp = NULL  # time-based matching finished (if any)

  M = rbind( M[, names(APS), with=FALSE ], APS )

  APS = NULL; gc()

  # M$uid = 1:nrow(M)  # seems to require an iid model for each obs for stability .. use this for iid
  M$AUID  = as.character(M$AUID)  # revert to factors -- should always be a character
  M$space = as.character( M$AUID)
 
  if (exists("tiyr", M)) {
    M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
    M$yr = aegis_floor( M$tiyr)
    M$time = as.character( M$yr )  # copy for INLA
    M$yr_factor = factor(M$yr)
  
    # do not sepraate out as season can be used even if not predicted upon
    ii = which( M$dyear > 1) 
    if (length(ii) > 0) M$dyear[ii] = 0.99 # cap it .. some surveys go into the next year

    M$dyri = discretize_data( M[["dyear"]], p$discretization[["dyear"]] )


    M$season = as.character( M$dyri )  # copy for INLA
  }


  if ( "substrate" %in% lookup ) {
    iM = which(!is.finite( APS[[pS$variabletomodel]] )) 
    if (length(iM) > 0 ) {
      
    }

  }

  return(M)
}
