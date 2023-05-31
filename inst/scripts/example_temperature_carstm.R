
# carstm example using temperature data subset


# prep input data (copied from stmv): 
# centered over Halifax, NS: bottemp[lon>-65 & lon< -62 & lat <45 &lat>43,]
# t=temperature (C); z=depth (m); tiyr=decimal year

fn = file.path( project.codedirectory("carstm", "inst", "extdata"), "aegis_spacetime_test.RDS")
bottemp = readRDS(fn)   
# bottemp = bottemp[lon>-65 & lon< -62 & lat <45 &lat>43,]

plot(lat ~ -lon, bottemp)
str(bottemp)
summary(bottemp)

plot(lon~lat, bottemp)
hist( bottemp$tiyr )  # decimal date

# required parameter settings:
p = list()
p$yrs = min(year(bottemp$date)):max(year(bottemp$date)) # 1980:2010

p$year.assessment = max(p$yrs)
 
# create/update library list
standard_libs = c( "colorspace", "lubridate",  "lattice", 
    "parallel", "sf", "GADMTools", "INLA" , "data.table" )

local_libs = c("aegis", "aegis.bathymetry", "aegis.coastline", 
    "aegis.polygons", "aegis.substrate", "aegis.temperature", "aegis.survey" ) 
    
p$libs = RLibrary ( c(standard_libs, local_libs) )

p$project_name = "test_ocean_bottom_temperatures_halifax"
p$data_root = file.path( "~", "test", p$project_name ) 
p$datadir  = file.path( p$data_root, "data" )
p$modeldir = file.path( p$data_root, "modelled" ) 

if ( !file.exists(p$data_root) ) dir.create( p$data_root, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )

 
p$variabletomodel = "t"
p$dimensionality="space-time-cyclic"
p$quantile_bounds =c(0.005, 0.995) # trim upper bounds (in posterior predictions)

# space resolution
p$aegis_proj4string_planar_km = projection_proj4string("utm20")
p$dres = 1/60/4 # resolution in angular units (degrees)
p$pres = 1  # spatial resolution in planar units (km)
p$lon0 = min( bottemp$lon )
p$lon1 = max( bottemp$lon )
p$lat0 = min( bottemp$lat )
p$lat1 = max( bottemp$lat )
p$psignif = 1  

p$nlons = trunc( diff(range(c(p$lon0,p$lon1)))/p$dres) + 1L
p$nlats = trunc( diff(range(c(p$lat0,p$lat1)))/p$dres) + 1L
corners = data.frame(lon=c(p$lon0,p$lon1), lat=c(p$lat0,p$lat1))
corners = lonlat2planar( corners, proj.type=p$aegis_proj4string_planar_km )
corners$plon = round( corners$plon, p$psignif)  # this matches the p$pres value of x km resolution
corners$plat = round( corners$plat, p$psignif)  # this matches the p$pres value of x km resolution
p$corners=corners

p$plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
p$plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
p$nplons = length(plons)
p$nplats = length(plats)
p$origin = c(min(p$corners$plon), min(p$corners$plat ))
p$gridparams = list( dims=c(p$nplons, p$nplats), origin=p$origin, res=c(p$pres, p$pres) ) # used for fast indexing and merging


p$year.assessment = max(p$yrs)
p$timezone="America/Halifax" 

# time resolution
p$ny = length(p$yrs)
p$nw = 10 # default value of 10 time steps number of intervals in time within a year for all temp and indicators
p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components
p$dyears = (c(1:p$nw)-1) / p$nw # intervals of decimal years... fractional year breaks
p$dyear_centre = p$dyears[ trunc(p$nw/2) ] + p$tres/2
p$prediction_dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01")) # used for creating timeslices and predictions  .. needs to match the values in aegis_parameters()
p$nt = p$nw*p$ny # i.e., seasonal with p$nw (default is annual: nt=ny)

# predictions at these time values (decimal-year), # output timeslices for predictions in decimla years, yes all of them here
tout = expand.grid( yr=p$yrs, dyear=1:p$nw, KEEP.OUT.ATTRS=FALSE )
p$prediction_ts = sort( tout$yr + tout$dyear/p$nw - p$tres/2 )# mid-points


p$inputdata_spatial_discretization_planar_km = p$pres / 10 # controls resolution of data prior to modelling (km )
p$inputdata_temporal_discretization_yr = 1/52  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling;; use 1/12 -- monthly or even 1/4.. if data density is low
p$dyear_discretization_rawdata = c( {c(1:365)-1}/365, 1)  # dyear_discretization_rawdata :: intervals of decimal years... fractional year breaks finer than the default 10 units (taking daily for now..) .. need to close right side for "cut" .. controls resolution of data prior to modelling


# areal units information
p$spatial_domain = "halifax"

p$areal_units_proj4string_planar_km =  p$aegis_proj4string_planar_km   # coord system to use for areal estimation and gridding for carstm
p$areal_units_type= "tesselation" 
p$areal_units_constraint_ntarget = length(p$yrs)   # n time slices req in each au
p$areal_units_constraint_nmin = 5    # n time slices req in each au
p$areal_units_resolution_km = 1   # starting resolution .. if using tesselation/ otherwise grid size ()
p$areal_units_overlay = "none"  
p$areal_units_timeperiod = "none"   # only relevent for groundfish polys

p$tus="yr" 
p$hull_alpha = 20 
p$fraction_todrop = 0.05 
p$fraction_cv = 1.0   # approx poisson (binomial)
p$fraction_good_bad = 0.9 
p$nAU_min = 30 


# carstm-specific parameters
p$project_class = "carstm"
p$carstm_model_label = "test_basic_form"
p$carstm_modelengine = "inla"   # {model engine}.{label to use to store}
p$carstm_inputs_prefilter = "aggregated" 
p$carstm_inputs_prefilter_n = 100  # only used for "sampled"


# create polygon  :: requires aegis, aegis.coastline, aegis.polygons


sppoly = try( areal_units( p=p, areal_units_directory=p$datadir, redo=FALSE ) )

if (inherits("try-error", sppoly)) {
    sppoly = areal_units( 
    p=p, 
    xydata=setDF(bottemp[,.(lon, lat, yr=floor(tiyr) ) ]), 
    areal_units_directory=p$datadir, 
    redo=TRUE )  # to force create
)

plot( sppoly[ "AUID" ] ) 
carstm_map( sppoly=sppoly, vn="au_sa_km2", map_mode="view" )  # interactive
 

crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
sppoly = st_transform(sppoly, st_crs(crs_lonlat))

bottemppts = st_as_sf( bottemp[,c("lon","lat")], coords=c("lon","lat"), crs=crs_lonlat )

# observations
bottemp$tag ="observations"
bottemp$time = year(bottemp$date)
bottemp$AUID = st_points_in_polygons(
    pts = bottemppts,
    polys = sppoly[, "AUID"],
    varname = "AUID"
)

depths = bottemp[ , .(z=median(z, na.rm=TRUE)), by=.(AUID) ]

APS = st_drop_geometry(sppoly)
setDT(APS)

APS$AUID = as.character( APS$AUID )
APS$tag ="predictions"
APS[, p$variabletomodel] = NA

APS = APS[ depths, on=.(AUID) ] 
   
n_aps = nrow(APS)
APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
names(APS)[ncol(APS)] = "tiyr"
APS$timestamp = lubridate::date_decimal( APS$tiyr, tz=p$timezone )
APS$time = trunc( APS$tiyr)  # year ("time*" is a keyword)
APS$dyear = APS$tiyr - APS$time 
 

vvv = intersect( names(APS), names(bottemp) )
M = rbind( bottemp[, vvv, with=FALSE ], APS[, vvv, with=FALSE ] )

APS = NULL; gc()

# M$uid = 1:nrow(M)  # seems to require an iid model for each obs for stability .. use this for iid
M$AUID  = as.character(M$AUID)  # revert to factors -- should always be a character

M$space = match( M$AUID, sppoly$AUID) 
M$space_time = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 

M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
M$time = trunc( M$tiyr)
M$time_space = M$time  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 

M$dyear = M$tiyr - M$time 
M$tiyr = NULL

# do not sepraate out as season can be used even if not predicted upon
ii = which( M$dyear > 1) 
if (length(ii) > 0) M$dyear[ii] = 0.99 # cap it .. some surveys go into the next year
 
M$dyri = discretize_data( M[["dyear"]], discretizations()[["dyear"]] )

cyclic_levels = factor(p$dyears + diff(p$dyears)[1]/2, ordered=TRUE )
M$cyclic = factor( as.character( M$dyri ), levels =levels(cyclic_levels) )   # copy for carstm/INLA


# "H" in formula are created on the fly in carstm ... they can be dropped in formula or better priors defined manually

formula = as.formula( paste(
    p$variabletomodel, ' ~ 1',
    ' + f( time, model="ar1",  hyper=H$ar1 ) ',   
    ' + f( cyclic, model="seasonal", scale.model=TRUE, season.length=10, hyper=H$iid  ) ',
    ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, hyper=H$bym2  ) ',
    ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
    ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group) ) '
) )


family = "gaussian"
 
require(carstm)
loadfunctions("carstm")

# takes about 15 minutes
res = carstm_model( 
    p=p, 
    sppoly=sppoly,
    posterior_simulations_to_retain=c("predictions", "random_spatial"), 
    nposteriors=1000,  # 1000 to 5000 would be sufficient to sufficiently sample most distributions: trade-off between file size and information content
    dimensionality="space-time-cyclic",
    # redo_fit=FALSE,  # if FALSE then reload fit and recompute posteriors 
    # redo_fit=TRUE,  # if TRUE then compute fit and compute posteriors 
    # args below are INLA options, passed directly to INLA
    formula=formula,
    family=family,
    data =M,  
    num.threads="6:2",  # adjust for your machine
    mc.cores=2,
    control.inla = list( strategy='laplace'  ),
    verbose=TRUE 
)    

# to load saved fit
# can be very large files .. slow 
fit = carstm_model( p=p, sppoly=sppoly, DS="carstm_modelled_fit")

summary(fit)
names(fit)
fit$summary$dic$dic
fit$summary$dic$p.eff

plot(fit)
plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )



# to load saved results summary
res = carstm_model( p=p, sppoly=sppoly, DS="carstm_modelled_summary")
( res$summary)

b0 = res$summary$fixed_effects["(Intercept)", "mean"]

ts =  res$random$time 
plot( mean ~ ID, ts, type="b", ylim=c(-2,2), lwd=1.5, xlab="year")
lines( quant0.025 ~ ID, ts, col="gray", lty="dashed")
lines( quant0.975 ~ ID, ts, col="gray", lty="dashed")


ts =  res$random$cyclic
plot( mean ~ID, ts, type="b", ylim=c(-1.5, 1.5), lwd=1.5, xlab="fractional year")
lines( quant0.025 ~ID, ts, col="gray", lty="dashed")
lines( quant0.975 ~ID, ts, col="gray", lty="dashed")



    

map_centre = c( (p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2   )
map_zoom = 7

# maps of some of the results
tmatch="2010"
umatch="0.75"  # == 0.75*12 = 9 (ie. Sept)  

tmout = carstm_map(  res=res, vn="predictions", tmatch=tmatch, umatch=umatch, 
    sppoly=sppoly,
    breaks=seq(-1, 9, by=1), 
    palette="-RdYlBu",
    plot_elements=c( "isobaths",  "compass", "scale_bar", "legend" ),
    tmap_zoom= c(map_centre, map_zoom),
    title=paste( "Bottom temperature predictions", tmatch, umatch)  
)
tmout

# persistent spatial effects
tmout = carstm_map(  res=res, vn=c( "random", "space", "combined" ), 
    sppoly=sppoly,
    breaks=seq(-5, 5, by=1), 
    palette="-RdYlBu",
    plot_elements=c( "isobaths",  "compass", "scale_bar", "legend" ),
    tmap_zoom= c(map_centre, map_zoom),
    title="Bottom temperature spatial effects (Celsius)"
)
tmout



# finished
