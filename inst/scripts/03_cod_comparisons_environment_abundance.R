
# ------------------------------------------------
# Atlantic cod comparison .. adding environmental variation


# ------------------------------------------------
# load data common environment and parameter setting
# source( system.file( "scripts", "00_cod_comparisons_data_environment.R", package = "carstm") )

# --------------------------------
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
# NOTE: the data selection is the same as in (01_cod_comparisons_basic_stranal.R)
p = carstm::carstm_parameters(
  id ="Atlantic cod summer standardtow",
  speciesname = "Atlantic_cod",
  groundfish_species_code = 10,   #  10= cod
  yrs = 1970:2017,
  trawlable_units = "towdistance"  # <<<<<<<<<<<<<<<<<<
  # trawlable_units = "standardtow"
  # trawlable_units = "sweptarea"
)



# --------------------------------
# parameter setting used to filter data via 'survey.db( DS="filter")'
# unlike stratanl, we do not need to remove strata until the last /aggregation step
p$selection=list(
  biologicals=list(
    spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code )
  ),
  survey=list(
    data.source="groundfish",
    yr = p$yrs,  # time frame for comparison specified above
    months=6:8,  # "summer"
    # dyear = c(150,250)/365, # alternate way of specifying season: summer = which( (x>150) & (x<250) ) , spring = which(  x<149 ), winter = which(  x>251 )
    settype = 1, # same as geartype in groundfish db
    gear = c("Western IIA trawl", "Yankee #36 otter trawl"),
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
)



# ------------------------------------------------
## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
# Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
# This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)

sppoly = areal_units( strata_type="stratanal_polygons",  proj4string_planar_km=p$internal.crs_planar, timeperiod="pre2014" )
sppoly$strata_to_keep = ifelse( as.character(sppoly$StrataID) %in% strata_definitions( c("Gulf", "Georges_Bank", "Spring", "Deep_Water") ), FALSE,  TRUE )


# ------------------------------------------------
# neighbourhood structure --- required to do areal unit spatial modelling

W.nb = neighbourhood_structure( sppoly=sppoly, strata_type="stratanal_polygons" )


# --------------------------------
# Get the data
p$selection$survey$strata_toremove = NULL  # emphasize that all data enters analysis initially ..

set = survey.db( p=p, DS="filter" )

# categorize Strata 
o = over( SpatialPoints( set[,c("plon", "plat")], sp::CRS(p$internal.crs) ), spTransform(sppoly, sp::CRS(p$internal.crs) ) ) # match each datum to an area
set$StrataID = o$StrataID
o = NULL
set = set[ which(!is.na(set$StrataID)),]

set$totno[which(!is.finite(set$totno))] = NA


# --------------------------------
# ensure we have some estimate of sweptarea and choose the appropriate
# one based upon which trawlable units we are using
ft2m = 0.3048
m2km = 1/1000
nmi2mi = 1.1507794
mi2ft = 5280
standardtow_sakm2 = (41 * ft2m * m2km ) * ( 1.75 * nmi2mi * mi2ft * ft2m * m2km )  # surface area sampled by a standard tow in km^2  1.75 nm
set$data_offset = switch( p$trawlable_units,
  standardtow =  rep(standardtow_sakm2, nrow(set)) , # "standard tow"
  towdistance = set$sa_towdistance,  # "sa"=computed from tow distance and standard width, 0.011801==),
  sweptarea = set$sa  # swept area based upon stand tow width and variable lenths based upon start-end locations wherever possible
)
set$data_offset[which(!is.finite(set$data_offset))] = median(set$data_offset, na.rm=TRUE )  # just in case missing data


# ------------------------------------------------
# update set with StrataID factor variables and a few other repeatedly used variables
set$StrataID = factor(set$StrataID, levels=levels(sppoly$StrataID))
set$yr_factor = factor(set$yr)
set$iid_error = 1:nrow(set) # for inla indexing
set$tag = "observations"


## --------------------------------
# construct meanweights matrix
weight_year = meanweights_by_strata( set=set, StrataID=as.character( sppoly$StrataID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
# weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
# weight_year = weight_year[ match(as.character(sppoly$StrataID), rownames(weight_year) )]



# adjust based upon RAM requirements and ncores
ncores = floor( aegis.env::ram_local( "ncores", ram_main=4, ram_process=6 ) / 2 )
inla.setOption(num.threads=ncores)
inla.setOption(blas.num.threads=ncores)



# RES = data.frame(yr=p$selection$survey[["yr"]]) # collect model comparisons
if (0) {
  fn = file.path( getwd(), "RES.rdata" )
  # save(RES, file=fn)
  # load(fn)
}



## ----------------------------------
# covariates of interest
covars = c("t", "tsd", "tmax", "tmin", "degreedays", "z",  "dZ", "ddZ" )

  # currently supported:
  # z = depth (m)
  # dZ = bottom slope (m/km)
  # ddZ = bottom curvature (m/km^2)
  # substrate.grainsize = mean grain size of bottom substrate (mm)
  # t = temperature (C) – subannual
  # tlb = temperature lower 95% bound (C) –subannual
  # tub = temperature upper 95% bound (C) –subannual
  # tmean = mean annual temperature
  # tsd = standard deviation of the mean annual temperature
  # tmin = minimum value of temperature in a given year – annual
  # tmax = maximum value of temperature in a given year – annual
  # tamplitude = amplitude of temperature swings in a year (tmax-tmin) – annual
  # degreedays = number of degree days in a given year – annual


# extract covariate means by strata
res = aegis_db_extract_by_polygon(
  sppoly=sppoly,
  vars=covars,
  spatial.domain=p$spatial.domain,
  yrs=p$yrs,
  dyear=0.6 # 0.6*12 months = 7.2 = early July
)

# extract covariates and supplent survey data via lookups
set = aegis_db_lookup(
  X=set,
  lookupvars=covars,
  xy_vars=c("lon", "lat"),
  time_var="timestamp"
)

#  good data
ok = which(
  is.finite(set$totno) &
  is.finite(set$t) &
  is.finite(set$z) &
  is.finite(set$data_offset) &
  set$StrataID %in% sppoly$StrataID[sppoly$strata_to_keep]
)


# NOTE: R-INLA uses NA differently than other packages
# — NA in the response means no likelihood contribution, i.e. response is unobserved
# — NA in a fixed effect means no contribution to the linear predictor, i.e. the covariate is set equal to zero
# — NA in a random effect f(...) means no contribution to the linear predictor

APS = aegis_prediction_surface( aegis_data=res$means  )
APS$yr = as.numeric( APS$year)
APS$totno = NA
APS$data_offset = 1  # force to be density n/km^2
APS$tag = "predictions"

varstokeep = c( "totno", "StrataID", "yr", "t", "z", "data_offset", "tag" )

M = rbind( set[ok, varstokeep], APS[,varstokeep] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better

M$yr_factor = factor( as.character(M$yr) )
M$StrataID  = factor( M$StrataID, levels=levels(sppoly$StrataID ))
M$strata  = as.numeric( M$StrataID)
M$year  = as.numeric( M$yr_factor)

dtemps = seq( min(M$t, na.rm=TRUE), max(M$t, na.rm=TRUE), length.out=11 )
ddepths = c(2.5, 5, 10, 20, 40, 80, 160, 320, 640 )

M$ti = as.numeric( as.character( cut( M$t, breaks=dtemps, labels=diff(dtemps)/2 + dtemps[-length(dtemps)], include.lowest=TRUE ) ))

M$zi = as.numeric( as.character( cut( M$z, breaks=ddepths, labels=diff(ddepths)/2 + ddepths[-length(ddepths)], include.lowest=TRUE ) ))

M$iid_error = 1:nrow(M) # for inla indexing for set level variation


# generic PC priors
m = log( {set$totno / set$data_offset}[ok] )
m[!is.finite(m)] = min(m[is.finite(m)])

H = carstm_hyperparameters( sd(m), alpha=0.5, median(m) )
# H$prec$prec.intercept = 1e-9



# ------------------------------------------------
# Model 6: as in Model 5 but with main effects .... this is a full factorial model
# NOTE this is rank-deficient ... some factor combinations have insufficient data

# missing covariates/data forces more fiddling with 'kk'

fit = glm(
  formula = totno ~ offset(log(data_offset)) + 1 + StrataID:yr_factor + StrataID + yr_factor ,
  family=poisson(link="log"),
  data=set[ok, ]
)

s = summary(fit)
AIC(fit)  # 359785

#  prediction surface as a data frame
aps=APS
aps$yr_factor = factor( aps$year, levels=p$yrs )
kk = which(aps$StrataID %in% unique(set$StrataID[ok]))
preds = predict( fit, newdata=aps[kk,], type="response", na.action=na.omit, se.fit=TRUE )

aps$predictions = NA
aps$predictions[kk] = preds$fit
aps$predictions.sd = NA
aps$predictions.sd[kk] = preds$se.fit

# reformat predictions into matrix form
out = reformat_to_matrix(
  input = aps$predictions,
  matchfrom = list( StrataID=aps$StrataID, yr_factor=aps$yr_factor),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=levels(set$yr_factor))
)

# convert numbers/km to biomass/strata
RES$glm_poisson_totno_factorial = colSums( {out * weight_year * sppoly$sa_strata_km2}[sppoly$strata_to_keep , ], na.rm=TRUE )

lines( glm_poisson_totno_factorial ~ yr, data=RES, lty=5, lwd=4, col="red", type="b")

# map means adjusted by temperature and depth . .. as depth does not change, time dynamics in maps due to temperature and
iy = match( as.character(sppoly$StrataID), aps$StrataID )
vn = "pred"
sppoly@data[,vn] = aps$predictions[iy]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# ------------------------------------------------
# Model 6:
# "INLA full factorial"	fixed	fixed	-	Poisson	6136	unstable	-

# full factorial upon totno   -- base model for comparison with GLM  ..

# NOTE ::: iid_error is req to stabilize solution as
# there are additional sources of variation that the
# factorial model does not fully account for


fit = inla(
  formula = totno ~ 1 + StrataID:yr_factor + StrataID + yr_factor + f(iid_error, model="iid", hyper=H$iid),
  family = "poisson",
  data=M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  # control.inla=list(strategy="laplace", dz=0.25, diff.logdens=9, restart=3, npoints=11, cutoff=1e-5),
  control.inla=list(correct=TRUE, correct.verbose=FALSE ),  # adding this will make it 3.5hrs long
  control.fixed= H$fixed,
  verbose=FALSE
)
fit$dic$dic  # 34897
fit$dic$p.eff # 6848

s = summary(fit)


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$INLA.full.factorial = colSums( {out * weight_year * sppoly$sa_strata_km2}[sppoly$strata_to_keep, ], na.rm=TRUE ) / standardtow_sakm2
INLA.full.factorial

lines( INLA.full.factorial ~ yr, data=RES, lty=1, lwd=8, col="blue")





# ------------------------------------------------
# Model 9: as in Model 6 but full factorial with covariates
# note this is rank-deficient

# problem .. missing values in some strata/years (as with stratanal)
# .. but this is more likely with more covariates
# this means we cannot get a consistent estimates ( as with stratanal)
# NOTE: this is a full interaction model with strata and years ...

# missing covariates/data forces more fiddling


fit = glm(
  formula = totno ~ offset(log(data_offset)) + 1 + StrataID:yr_factor + StrataID + yr_factor +t + z,
  family=poisson(link="log"),
  data=set[ok, ]
)

s = summary(fit)
AIC(fit)  # 348154

aps = APS
aps$yr_factor = factor( aps$year, levels=p$yrs)
kk = which(aps$StrataID %in% unique(set$StrataID[ok]))
preds = predict( fit, newdata=aps[kk,], type="response", na.action=na.omit, se.fit=TRUE )

aps$predictions = NA
aps$predictions[kk] = preds$fit
aps$predictions.sd = NA
aps$predictions.sd[kk] = preds$se.fit

# reformat predictions into matrix form
out = reformat_to_matrix(
  input = aps$predictions,
  matchfrom = list( StrataID=aps$StrataID, yr_factor=aps$yr_factor),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=levels(set$yr_factor))
)

# convert numbers/km to biomass/strata
RES$glm_poisson_totno_env = colSums( {out * weight_year * sppoly$sa_strata_km2}[sppoly$strata_to_keep , ], na.rm=TRUE )

lines( glm_poisson_totno_env ~ yr, data=RES, lty=5, lwd=4, col="green", type="b")

# map means adjusted by temperature and depth . .. as depth does not change, time dynamics in maps due to temperature and
iy = match( as.character(sppoly$StrataID), aps$StrataID )
vn = "pred"
sppoly@data[,vn] = aps$predictions[iy]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


# Bottom line: an additive, linear, fixed effect model requires a lot of handwaving and fiddling
# esp when missing data and linear extrapolation that might not make sense at the extremes




# ------------------------------------------------
# Model 9b: as in Model 9 but GAM .. occasionally unstable .. optimizer is finding a flat area
# rank deficient

fit = mgcv::gam(
  formula = totno ~ offset(log(data_offset)) + 1 + StrataID:yr_factor + StrataID + yr_factor + s(t, bs="tp", k=3)  + s(z, bs="tp", k=3),
  family=poisson(link="log"),
  data=set[ok, ]
)

s = summary(fit)
AIC(fit)  # 340300

aps = APS
aps$yr_factor = factor( aps$year, levels=levels(set$yr_factor))
kk = which(aps$StrataID %in% unique(set$StrataID[ok]))
preds = predict( fit, newdata=aps[kk,], type="response", na.action=na.omit, se.fit=TRUE )

aps$predictions = NA
aps$predictions[kk] = preds$fit
aps$predictions.sd = NA
aps$predictions.sd[kk] = preds$se.fit

# reformat predictions into matrix form
out = reformat_to_matrix(
  input = aps$predictions,
  matchfrom = list( StrataID=aps$StrataID, yr_factor=aps$yr_factor),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=levels(set$yr_factor))
)

# convert numbers/km to biomass/strata
RES$gam_poisson_totno_env = colSums( {out * weight_year * sppoly$sa_strata_km2}[sppoly$strata_to_keep , ], na.rm=TRUE )

lines( gam_poisson_totno_env ~ yr, data=RES, lty=5, lwd=4, col="green", type="b")



# ------------------------------------------------
# Model 11:
# "INLA Envir 0"  -- base model using inla + envir --- mimicking the GAM model as much as possible


# simple factorial with totno and poisson
# improvement upon Model 6b .. INLA imputes missing data given proper data model .. less fiddling ..
# random effects.. `fewer params`


# NOTE: R-INLA uses NA differently than other packages
# — NA in the response means no likelihood contribution, i.e. response is unobserved
# — NA in a fixed effect means no contribution to the linear predictor, i.e. the covariate is set equal to zero
# — NA in a random effect f(...) means no contribution to the linear predictor


# interaction-only model does not really make sense here .. esp as there are other model components
fit = inla(
  formula =
    totno ~ 1 + offset( log( data_offset) )
      + StrataID:yr_factor + StrataID + yr_factor
      + f(iid_error, model="iid", hyper=H$iid)
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2),
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  control.fixed= H$fixed,
  verbose=TRUE
)

s = summary(fit)
s$dic$dic  #  34988
s$dic$p.eff # 6894

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$INLA.Envir.0 = colSums( {out * weight_year * sppoly$sa_strata_km2}[sppoly$strata_to_keep, ], na.rm=TRUE )

lines( INLA.Envir.0 ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )





# ------------------------------------------------
# Model 11a:
"INLA Envir 1"	rw2: temp+depth 	Poisson	5960	33652	4


# simple factorial with totno and poisson
# improvement upon Model 6b .. INLA imputes missing data given proper data model .. less fiddling ..
# random effects.. `fewer params`


# NOTE: R-INLA uses NA differently than other packages
# — NA in the response means no likelihood contribution, i.e. response is unobserved
# — NA in a fixed effect means no contribution to the linear predictor, i.e. the covariate is set equal to zero
# — NA in a random effect f(...) means no contribution to the linear predictor

if(0) {
  # alter priors for fixed effects
  H$fixed$prec = list( prior="pc.prec", param=c(1, 0.5) )  # NOTE: pc.priors are on sd scale ..
  H$fixed$prec.intercept = 1
}


# interaction-only model does not make sense here .. esp as there are other model components
fit = inla(
  formula =
    totno ~ 1 + offset( log( data_offset) )
      + f(strata, model="iid", group=year, hyper=H$iid)
      + f(iid_error, model="iid", hyper=H$iid)
      + f(year, model="iid", hyper=H$iid)
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$rw2)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$rw2),
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  control.fixed= H$fixed,
  verbose=TRUE
)

s = summary(fit)
s$dic$dic  #  33687
s$dic$p.eff # 5969

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$INLA.Envir.1 = colSums( {out * weight_year * sppoly$sa_strata_km2}[sppoly$strata_to_keep, ], na.rm=TRUE )

lines( INLA.Envir.1 ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")




dev.new(width=11, height=7)
col = c("slategray", "turquoise", "darkorange", "blue", "green", "darkred", "cyan", "darkgreen", "slateblue"  )
pch = c(20, 21, 22, 23, 24, 25, 26, 27, 20)
lty = c(1, 3, 4, 5, 6, 7, 1, 3, 4 )
lwd = c(4, 4, 4, 4, 4, 4, 4, 4, 4 )
type =c("l", "l", "l", "l", "l", "l", "l", "l", "l")

legend=c("Standard tow stratanal", "GLM factorial", "GLM Envir", "GAM Envir")

plot( stratanal_towdistance  ~ yr, data=RES, lty=lty[1], lwd=lwd[1], col=col[1], pch=pch[1], type=type[1], ylim=c(0,0.3e9), xlab="Year", ylab="kg")
lines( glm_poisson_totno_factorial ~ yr, data=RES, lty=lty[2], lwd=lwd[2], col=col[2], pch=pch[2], type=type[2])
lines( glm_poisson_totno_env ~ yr, data=RES, lty=lty[3], lwd=lwd[3], col=col[3], pch=pch[3], type=type[3])
lines( gam_poisson_totno_env ~ yr, data=RES, lty=lty[4], lwd=lwd[4], col=col[4], pch=pch[4], type=type[4])

ii = 1:4
legend("topright", legend=legend[ii], lty=lty[ii], col=col[ii], lwd=lwd[ii] )


dev.new(width=6, height=4)
hist( RES$INLA.Envir.1 / RES$stratanal_towdistance, breaks=20 )


cor( RES[, c("stratanal_towdistance", "glm_poisson_totno_factorial", "INLA.Envir.1")])

plot( RES[, c("stratanal_towdistance", "glm_poisson_totno_factorial", "INLA.Envir.1")])




# ---- bias in station selection:

set$strata_year = paste( set$StrataID, set$yr, sep=".")
zz = applyMean( set[, c("strata_year", "z")]  )
tt = applyMean( set[, c("strata_year", "t")]  )

APS$strata_year = paste( APS$StrataID, APS$yr, sep=".")
APS = merge( APS, zz, by="strata_year", all.x=TRUE, all.y=FALSE, suffixes=c("", ".set") )
APS = merge( APS, tt, by="strata_year", all.x=TRUE, all.y=FALSE, suffixes=c("", ".set") )

APS$z_diff = APS$z - APS$z.set
APS$t_diff = APS$t - APS$t.set

# sampling bias between survey and strata
dev.new(); plot( z ~ z.set, APS ); abline(0,1)
dev.new(); plot( t ~ t.set, APS ); abline(0,1)

dev.new(); plot( z_diff ~ yr, APS, pch=20, cex=0.75, col="slategray"); abline(h=0); out=data.frame(yr=p$yrs); out$zz=predict( loess(z_diff ~ yr, APS, span=0.05 ), newdata=out, se=FALSE); lines(zz~yr, out, lwd=4, col="red")
dev.new(); plot( t_diff ~ yr, APS, pch=20, cex=0.75, col="slategray"); abline(h=0); out=data.frame(yr=p$yrs); out$zz=predict( loess(t_diff ~ yr, APS, span=0.05 ), newdata=out, se=FALSE); lines(zz~yr, out, lwd=4, col="red")

APS$strata = as.numeric( APS$StrataID )
dev.new(); plot( z_diff ~ strata, APS, pch=20, cex=0.75, col="slategray"); abline(h=0); out=data.frame(strata=sort(unique(APS$strata))); out$zz=predict( loess(z_diff ~ strata, APS, span=0.05 ), newdata=out, se=FALSE); lines(zz~strata, out, lwd=4, col="red")
dev.new(); plot( t_diff ~ strata, APS, pch=20, cex=0.75, col="slategray"); abline(h=0); out=data.frame(strata=sort(unique(APS$strata))); out$zz=predict( loess(t_diff ~ strata, APS, span=0.05 ), newdata=out, se=FALSE); lines(zz~strata, out, lwd=4, col="red")


# ### end
