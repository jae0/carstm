
# ------------------------------------------------
# Atlantic cod comparison of naive strata-based averages
# This replicates standard groundfish strata-based estimation of means and totals
# "standard" random-stratified estimation functions (based on stratanal and bootstrap estimation techniques )
# derived from numerous authors .. this code has been optimized a bit more for clarity and speed by JSC

# ------------------------------------------------
# load data environment

RES = data.frame(yr=1970:2017)  # collect model comparisons in this data frame
if (0) {
  fn = file.path( getwd(), "RES.rdata" )
  # save(RES, file=fn)
  # load(fn)
}

yrs =2000:2018

for (tu in c( "standardtow", "towdistance", "sweptarea") ) {


    # construct basic parameter list defining the main characteristics of the study and some plotting params
    p = carstm::carstm_parameters(
      label ="Atlantic cod summer standardtow",
      speciesname = "Atlantic_cod",
      groundfish_species_code = 10,   #  10= cod
      yrs = yrs,
      polygon_source = "pre2014",   # "pre2014" for older
      areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia"),  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
      trawlable_units = tu
    )



    # specific selection params required for survey.db(DS="filter") data selection mechanism
    p = aegis.survey::survey_parameters(
      p=p,
      selection=list(
        biologicals=list(
          spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code )
        ),
        survey=list(
          data.source="groundfish",
          yr = p$yrs,      # time frame for comparison specified above
          months=6:8,
          # dyear = c(150,250)/365, #  summer = which( (x>150) & (x<250) ) , spring = which(  x<149 ), winter = which(  x>251 )
          settype = 1,
          gear = c("Western IIA trawl", "Yankee #36 otter trawl"),
          strata_toremove=c("Gulf", "Georges_Bank", "Spring", "Deep_Water"),  # <<<<< strata to remove from standard strata-based analysis
          polygon_enforce=TRUE,
          ranged_data="dyear"
        )
    )



    if (0) {
      # used for DEBUGGING:  access via strata_dataformat .. ie. directly from groundfish.db .. gscat
      # gscat does not have the data corrections due to miscoding etc, vessel-species "catchability" corrections, etc
      # .. but the totals are nearly identical to survey.db access
      p$selection$biologicals$spec = groundfish_species_code
      p$selection$biologicals$spec_bio = NULL

      set = strata_dataformat( p=p )   # return values in kg or no per set
      # dim(set) # [1] 1682   45
      # sum(set$totwgt) # [1] 9683
      # sum(set$totno)  # [1] 15261
    }


    # ------------------------------------------------
    # NOTE polygon areas of strata are predetermined in GSSTRATUM and used by "stratanal"
    # This behaviour is mimicked in survey.db .. although it does not have to be
    # .. The CAR variation computes surfaces directly from polygons
    # Using this approach is better in that there is more filter control
    # Results for the basic test cases are essentially identical to "stratanal" (via "strata_dataformat", above)
    # but faster and more QA/QC done on the input data
    set = aegis.survey::survey.db( p=p, DS="filter", add_groundfish_strata=TRUE )   # return values in kg or no per set

    if (0) {
      # the above merges based upon AUID's designated in groundfish tables.  Alternatively one can use positions directly:
      set = survey.db( p=p, DS="filter" )
      # categorize Strata
      sppoly = areal_units( areal_units_source="stratanal_polygons",  areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km, timeperiod="pre2014" )
      sppoly$strata_to_keep = ifelse( as.character(sppoly$AUID) %in% strata_definitions( c("Gulf", "Georges_Bank", "Spring", "Deep_Water") ), FALSE,  TRUE )

      o = over( SpatialPoints( set[,c("lon", "lat")], sp::CRS(projection_proj4string("lonlat_wgs84")) ), spTransform(sppoly, sp::CRS(projection_proj4string("lonlat_wgs84")) ) ) # match each datum to an area
      set$AUID = o$AUID
      o = NULL
      set = set[ which(!is.na(set$AUID)),]
    }


    #  dim(set) # [1] 1684   39
    #  sum(set$totwgt) # [1] 9691
    #  sum(set$totno) # [1] 15263
    #  sum(set$totwgt_adjusted) # [1] 242838
    #  sum(set$totno_adjusted) # [1] 380823


    # compute no of trawlable units for each stratum

    # "strat" and "nh" are used by stratanal
    set$strat = set$AUID

    ft2m = 0.3048
    m2km = 1/1000
    nmi2mi = 1.1507794
    mi2ft = 5280
    standardtow_sakm2 = (41 * ft2m * m2km ) * ( 1.75 * nmi2mi * mi2ft * ft2m * m2km )  # surface area sampled by a standard tow in km^2  1.75 nm
    # = 0.0405 and NOT 0.011801 .. where did this come from?
    # set up trawlable units used in stratanal
    set$nh = switch( p$trawlable_units,
      sweptarea = as.numeric(set$au_sa_km2) * set$cf_cat, # convert strata area to trawlable units 41ft by 1.75 nm, divide area by sweptarea
      standardtow = as.numeric(set$au_sa_km2) / standardtow_sakm2, # convert strata area to trawlable units 41ft by 1.75 nm, divide area by 0.011801
      towdistance = as.numeric(set$au_sa_km2) / set$sa_towdistance # convert strata area to trawlable units 41ft by 1.75 nm, divide area by 0.011801
    )
    i = which(!is.finite(set$nh))
    set$nh[i] = as.numeric(set$au_sa_km2[i]) / standardtow_sakm2  # override missing with "standard set" .. also if no trawlable_units set



    # ------------------------------------------------
    # Random stratified estimates from a faster variation of Michelle's code
    results_basic = strata_timeseries(
      set = set,
      speciesname = p[["label"]],
      yr = p$yrs,
      variable ="totwgt",
      alpha.t = 0.05, # confidence interval eg. 0.05 = 95%, 0.1 = 90%
      alpha.b = 0.05,
      nresamp = 1000,
      prints=TRUE
    )

    (results_basic)

    RES[,paste("stratanal", tu, sep="_")] = results_basic$pop.total[match(results_basic$year, RES$yr)]
    # plot( stratanal ~ yr, data=RES, lty=5, lwd=4, col="red", type="b", ylim=c(0,8e8))
    # lines ( stratanal ~ yr, data=RES, lty=5, lwd=4, col="red", type="b", ylim=c(0,8e8))

}

dev.new(width=11, height=7)
col = c("slategray", "turquoise", "darkorange" )
pch = c(20, 21, 22)
lty = c(1, 3, 4 )
lwd = c(4, 4, 4)
type =c("l", "l", "l")
plot( stratanal_standardtow  ~ yr, data=RES, lty=lty[1], lwd=lwd[1], col=col[1], pch=pch[1], type=type[1], ylim=c(0,2.6e8), xlab="Year", ylab="kg")
lines( stratanal_towdistance ~ yr, data=RES, lty=lty[2], lwd=lwd[2], col=col[2], pch=pch[2], type=type[2])
lines( stratanal_sweptarea ~ yr, data=RES, lty=lty[3], lwd=lwd[3], col=col[3], pch=pch[3], type=type[3])
legend("topright", legend=c("Standard tow", "Length adjusted", "Length & width adjusted"), lty=lty, col=col, lwd=lwd )


dev.new(width=6, height=4)
hist( RES$stratanal_towdistance / RES$stratanal_standardtow, breaks=20 )

dev.new(width=6, height=4)
hist( RES$stratanal_sweptarea / RES$stratanal_standardtow, breaks=20 )


cor( RES[, c("stratanal_standardtow", "stratanal_towdistance", "stratanal_sweptarea")])

plot( RES[, c("stratanal_standardtow", "stratanal_towdistance", "stratanal_sweptarea")])



# ------------------------------------------------
# these are Michelle's results: (base access of gcat without correction factors for boat, species, etc)
        speciesname year pop.total variable orig.mean boot.mean var.boot.mean lower.ci upper.ci   length
2.5%   COD ATLANTIC 2017  14593959 totwgt_sd    3.4420    3.4258       3.61840   3.3099   3.5451 0.235210 21313 0.81003
2.5%34 COD ATLANTIC 2016  27531380 totwgt_sd    6.4932    6.3838      23.36500   6.0890   6.6869 0.597900 20779 0.89443
2.5%33 COD ATLANTIC 2015   8915342 totwgt_sd    2.1027    2.0970       0.17429   2.0716   2.1232 0.051683 24031 0.71914
2.5%32 COD ATLANTIC 2014  28570078 totwgt_sd    6.7382    6.8005      13.48700   6.5766   7.0328 0.456180 20416 0.88363
2.5%31 COD ATLANTIC 2013  12550459 totwgt_sd    2.9600    2.9837       1.13150   2.9189   3.0504 0.131470 24549 0.76574
2.5%30 COD ATLANTIC 2012   9538831 totwgt_sd    2.2497    2.2245       0.37251   2.1873   2.2630 0.075729 22789 0.76290
2.5%29 COD ATLANTIC 2011  35724538 totwgt_sd    8.4256    8.4033      20.51000   8.1265   8.6906 0.564150 24609 0.79815
2.5%28 COD ATLANTIC 2010  44532221 totwgt_sd   10.5030   10.3040      43.71900   9.9038  10.7220 0.817780 28273 0.83744

# ------------------------------------------------
# these are with "standard tow" assumptions:
                 speciesname year pop.total variable orig.mean boot.mean var.boot.mean lower.ci upper.ci   length  dwao    gini lower.ci.gini upper.ci.gini mean.3.yr median median.50
2.5%  Cod summer standardtow 2017  14863259   totwgt    3.5136    3.4339       3.80320   3.3145   3.5572 0.242740 21451 0.77906
2.5%7 Cod summer standardtow 2016  21430734   totwgt    5.0662    5.1501      11.10100   4.9457   5.3587 0.412980 20681 0.84814
2.5%6 Cod summer standardtow 2015   8723439   totwgt    2.0622    2.0598       0.16459   2.0347   2.0851 0.050439 23705 0.65717
2.5%5 Cod summer standardtow 2014  27156331   totwgt    6.4197    6.4166      13.37500   6.1928   6.6459 0.453100 20786 0.83574
2.5%4 Cod summer standardtow 2013  12288438   totwgt    2.9050    2.9377       1.12560   2.8729   3.0045 0.131620 24411 0.73470
2.5%3 Cod summer standardtow 2012   9105517   totwgt    2.1525    2.1576       0.38184   2.1198   2.1957 0.075942 22465 0.73825
2.5%2 Cod summer standardtow 2011  34542306   totwgt    8.1657    8.1041      20.12700   7.8330   8.3863 0.553380 24584 0.77513
2.5%1 Cod summer standardtow 2010  42903020   totwgt   10.1420   10.4090      42.51400  10.0130  10.8180 0.805170 28360 0.82533

# ------------------------------------------------
# towed distance
                 speciesname year pop.total variable orig.mean boot.mean var.boot.mean lower.ci upper.ci   length  dwao    gini
2.5%  Cod summer towdistance 2017  15025543   totwgt    3.4420    3.4205       3.76150   3.3030   3.5434 0.240360 22219 0.77524
2.5%7 Cod summer towdistance 2016  23732430   totwgt    5.5022    5.4666      14.00500   5.2375   5.7001 0.462570 21266 0.84436
2.5%6 Cod summer towdistance 2015   8798976   totwgt    2.0690    2.0732       0.15893   2.0490   2.0982 0.049255 23901 0.65307
2.5%5 Cod summer towdistance 2014  28503738   totwgt    6.5806    6.6708      12.70500   6.4516   6.8923 0.440700 21401 0.83786
2.5%4 Cod summer towdistance 2013  12434510   totwgt    2.8474    2.8555       1.02700   2.7944   2.9199 0.125560 25127 0.73675
2.5%3 Cod summer towdistance 2012   9340895   totwgt    2.1474    2.1441       0.37216   2.1066   2.1824 0.075768 23436 0.73386
2.5%2 Cod summer towdistance 2011  35721843   totwgt    8.1870    8.0592      21.59400   7.7759   8.3539 0.578040 25474 0.77066
2.5%1 Cod summer towdistance 2010  43790809   totwgt    9.9849    9.8188      43.78700   9.4161  10.2420 0.825920 29684 0.81764

# ------------------------------------------------
# sweptarea
                speciesname year pop.total variable orig.mean boot.mean var.boot.mean lower.ci upper.ci   length  dwao
2.5%  Cod summer sweptarea 2017  14584703   totwgt    3.5108    3.5984       4.12750   3.4741   3.7256 0.251540 20961
2.5%7 Cod summer sweptarea 2016  20677264   totwgt    5.0506    5.0506      11.32400   4.8465   5.2649 0.418400 19806
2.5%6 Cod summer sweptarea 2015   7397592   totwgt    1.9220    1.9240       0.13729   1.9011   1.9471 0.045981 20763
2.5%5 Cod summer sweptarea 2014  25264103   totwgt    6.5155    6.5552      14.94200   6.3205   6.7999 0.479430 18288
2.5%4 Cod summer sweptarea 2013  10290871   totwgt    2.7651    2.7570       1.12270   2.6915   2.8227 0.131200 20251
2.5%3 Cod summer sweptarea 2012   8839376   totwgt    2.1085    2.1188       0.36425   2.0816   2.1566 0.074987 22206
2.5%2 Cod summer sweptarea 2011  34866336   totwgt    8.4871    8.0913      24.96300   7.7888   8.4061 0.617280 23648
2.5%1 Cod summer sweptarea 2010  45648142   totwgt   10.7420   10.6010      52.95200  10.1560  11.0630 0.906710 28503


# TODO basic corelations and plots, summarizing the above


set$strata_year = paste( set$AUID, set$yr, sep=".")
nn = applySummary( set[, c("strata_year", "totno")]  )

V = expand.grid( AUID=levels(set$AUID), yr=sort( unique(set$yr) ) )
V$strata_year = paste( V$AUID, V$yr, sep=".")
V = merge( V, nn, by="strata_year", all.x=TRUE, all.y=FALSE, suffixes=c("", ".totno") )

dev.new(); plot( log(totno.mean) ~ log(totno.sd), V ); abline(0,1)


### end
