
# testing diseasemapping

install.packages("diseasemapping", repos="http://R-Forge.R-project.org")
install.packages("mapmisc", repos="http://R-Forge.R-project.org")

require(diseasemapping)
require(mapmisc)

set.seed(0)

# data files also saved to local directory ~/projects/diseasemapping
load( url("http://pbrown.ca/data/nuts.RData") )  # xFinal
euroNuts = readRDS( url( "http://pbrown.ca/data/euroNuts.rds" ) )


# copied from 
euroNutsSub = euroNuts[!euroNuts$CNTR_CODE %in% c("CY",'TR','IE'), c('NUTS_ID','NUTS_NAME','CNTR_CODE')]
# euroNutsSub@data[1:2,]


xSub = xFinal[xFinal$geo %in% euroNutsSub$NUTS_ID, ]
xSub$logPop = log(xSub$pop)
theRates = glm(deaths ~ offset(logPop) + 0+age:sex, family='poisson', data=xSub)
xSub$expected =  predict(theRates, xSub[,c('logPop','age','sex')], type='response')

euroMelt = reshape2::melt(xSub, 
  id.vars = c('geo','gdp','pctEmp','manure','year', 'sex'),
  measure.vars = c('expected','deaths')
)

MO = reshape2::dcast(
  euroMelt,
  geo + year + gdp + pctEmp + manure ~ variable + sex,
  fun.aggregate = sum
)

MO$smr = MO$deaths_F / MO$expected_F 
MO$logExpected_F = log(MO$expected_F)
MO$gdpTrans = log(MO$gdp / median(MO$gdp, na.rm=TRUE))/log(2)
MO$yearFac = factor(MO$year)
MO$yearInt = as.integer(MO$yearFac)

MO[1:3, c('year','geo','deaths_F', 'expected_F','logExpected_F','gdp')]



#Neighbourhood matrix 
adjMat = spdep::poly2nb(euroNutsSub, row.names=euroNutsSub$NUTS_ID)


if (0) {

		ls()

		# 2 objects nuts and xFinal
		str( xFinal)  # data.frame

		class( nuts) # list
		length( nuts) # length 3
		str( nuts[[1]]@data)
		str( nuts[[2]]@data)
		str( nuts[[3]]@data)

		class( euroNuts )
		str(euroNuts@data)

		euroBg = mapmisc::openmap(euroNutsSub, path='stamen-watercolor')


		smrCol = mapmisc::colourScale(
			MO$smr, breaks=c(0, 0.4, 0.6, 0.8, 1, 1.5, 2, 7.5), 
			style = 'fixed',
			col='Spectral', rev=TRUE)

		obsCol = mapmisc::colourScale(
			MO$deaths_F, breaks=10, style='equal', dec=-1, digits=1,
			transform='sqrt',
			col='Spectral', rev=TRUE)

		obsCol$breaks[length(obsCol$breaks)] = ceiling(max(c(
			obsCol$breaks[length(obsCol$breaks)],
			max(MO$expected_F, na.rm=TRUE)))/10)*10

		expCol = mapmisc::colourScale(
			MO$expected_F, breaks=obsCol$breaks, style='fixed',
			col=obsCol$col)

		euroSetup = function() {
			mapmisc::map.new(euroNutsSub, buffer = c(0.01,0.01,-0.1,-1)*1000*1000)
			raster::plotRGB(euroBg, add=TRUE)
		}

		Syears= sort(unique(MO$year))

		for(D in Syears) {
			euroSetup()
			thisYear = which(MO$year == D)
			geoHere = match(MO[thisYear, 'geo'], euroNutsSub$NUTS_ID)
			sp::plot(euroNutsSub, col=obsCol$plot[thisYear][geoHere], border=T, add=TRUE)

			euroSetup()
			thisYear = which(MO$year == D)
			geoHere = match(MO[thisYear, 'geo'], euroNutsSub$NUTS_ID)


			knitr::kable(fitStOD$parameters$summary, digits=2)

			qCols = paste0(c(0.5, 0.025, 0.975), 'quant')
			myTable = rbind(
				exp(fitStOD$parameters$summary[-(7:8),qCols]),
				fitStOD$parameters$summary[7:8,qCols],
				sdW = 1/sqrt(fitStOD$inla$summary.hyperpar[
					'Precision for geo', qCols[c(1,3,2)]]))
			knitr::kable(myTable, digits=2)




			dim(fitStOD$data)
			length(euroNutsSub)
			names(fitStOD$data)[1:9]



			randomCol = mapmisc::colourScale(exp(fitStOD$data$random.0.5quant),
				breaks = 9, digits=1.5, col='Spectral', rev=TRUE, 
				style='equal', transform='log')
			mapmisc::map.new(euroNutsSub, buffer = c(0,0,0,-1000*1000))
			sp::plot(euroNutsSub, col=randomCol$plot, border=NA, add=TRUE)
			mapmisc::legendBreaks('topright', randomCol, cex=0.7, bty='n')


			levels(fitStOD$inla$.args$data$yearFac)[1]


			fittedCol = mapmisc::colourScale(
				exp(fitStOD$data$fitted.0.5quant),
				breaks = 9, digits=1.5, 
				col='Spectral', rev=TRUE, 
				style='equal', transform='log')
			mapmisc::map.new(euroNutsSub,
				buffer = c(0,0,0,-1000*1000))
			sp::plot(euroNutsSub, 
				col=fittedCol$plot, 
				border=NA, add=TRUE)
			mapmisc::legendBreaks('topleft', 
				fittedCol, cex=0.8, bty='n')
			sp::plot(euroNutsSub, col=smrCol$plot[thisYear][geoHere], border=T, add=TRUE)
		}

		tSeq = seq(0,3, len=21)
		z = cbind(0.1*rnorm(length(tSeq)), rnorm(length(tSeq)))

		matplot(tSeq, z, lty=1, lwd=2, xlab='time', ylab='rw0', type='l')

		z1 = apply(rbind(c(1,-1), z), 2, cumsum)[-1,]
		matplot(tSeq, z1, lty=1, lwd=2, xlab='time', ylab='rw1', type='l')

		z2 = apply(rbind(c(-6,2), z1), 2, cumsum)[-1,]
		matplot(tSeq[1:nrow(z2)], z2, lty=1, lwd=2, xlab='time', ylab='rw2', type='l' )


}



fitStOD = diseasemapping::bym(
  deaths_F ~ offset(logExpected_F) + yearFac + gdpTrans +
    f(geo, model='iid', replicate=yearInt, 
      prior='pc.prec', param=c(log(1.5), 0.5)),
  data = MO, region.id = 'geo', adjMat = adjMat, family='poisson',
  prior = list(sd = c(log(1.5), 0.5), propSpatial =c(0.5, 0.5)))


fitStFull = INLA::inla(
  deaths_F ~ offset(logExpected_F) + yearFac + gdp + 
    f(space, model='bym2', graph=NB_graph, 
      hyper = list(
        theta1 = list(prior='pc.prec', param=c(log(1.5), 0.5)),
        theta2 = list(prior = 'pc', param = c(0.5, 0.5))) ) +
    f(space_time, model='bym2', graph=NB_graph, 
      hyper = list(
        theta1 = list(prior='pc.prec', param=c(log(1.5), 0.5)),
        theta2 = list(prior = 'pc', param = c(0.1, 0.5))),
      group=yearInt, control.group = list(model='ar1', hyper = 
        list(theta=list(prior='pccor0', param=c(0.1, 0.5)))) ),
  data = MO, family='poisson',
#  control.inla = list(strategy='gaussian', int.strategy='eb'),
  control.compute = list(config=TRUE), 
  control.mode = list(theta = c(3.52, 0.38, 4.9, 0.17, 3.5),  restart=TRUE), 
  verbose=TRUE )











