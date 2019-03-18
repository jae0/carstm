strata_timeseries = function( set, ... ) {

  params = list(...)

  res = data.frame()

  for (yr in params[["yr"]]){
    i = which(set$yr == yr)
    yhi = split(set[ i, params[["variable"]] ], set[i,"strat"] ) #split the variable by strata
    nh = as.vector(sapply(yhi, length)) #numer of tows per strata
    nhws = sapply(yhi, function(x) length(x [x > 0])) #calculate the number of samples > 1 in each strata
    Nh = sapply( split(set[ i, "nh" ], set[i,"strat"] ), mean, na.rm=TRUE) #split the variable by strata to get trawalable units in each strata
    Wh = Nh / sum(Nh) #strata percent of the total area in trawlable units

    #Calculate Stratified Estimates and Confidence intervals
    #-------------------------------------------------------------------------------------
    yh = as.vector(sapply(yhi, mean)) #mean of variable for each strata
    yst = sum(Wh * yh, na.rm = TRUE) #sum of the mean of the variable for each strata, multiplied by percent area of each strata

    sh = as.vector(sapply(yhi, var)) #calculate variance of each variable, per strata
    se.yst = sqrt(sum((((Nh * (Nh - nh))/sum(Nh)^2) * sh)/nh, na.rm = TRUE)) #calculate standard error
    ah = (Nh * (Nh - nh))/nh #
    df.yst = (sum(ah * sh, na.rm = TRUE)^2)/(sum(((ah * sh)^2)/(nh - 1), na.rm = TRUE)) #degrees of freedom

    #Calculate the confidence interval, based on the t distribution, rather than the normal distribution
    #qt fuction calculates the value for the 95% confidence interval by looking up the t distribution,
    #based on the degrees of freedom calculated above. This is multiplied by the standard error calculated above
    #formulas are: Lower limit = M - (tCL)(sM), Upper limit = M + (tCL)(sM)
    #------------------------------------------------------------------------------------
    ci.yst = yst + (c(qt(params[["alpha.t"]]/2, df.yst), -qt(params[["alpha.t"]]/2, df.yst)) * se.yst) #confidence interval

    #Calculate Design Weighted Area Occupied
    dwao = sum(Wh*(nhws / nh)) * sum(Nh) * 0.011801

    #Calculate Gini Index  .. Needed? JC?
    gi = NA
    gi = gini(x=yh, y=Nh)

    #Calculate the total population
    Yst = yst * sum(Nh)
#    Yst = sum( yst * Nh)

    #Use mirror match method (BWR) to calculate confidence intervals in Bias Corrected (BC) method to calculate confidence interval
    #-------------------------------------------------------------------------
    # call = match.call(expand = FALSE)

    nrs = 1:params[["nresamp"]]
    out = matrix(0, nrow=params[["nresamp"]]+1, ncol=3, dimnames = list(c("Actual", nrs), c("Mean", "Variance",'gini')))
    out[1, ] = c(yst, (se.yst)^2, gi)

    kh = (nh - 1)/(1 - nh/Nh)
    ph = ((1/kh) - (1/ceiling(kh)))/((1/floor(kh)) - (1/ceiling(kh)))


    for (i in nrs) {
      yhib = bwr.boot(yhi, kh, ph, sample, replace = TRUE, simplify = FALSE)
      yhib[nh == 1] = yhi[nh == 1]
      nhws = sapply(yhib, FUN = function(x) sum(x > 0))
      n_yhib = as.vector(sapply(yhib, length))
      var_yhib = as.vector(sapply(yhib, var))
      mean_yhib = as.vector(sapply(yhib, mean))
      out[i + 1, ] = c(
        sum(Wh * mean_yhib, na.rm = TRUE),  # mean`
        sum( (((Nh * (Nh - n_yhib))/sum(Nh)^2) * var_yhib ) / n_yhib, na.rm = TRUE),  # variance
        gini(x = mean_yhib, y = Nh)  # gini
      )
    }

    orig.mean  = out[1, 1]
    orig.var   = out[1, 2]
    boot.means = out[nrs+1, 1]
    boot.vars  = out[nrs+1, 2]
    gi =         out[nrs+1, 3]
    # call = call
    # method = method

    #Summary of Bootstrapped means
    #-------------------------------------------------------------------------------------------------------
    options(digits = 4)
    boot.est = mean(boot.means) #bootstrap mean
    gini.mean = mean(gi)
    ci.boot=list()

    loc.bc = sum(boot.means < boot.est)
    lim.bc = sort(boot.means)[c(loc.bc, loc.bc + 1)]
    z0 = (loc.bc + ((boot.est - lim.bc[1])/(lim.bc[2] - lim.bc[1])))
    z0 = qnorm(z0/length(boot.means))
    probs.z0 = pnorm(qnorm(c(params[["alpha.b"]]/2, (1 - params[["alpha.b"]]/2), 0.5)) + 2 * z0)

    ci.boot.meanBCA = BCa(boot.means,0.01,alpha=c(0.025,0.975),mean)
    ci.boot[[1]] = ci.boot.mean = c(ci.boot.meanBCA[4:5],ci.boot.meanBCA[3])  #ci.boot.mean = quantile(boot.means, probs = probs.z0)
    ci.boot.gini = quantile(gi, probs = c(params[["alpha.b"]]/2, (1 - params[["alpha.b"]]/2), 0.5), na.rm=T)


    #Print out the yearly estimates and write them to a data frame
    #--------------------------------------------------------------------------------------------------------
    options(digits = max(options()$digits - 5, 5))

    if (params[["prints"]]) {cat(
      "\n",
      "Pop Total =", format(Yst), "\n",
      "Original Mean =", format(orig.mean), "\n",
      "Year =", yr, "\n",
      "Original Variance =", format(orig.var), "\n",
      "Number of bootstraps = ", length(boot.means), "\n",
      "Bootstrap Mean=", format(boot.est), "\n",
      "Variance of Bootstrap Mean=", format(var(boot.means)), "\n",
      "CI Method=", c(params[["CI.method"]]), "\n",
      "CI's for alpha=", params[["alpha.b"]], "are ", format(ci.boot.mean[1:2]), "\n",
      "Length =", format(ci.boot.mean[2] - ci.boot.mean[1]), "\n",
      "Shape=", format(log((ci.boot.mean[2] - ci.boot.mean[3])/(ci.boot.mean[3] - ci.boot.mean[1]))), "\n",
      "Resample Method = ", params[["method"]], "\n")
    }

    ci.boot[[2]] = ci.boot.gini

    n.row = data.frame(
      speciesname=params[["speciesname"]],
      year = as.numeric(yr),
      pop.total = Yst,
      variable = params[["variable"]],
      orig.mean = as.numeric(format(orig.mean)),
      boot.mean = as.numeric(format(boot.est)),
      var.boot.mean = as.numeric(format(var(boot.means))),
      lower.ci = as.numeric(format(ci.boot.mean[1])),
      upper.ci = as.numeric(format(ci.boot.mean[2])),
      length = as.numeric(format(ci.boot.mean[2] - ci.boot.mean[1])),
      dwao = dwao,
      gini = gini.mean,
      lower.ci.gini = format(ci.boot.gini[1]),
      upper.ci.gini = format(ci.boot.gini[2])
    )

    res = rbind(n.row, res )

  }

  library(zoo)
  #Calculate 3 yr Mean
  res$mean.3.yr = zoo::rollapply(res$boot.mean, 3, mean, fill=NA, align="right")

  #Calculate Mean Lines
  res$median = median(res$boot.mean)
  res$median.50 = median(res$boot.mean) * 0.5
  #res$gm.40 = geometric.mean (res$boot.mean) * 0.4

  print(res)

  return(res)
}
