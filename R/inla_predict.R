
inla_predict = function( fit, newdata ) {

    theta = fit$mode$theta
    formula = deparse(fit$.args$formula)
    
    newfit = inla(
        formula = formula,
        data = newdata,
        control.mode = list(theta = fit$mode$theta, restart = FALSE)
    )

    return(newfit)
}
