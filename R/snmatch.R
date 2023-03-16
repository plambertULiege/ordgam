## #######################################################################################
## GOAL: Fit a skew-Normal distribution to a density evaluated on a grid covering its support
## ---------------------------------------------------------------------------------------
## INPUT:
##   x: vector with the grid points covering the density support
##   y: density values on the grid <x> (possibly up to a multiplicative constant)
##   log: indicates if the log of the density is provided in <y> ; FALSE by default
##
## OUTPUT: list with
##   dp: skew-Normal parameters for the fitted density
##   fitted.moments: mean, variance, skewness, kurtosis of the fitted skew-t
##   target.moments: mean, variance, skewness, kurtosis of the target density approximated by quadrature
## #######################################################################################
## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2022)
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' @keywords internal
snmatch = function(x,y,log=FALSE){
    if (log) y = exp(y-max(y))
    cst = 1 / tail(cubicBsplines::trapeze(x,y),1)
    xg = x ; yg = y * cst
    ##
    ## Target Mean, Variance, Skewness, Kurtosis
    m1 = tail(cubicBsplines::trapeze(xg, xg*yg),1)
    m2 = tail(cubicBsplines::trapeze(xg, (xg-m1)^2*yg),1)
    g1 = tail(cubicBsplines::trapeze(xg, (xg-m1)^3*yg),1) / m2^1.5
    target.moments = c(m1,m2,g1)
    names(target.moments) = c("mu","sigma2","gam1")
    ##
    kappa2 = pi * (abs(g1) / (sqrt(2)*(4-pi)))^1.5
    psi2 = kappa2 / (1 + 2*kappa2/pi)
    psi = sign(g1) * sqrt(psi2)
    alpha = psi / sqrt(1-psi2)
    omega = sqrt(m2 / (1 - 2*psi2/pi))
    xi = m1 - omega * sqrt(2/pi) * psi
    ##
    dp = c(xi=xi, omega=omega, alpha=alpha)
    ##
    dp2moments = function(dp){
        dp = unname(dp)
        xi = dp[1] ; omega = dp[2] ; alpha = dp[3]
        bn = sqrt(2/pi)
        delta = alpha/sqrt(1+alpha^2)
        ##
        mu = xi + omega * bn * delta
        sz = sqrt(1 - (bn*delta)^2)
        sigma2 = omega^2 * sz^2
        gam1 = .5 * (4-pi) * omega^3 * bn^3 * delta^3 / sigma2^1.5
        ##
        ans = c(mu,sigma2,gam1)
        names(ans) = c("mu","sigma2","gam1")
        return(ans)
    }
    ## Prepare output
    ans = list(dp=dp,
               fitted.moments=dp2moments(dp),
               target.moments=target.moments)
    return(ans)
}
