## #######################################################################################
## GOAL: Fit a skew-t distribution to a density evaluated on a grid covering its support
## ---------------------------------------------------------------------------------------
## INPUT:
##   x: vector with the grid points covering the density support
##   y: density values on the grid <x> (possibly up to a multiplicative constant)
##   log: indicates if the log of the density is provided in <y> ; FALSE by default
##   maxiter: maximum number of iterations for the fitting algorithm
##   tol: tolerance in the moment matching procedure
##
## OUTPUT: list with
##   dp: skew-t parameters for the fitted density
##   fitted.moments: mean, variance, skewness, kurtosis of the fitted skew-t
##   target.moments: mean, variance, skewness, kurtosis of the target density approximated
##                     by quadrature
##   niter: required number of iterations
##   maxiter: see INPUT
##   tol: see INPUT
##   converged: convergence indicator
## #######################################################################################
## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2022)
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' @keywords internal
stmatch = function(x,y,log=FALSE,maxiter=100,tol=1e-4){
    if (log) y = exp(y-max(y))
    cst = 1 / tail(cubicBsplines::trapeze(x,y),1)
    xg = x ; yg = y * cst
    ##
    ## Target Mean, Variance, Skewness, Kurtosis
    m1 = tail(cubicBsplines::trapeze(xg, xg*yg),1)
    m2 = tail(cubicBsplines::trapeze(xg, (xg-m1)^2*yg),1)
    g1 = tail(cubicBsplines::trapeze(xg, (xg-m1)^3*yg),1) / m2^1.5
    g2 = tail(cubicBsplines::trapeze(xg, (xg-m1)^4*yg),1) / m2^2 - 3
    target.moments = c(m1,m2,g1,g2)
    names(target.moments) = c("mu","sigma2","gam1","gam2")
    ##
    bn.fun = function(nu) sqrt(nu/pi) * exp(lgamma(.5*(nu-1))-lgamma(.5*nu))
    sz.fun = function(alpha,nu) sqrt(nu/(nu-2) - (bn.fun(nu)*alpha/sqrt(1+alpha^2))^2)
    ## Starting values for xi, omega, alpha, nu
    xi = m1 ## Location
    alpha = 0 ## No skewness
    nu = 10    ## d.f.
    sz = sz.fun(alpha,nu)
    omega = sqrt(m2) / sz ## Dispersion
    dp = c(xi=xi, omega=omega, alpha=alpha, nu=nu)
    ##
    dp2moments = function(dp){
        dp = unname(dp)
        xi = dp[1] ; omega = dp[2] ; alpha = dp[3] ; nu = dp[4]
        bn = sqrt(nu)*exp(lgamma(.5*(nu-1))-lgamma(.5*nu)) / sqrt(pi)
        delta = alpha/sqrt(1+alpha^2)
        ##
        mu = xi + omega * bn * delta
        sz = sqrt(nu/(nu-2) - (bn*delta)^2)
        sigma2 = omega^2 * sz^2
        gam1 = bn*delta / sz^3 * (nu*(3-delta^2)/(nu-3) - 3*nu/(nu-2) + 2*(bn*delta)^2)
        gam2 = (1/sz^4) * (3*nu^2/(nu-2)/(nu-4) - 4*(bn*delta)^2*nu*(3-delta^2)/(nu-3) +6*(bn*delta)^2*nu/(nu-2) -3*(bn*delta)^4) - 3
        ##
        ans = c(mu,sigma2,gam1,gam2)
        names(ans) = c("mu","sigma2","gam1","gam2")
        return(ans)
    }
    ##
    update.dp = function(dp){
        xi = dp[1] ; omega = dp[2] ; alpha = dp[3] ; nu = dp[4]
        ## Update xi
        xi <- c(m1 - omega * bn.fun(nu) * alpha/sqrt(1+alpha^2))
        ##
        ## Update omega
        omega <- c(sqrt(m2) / sz.fun(alpha,nu))
        ##
        ## ## Update alpha & nu
        ## foo = function(coef){
        ##     alpha = coef[1] ; nu = abs(coef[2])
        ##     delta = alpha/sqrt(1+alpha^2)
        ##     bn = bn.fun(nu)
        ##     sz = sz.fun(alpha,nu)
        ##     ## alpha
        ##     ans1 = (g1 - bn*delta / sz^3 * (nu*(3-delta^2)/(nu-3) - 3*nu/(nu-2) + 2*(bn*delta)^2))^2
        ##     ## nu
        ##     temp = (1/sz^4) * (3*nu^2/(nu-2)/(nu-4) - 4*(bn*delta)^2*nu*(3-delta^2)/(nu-3) +6*(bn*delta)^2*nu/(nu-2) -3*(bn*delta)^4) - 3
        ##     ans2 = (g2 - temp)^2
        ##     ##
        ##     ans = sqrt(ans1 + ans2)
        ##     return(ans)
        ## }
        ## coef <- nlm(foo, c(alpha,nu), iterlim=5)$est
        ## alpha = coef[1] ; nu = abs(coef[2])
        ##
        ## Update alpha
        foo1 = function(alpha){
            delta = alpha/sqrt(1+alpha^2)
            bn = bn.fun(nu)
            sz = sz.fun(alpha,nu)
            ans = abs(g1 - bn*delta / sz^3 * (nu*(3-delta^2)/(nu-3) - 3*nu/(nu-2) + 2*(bn*delta)^2))
            return(ans)
        }
        alpha <- stats::nlm(foo1,alpha,iterlim=5)$est
        ##
        ## Update nu
        foo2 = function(nu){
            delta = alpha / sqrt(1+alpha^2)
            bn = bn.fun(nu)
            sz = sz.fun(alpha,nu)
            temp = (1/sz^4) * (3*nu^2/(nu-2)/(nu-4) - 4*(bn*delta)^2*nu*(3-delta^2)/(nu-3) +6*(bn*delta)^2*nu/(nu-2) -3*(bn*delta)^4) - 3
            ans = abs(g2 - temp)
            return(ans)
        }
        nu <- stats::nlm(foo2, nu, iterlim=5)$est
        ##
        ## Update dp
        dp <- c(xi,omega,alpha,nu)
        names(dp) = c("xi","omega","alpha","nu")
        return(dp)
    } ## End update.dp
    ##
    ## Repeat updates till convergence
    converged = FALSE
    niter = 0
    ## for (k in 1:maxiter){
    while(!converged){
        niter = niter + 1
        dp = update.dp(dp)
        converged = all(abs(target.moments - dp2moments(dp)) < tol)
        if (niter >= maxiter) break
    }
    ## Prepare output
    ans = list(dp=dp,
               fitted.moments=dp2moments(dp),
               target.moments=target.moments,
               niter=niter, maxiter = maxiter, tol=tol,
               converged=converged)
    return(ans)
}
