#' Skew-t approximation to a density observed on a sparse grid
#'
#' @param x Vector containing a grid of values on the density support and covering the posterior mode.
#' @param lfx Log density values on the grid \code{x} (possibly up to an additive constant)
#'
#' @return A list containing
#' \itemize{
#' \item{\code{dp} : \verb{ }}{Parameters of the approximating skew-t density.}
#' \item{\code{fitted.moments} : \verb{ }}{Mean, variance, skewness, kurtosis of the approximating skew-t.}
#' }
#'
#'
#' @references
#' Lambert, P. and Gressani, O. (2022) Penalty parameter selection and
#' asymmetry corrections to Laplace approximations in Bayesian P-splines models.
#' arXiv:2210.01668.
#'
#' @seealso \code{\link{SNapprox}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @examples
#' library(ordgam)
#'
#' ## Density to be approximated by a Skew-t
#' dtarget = function(x) dgamma(x,10,2)
#' curve(dtarget(x),0,15,lwd=2,ylab="Density")
#'
#' ## Values of the target density on a sparse grid
#' ngrid = 6 ## Sparse grid size
#' xgrid = seq(2,8,length=ngrid) ## Grid
#' lfx = log(dtarget(xgrid)) ## Log values
#'
#' ## Skew-t approximation
#' dp = ordgam::STapprox(xgrid,lfx)$dp
#' curve(sn::dst(x,dp=dp),add=TRUE,lwd=2,lty=2,col=2)
#' points(xgrid,exp(lfx))
#' legend("topright",legend=c("Target density","Skew-t approx."),
#'        col=1:2,lty=1:2,lwd=2,bty="n")
#'
#' @export
STapprox = function(x,lfx){
  dp0 = c(x[which.max(lfx)], sqrt(-.5/lm(lfx ~ x + I(x^2))$coef[3]), 0, 20)
  ##
  lossfun = function(dp){
    dp[2] = abs(dp[2])
    ## log of the approximating skew-t density
    lf = sn::dst(x,dp=dp,log=TRUE) ; lf = lf-max(lf)
    ## log of the target density
    lg = lfx ; lg = lg-max(lg)
    ## Symmetric KL divergence
    ans = sum(exp(lf)*(lf-lg) + exp(lg)*(lg-lf))
    return(ans)
  }
  dp.hat = unname(nlminb(dp0,lossfun)$par)
  dp.hat[2] = abs(dp.hat[2])
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
  ans = list(dp=dp.hat,
             fitted.moments=dp2moments(dp.hat))
}
