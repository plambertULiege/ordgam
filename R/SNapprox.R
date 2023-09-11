#' Skew-Normal approximation to a density evaluated on a sparse grid
#'
#' @param x Vector containing a grid of values on the density support and covering the posterior mode.
#' @param lfx Log density values on the grid \code{x} (possibly up to an additive constant)
#'
#' @return A list containing
#' \itemize{
#' \item{\code{dp} : \verb{ }}{Parameters of the approximating skew-Normal density.}
#' \item{\code{fitted.moments} : \verb{ }}{Mean, variance, skewness, kurtosis of the approximating skew-Normal.}
#' }
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling. doi:10.1177/1471082X231181173. Preprint: arXiv:2210.01668.
#'
#' @seealso \code{\link{STapprox}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @examples
#' library(ordgam)
#'
#' ## Density to be approximated by a Skew-Normal
#' dtarget = function(x) dgamma(x,10,2)
#' curve(dtarget(x),0,15,lwd=2,ylab="Density")
#'
#' ## Values of the target density on a sparse grid
#' ngrid = 6 ## Sparse grid size
#' xgrid = seq(2,8,length=ngrid) ## Grid
#' lfx = log(dtarget(xgrid)) ## Log values
#'
#' ## Skew-Normal approximation
#' dp = ordgam::SNapprox(xgrid,lfx)$dp
#' curve(sn::dsn(x,dp=dp),add=TRUE,lwd=2,lty=2,col=2)
#' points(xgrid,exp(lfx))
#' legend("topright",legend=c("Target density","Skew-Normal approx."),
#'        col=1:2,lty=1:2,lwd=2,bty="n")#' library(ordgam)
#' @export
#'
SNapprox = function(x,lfx){
  dp0 = c(x[which.max(lfx)], sqrt(-.5/lm(lfx ~ x + I(x^2))$coef[3]), 0)
  ##
  lossfun = function(dp){
    dp[2] = abs(dp[2])
    ## log of the approximating skew-t density
    lf = sn::dsn(x,dp=dp,log=TRUE) ; lf = lf-max(lf)
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
  ##
  ans = list(dp=dp.hat,
             fitted.moments=dp2moments(dp.hat))
}
