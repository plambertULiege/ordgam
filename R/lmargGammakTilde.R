## Reparametrize p(gamma| lambda,data) --> p(gamma.tilde| lambda,data)
## ----------------------------------
#' Marginal posterior density function for a remapped non-penalized parameter in an ordgam model
#'
#' @param gamtk Remapped parameter value at which the marginal log posterior density for <gamma.tilde[k]> must be evaluated.
#' @param k Targetted component in the vector of remapped non-penalized parameters <gamma.tilde>.
#' @param model An \code{\link{ordgam.object}}.
#'
#' @return Log of p(gamma.tilde[k] | lambda,data)
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling. doi:10.1177/1471082X231181173. Preprint: arXiv:2210.01668.
#'
#' @seealso \code{\link{ordgam}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @examples
#' library(ordgam)
#' data(freehmsData)
#' mod = ordgam(freehms ~ s(eduyrs) + s(age), data=freehmsData, descending=TRUE,
#'              lambda0=c(192,18),select.lambda=FALSE)
#' ngamma = with(mod, nalpha+nfixed) ## Number of non-penalized parms
#' k = 1 ## Focus on gamma.tilde[1]
#' x.grid = seq(-4,4,length=7) ## Grid of values for gamma.tilde[k]
#' lfy.grid = ordgam::lmarg.gammaTilde(x.grid,k=k,mod) ## log p(gamma.tilde[k] | D) on the grid
#' gamt.ST = ordgam::STapprox(x.grid,lfy.grid)$dp ## Approximate using a skew-t
#' ## Plot the estimated marginal posterior for <gamma.tilde[k]>
#' xlab = bquote(tilde(gamma)[.(k)])
#' ylab = bquote(p(tilde(gamma)[.(k)]~ "|"~lambda~","~D))
#' xlim = sn::qst(c(.0001,.9999),dp=gamt.ST)
#' curve(sn::dst(x,dp=gamt.ST),xlim=xlim,
#'       xlab=xlab,ylab=ylab,col="blue",lwd=2,lty=1)
#'
#' @export
lmarg.gammaTilde = function(gamtk, k, model){
    ## Function giving log p(gamma.tilde[k]|lambda,data)
    ## NOTE: the argument <gamtk> can be a vector of values
    ngamma = with(model, nalpha+nfixed)
    log.joint = lpost.gamma(model)
    ##
    gam.hat = model$theta[1:ngamma]
    Sig = model$Sigma.theta[1:ngamma,1:ngamma]
    ##
    sv = svd(Sig)
    V = sv$u
    Vk = V[,k]
    omega = sv$d
    ##
    J = length(gamtk) ; ans = numeric(J)
    for (j in 1:J){ ## Loop over values of gam.tilde[k] in <gamtk>
        ## Approximate p(gamma.tilde[k] | lambda,data)
        gam.j = gam.hat + gamtk[j] * sqrt(omega[k]) * Vk
        ans[j] = c(log.joint(gam.j))
    }
    return(ans)
}
