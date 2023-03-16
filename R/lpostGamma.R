#' Posterior density function for the non-penalized parameters in an ordgam model
#'
#' @param model An \code{\link{ordgam.object}}
#'
#' @return Log joint posterior density function for the non-penalized regression parameters.
#'
#' @references
#' Lambert, P. and Gressani, O. (2022) Penalty parameter selection and
#' asymmetry corrections to Laplace approximations in Bayesian P-splines models.
#' arXiv:2210.01668.
#'
#' @seealso \code{\link{ordgam}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @examples
#' library(ordgam)
#' data(freehmsData)
#' mod = ordgam(freehms ~ s(eduyrs) + s(age), data=freehmsData, descending=TRUE)
#' print(mod$theta) ## Model regression parameters
#' gam.hat = mod$theta[1:4] ## Non-penalized parameter estimates
#' ordgam::lpost.gamma(mod)(gam.hat)
#'
#' @export
lpost.gamma = function(model){
    ## Function returning the function log p(gamma|data)
    ##    for a given <ordgam> model objectp
    ngamma = with(model, nalpha+nfixed)
    zeta = model$theta ; nzeta = length(zeta)
    ## Special case: No penalized coef --> return p(zeta | D)
    if (ngamma == nzeta) return(model$lpost.fun)
    ## Otherwise: approximate  p(gamma | lambda, D)
    S.zeta = model$Sigma.theta
    ev.zeta = svd(S.zeta)$d ## Eigenvalues of Sigma.zeta
    lambda = model$lambda ## Penalty parameters
    lpost.fun = model$lpost.fun ## Log p(zeta,lambda | data)
    ##
    id1 = 1:ngamma ; id2 = (ngamma+1):nzeta
    S2.1 = S.zeta[id2,id2] - S.zeta[id2,id1] %*% solve(S.zeta[id1,id1],S.zeta[id1,id2])
    ev2.1 = svd(S2.1)$d
    S.1 = S.zeta[id1,id1]
    V = svd(S.1)$U
    ##
    lpgamma = function(gamma){ ## Log p(gamma | lambda,data)
        m2.1 = zeta[id2] + S.zeta[id2,id1] %*% solve(S.zeta[id1,id1], gamma-zeta[id1])
        theta = c(gamma, m2.1)
        ## Evaluate VarCov at (gamma,m2.1)
        S.zeta = attr(lpost.fun(theta,lambda),"Sigma")
        id1 = 1:ngamma ; id2 = (ngamma+1):nzeta
        S2.1 = S.zeta[id2,id2] - S.zeta[id2,id1] %*% solve(S.zeta[id1,id1],S.zeta[id1,id2])
        ev2.1 = svd(S2.1)$d
        ##
        ans = c(lpost.fun(theta,lambda,gradient=FALSE,Hessian=FALSE) + .5*sum(log(ev2.1[ev2.1>1e-6])))
        attr(ans,"Sigma") = S.zeta[id1,id1]
        return(ans)
    }
    ##
    return(lpgamma)
}
