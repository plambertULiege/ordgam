#' Fit a proportional odds model for ordinal data
#'
#' @param y Vector containing the ordinal response (coded using integers in 1:nc).
#' @param nc (optional) Maximum value of \code{y}.
#' @param Xcal Design matrix (excluding intercept columns).
#' @param descending Logical indicating if the odds of the response taking a value in the upper scale should be preferred over values in the lower scale.
#' @param prior (optional) List giving the 'mean' and 'Prec'(ision) of the regression parameters.
#' @param theta0 (Optional) Vector containing starting values for the regression parameters.
#' @param ci.level Confidence levels of the computed credible intervals for the regression parameters.
#'
#' @return An object of class \link{ordregr.object}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling (in press). Preprint: arXiv:2210.01668.
#'
#' @seealso \code{\link{ordgam}}, \code{\link{ordregr.object}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @examples
#' library(ordgam)
#' data(freehmsData)
#' Xcal = with(freehmsData, cbind(gndr,eduyrs,age))
#' mod = ordregr(y=freehmsData$freehms, Xcal=Xcal, descending=TRUE)
#' print(mod)
#'
#' @export
#'
ordregr = function(y, nc=NULL, Xcal=Xcal, descending=FALSE,
                  prior=list(mean=NULL,Prec=NULL), theta0=NULL, ci.level=.95){
    cl <- match.call()
    n = length(y)
    if (is.null(nc)) nc = length(unique(y))
    nalpha = nc-1 ; nbeta = ncol(Xcal)
    if (is.null(nbeta)) nbeta = 0
    ## Initial values
    if (is.null(theta0)){
        tab = table(y)
        alpha0 = qlogis(head(cumsum(tab/sum(tab)),-1))
        beta0 = rep(0,nbeta)
        theta0 = c(alpha0,beta0)
        ntheta = length(theta0)
    } else {
        ntheta = length(theta0)
        if (ntheta != (nalpha + nbeta)){
            cat("<theta0> should be of length",ntheta,"\n")
            return(NULL)
        }
        alpha0 = theta0[1:(nc-1)]
        beta0 = tail(theta0,nbeta)
        if (!all(diff(alpha0) > 0)){
            cat("The first",nalpha,"entries in <theta0> should be increasing !\n")
            return(NULL)
        }
    }
    names(theta0) = c(1:nalpha,colnames(Xcal))
    ## L2-norm
    L2norm = function(x) sqrt(sum(x^2))
    ## Generic Newton-Raphson algorithm
    ## --------------------------------
    NewtonRaphson = function(g, theta, tol=1e-2, itermax=20, verbose=FALSE){
        theta.cur = theta
        obj.cur = g(theta.cur,Dtheta=TRUE)
        g.start = obj.cur$g ## Function value at the start
        ok = (L2norm(obj.cur$grad) < tol) ## Stopping rule
        ## ok = all(abs(obj.cur$grad) < tol) ## Stopping rule
        iter = 0
        while(!ok){
            iter = iter + 1
            theta.cur = obj.cur$theta
            dtheta = c(obj.cur$dtheta)
            step = 1 ; nrep = 0
            repeat { ## Repeat step-halving directly till improve target function
                nrep = nrep + 1
                if (nrep > 5) break
                theta.prop = theta.cur + step*dtheta ## Update.theta
                obj.prop = g(theta.prop,Dtheta=TRUE)
                if (obj.prop$g >= obj.cur$g) break
                step = .5*step
            }
            obj.cur = obj.prop
            ok = (L2norm(obj.cur$grad) < tol) ## Stopping rule
            ## ok = all(abs(obj.cur$grad) < tol) ## Stopping rule
            if (iter > itermax) break
        }
        if (verbose) cat(obj.cur$g," (niter = ",iter,") - grad = ",L2norm(obj.cur$grad),"\n",sep="")
        ## if (verbose) cat(obj.cur$g," (niter = ",iter,")\n",sep="")
        ans = list(val=obj.cur$g, val.start=g.start, theta=obj.cur$theta, grad=obj.cur$grad, Hessian=obj.cur$Hessian, iter=iter)
        return(ans)
    } ## End NewtonRaphson
    ##
    ## Function to be optimized
    g.fun = function(theta, Dtheta=TRUE){
        temp = ordregr_lpost(y=y, nc=nc, Xcal=Xcal, theta=theta, descending=descending,
                            prior=prior, gradient=Dtheta, Hessian=Dtheta)$lpost
        ans = list(g=c(temp), theta=theta,
                   dtheta=attr(temp,"dtheta"), grad=attr(temp,"grad"),
                   Hessian=attr(temp,"Hessian"))
        return(ans)
    }
    ## N-R for <theta>
    grad.tol = 1e-3
    obj.NR = NewtonRaphson(g=g.fun,theta=theta0,tol=grad.tol)
    ## Compute Hessian with and without prior
    obj.llik = ordregr_lpost(y=y, nc=nc, Xcal=Xcal, theta=obj.NR$theta,
                             descending=descending,
                             gradient=TRUE, Hessian=TRUE)$llik
    obj.NR$llik = c(obj.llik)
    obj.NR$Hessian0 = attr(obj.llik,"Hessian") ## Hessian without prior
    obj.NR$Sigma.theta = with(obj.NR, MASS::ginv(-Hessian)) ## Var-Cov with prior
    obj.NR$ED.full = with(obj.NR, rowSums(t(Sigma.theta) * (-Hessian0))) ## Effective dim
    names(obj.NR$ED.full) = names(theta0)
    ##
    ## Report parameter estimates with their se's, z-score, Pval
    fun = function(est,se,ci.level=.95){
        alpha = 1 - ci.level
        z.alpha = qnorm(1-.5*alpha)
        mat = cbind(est=est,se=se,
                    low=est-z.alpha*se, up=est+z.alpha*se,
                    "Z"=est/se,
                    "Pval"=1-pchisq((est/se)^2,1))
        attr(mat,"ci.level") = ci.level
        return(mat)
    } ## End fun
    ##
    ans = obj.NR
    ans$se.theta = sqrt(diag(MASS::ginv(-obj.NR$Hessian)))
    ans$theta.mat = with(ans, fun(est=theta,se=se.theta,ci.level=ci.level))
    ans$nc = nc ; ans$nalpha = nalpha ; ans$nbeta = nbeta ; ans$nfixed=nbeta
    ans$ci.level = ci.level
    ##
    ans$n=n
    ans$call=cl
    ans$descending=descending
    ## Check if prior$mean & prior$Prec are both specified for <beta>
    use.prior = !any(unlist(lapply(prior,is.null)))
    ans$use.prior = use.prior
    ##
    ev = svd(-obj.NR$Hessian)$d
    ans$lpost = ans$val
    ans$levidence = ans$lpost -.5*sum(log(ev[ev>1-6])) ## Log marginal likelihood
    ans$AIC  = -2*ans$llik + 2*sum(ans$ED.full)
    ans$BIC  = -2*ans$llik + sum(ans$ED.full)*log(n)
    ##
    class(ans) = "ordregr"
    return(ans)
} ## End ordregr
