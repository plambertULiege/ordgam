#' Fit a proportional odds model for ordinal data specified using a formula with (potentially) additive terms (for given fixed penalty parameters).
#'
#' @param formula A model formula
#' @param data A data frame containing a column 'y' with the ordinal response (taking integer values) besides the covariates.
#' @param nc (optional) Number of categories for the ordinal response.
#' @param K Number of B-splines to model each additive term (Default: 10).
#' @param pen.order Penalty order (Default: 2).
#' @param descending Logical indicating if the odds of the response taking a value in the upper scale should be preferred over values in the lower scale.
#' @param Wood.test Logical indicating if the significance test described by S. Wood (Biometrika 2013)  should be used to test for additive term significance. When FALSE, an approximate Chi2 test based on the effective degrees of freedom for each additive term is reported.
#' @param theta0 (Optional) Vector containing starting values for the regression parameters
#' @param lambda0 Vector of penalty parameters for the additive terms (Default: 10 for each additive term).
#' @param ci.level Confidence levels of the computed credible intervals for the regression parameters.
#'
#' @return an object of type \code{\link{ordgam.object}}.
#'
#' @references
#' Lambert, P. and Gressani, 0. (2022) Penalty parameter selection and
#' asymmetry corrections to Laplace approximations in Bayesian P-splines models.
#' arXiv:2210.01668.
#'
#' @seealso \code{\link{ordregr}}, \code{\link{ordgam.object}}.
#'
#' @examples
#' library(ordgam)
#' data(freehmsData)
#' mod = ordgam(freehms ~ gndr + s(eduyrs) + s(age),
#'              data=freehmsData, descending=TRUE)
#' print(mod)
#' plot(mod)
#'
#' @export
ordgam = function(formula, data, nc=NULL, K=10, pen.order=2,
                  descending=TRUE, Wood.test=TRUE,
                  theta0=NULL, lambda0=NULL, ci.level=.95){
    cl <- match.call()
    ## Design matrix, response & other objects for the given <formula> and <data>
    ## --------------------------------------------------------------------------
    regr = DesignFormula(formula=formula, data=data,
                         K=K, pen.order=pen.order, nointercept=TRUE)
    ##
    y = regr$y
    n = nrow(regr$Xcal) ## Number of data entries
    Xcal = regr$Xcal ## Design matrix (omitting intercept columns)
    nbeta = ncol(Xcal) ## Number of regression parameters (omitting intercepts)
    JJ = regr$J ## Number of additive terms
    KK = regr$K ## Number of B-splines per additive term
    nfixed = regr$nfixed ## Number of 'non-spline' regression parameters (excluding intercepts)
    lambda.lab = regr$lambda.lab ## Penalty labels
    Pd.x = regr$Pd.x ## Basic penalty matrix for an additive term
    Z = regr$Z ## Design matrix associated to 'non-spline' regression parameters
    regr.lab = colnames(regr$Xcal) ## Labels of the regression parms
    addregr.lab = regr$additive.lab ## Labels of the additive terms
    qq = ncol(regr$Xcal) ## Total number of regression and spline parameters
    if (is.null(nc)) nc = length(unique(y)) ## Number of categories for the response
    nalpha = nc-1 ## Number of intercepts
    ##
    ## Smoothness prior
    ## ----------------
    if (JJ > 0){
        if (is.null(lambda0)) lambda0 = rep(10,JJ)
        beta0 = rep(0,nbeta)
        Prec.beta = Pcal.fun(nfixed=nfixed, lambda=lambda0, Pd.x)
        prior = list(mean=beta0,Prec=Prec.beta)
    } else {
        prior = list(mean=NULL,Prec=NULL)
    }
    ## Fit Proportional Odds model
    ## ... with penalty for additive terms
    ## -----------------------------------
    obj = ordregr(y=y, nc=NULL, Xcal=Xcal, descending=descending,
                  prior=prior, theta0=theta0, ci.level=ci.level)
    obj$y = y
    obj$regr = regr
    obj$nfixed = nfixed ## Number of non-additive regression parms
    ##
    ## Prepare some extra output
    ## -------------------------
    ## Function testStat implements Wood (2013) Biometrika 100(1), 221-228
    ## Goal: evaluate H0: fj = Xt %*% beta = 0  when  (beta | D) ~ N(p,V)
    ## The type argument specifies the type of truncation to use.
    ## on entry `rank' should be an edf estimate
    ## 0. Default using the fractionally truncated pinv.
    ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
    ## res.df is residual dof used to estimate scale. <=0 implies
    ## fixed scale.
    Tr.test = function(p,Xt,V,edf){
        ans <- testStat(p, Xt, V, min(ncol(Xt), edf), type = 0, res.df = -1)
        ## ans <- mgcv:::testStat(p, Xt, V, min(ncol(Xt), edf), type = 0, res.df = -1)
        return(ans)
    }
    ## End Wood.test
    ##
    ## Calculate EDF, Chi2 and Pval for additive terms
    if (JJ > 0){
        ED = Chi2 = Pval = numeric(JJ)
        for (j in 1:JJ){
            idx = (nalpha + nfixed) + (j-1)*KK + (1:KK)
            beta.j = obj$theta[idx]
            ED[j] = sum(obj$ED.full[idx])
            if (Wood.test){
                ## Begin Wood
                ngrid = 200
                knots.x = regr$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
                pen.order = regr$pen.order
                x.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
                ## Centered B-spline basis
                cB = centeredBasis.gen(x.grid,knots=knots.x,cm=NULL,pen.order)$B
                bele = Tr.test(beta.j,Xt=cB,V=obj$Sigma[idx,idx],ED[j])
                Chi2[j] = bele$stat ; Pval[j] = bele$pval
                ## cat("Wood:",Chi2[j],Pval[j],ED1[j],"\n")
                ## End Wood
                } else {
                    ## "Naive" method
                    Chi2[j] = sum(beta.j * c(solve(obj$Sigma[idx,idx]) %*% beta.j))
                    Pval[j] = 1 - pchisq(Chi2[j],ED[j])
                    ## cat("Naive:",Chi2[j],Pval[j],ED1[j],"\n")
                    ## End Naive
                }
            }
        ED = cbind(ED,Chi2,Pval)
        rownames(ED) = paste("f(",regr$additive.lab,")",sep="")
        if (Wood.test){
            colnames(ED) = c("edf","Tr","Pval")
        } else {
            colnames(ED) = c("edf","Chi2","Pval")
        }
        obj$ED=ED
        obj$Wood.test = Wood.test
    }
    ## Function computing the log-posterior and its attributes
    ##   for given <theta> & <lambda>
    lpost.fun = function(theta,lambda,gradient=TRUE,Hessian=TRUE){
        if (JJ >0){
            beta0 = rep(0,nbeta)
            Prec.beta = Pcal.fun(nfixed=nfixed, lambda=lambda, Pd.x)
            prior = list(mean=beta0,Prec=Prec.beta)
        } else {
            prior = NULL
        }
        ans = ordregr_lpost(y=y, nc=nc, Xcal=Xcal, theta=theta, descending=descending,
                            prior=prior, gradient=gradient, Hessian=Hessian)$lpost
        return(ans)
    }
    ##
    obj$lpost.fun = lpost.fun
    obj$lambda = lambda0
    ##
    obj$call = cl
    obj$formula = formula
    ans = obj
    class(ans) = c("ordregr","ordgam")
    return(ans)
}
## End ordgam
