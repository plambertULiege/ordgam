#' Fit of an additive proportional odds model for ordinal data using Laplace approximations and P-splines
#'
#' @param formula A model formula
#' @param data A data frame containing a column 'y' with the ordinal response (taking integer values) besides the covariates.
#' @param nc (optional) Number of categories for the ordinal response.
#' @param K Number of B-splines to model each additive term (Default: 10).
#' @param pen.order Penalty order (Default: 2).
#' @param descending Logical indicating if the odds of the response taking a value in the upper scale should be preferred over values in the lower scale (Default: TRUE).
#' @param select.lambda Logical indicating if the penalty parameters should be tuned (Default: TRUE).
#' @param lambda.family Prior for <lambda>. Possible choices are "none", "dgamma", "BetaPrime" or "myprior" for a user specified function for the prior of <lambda>.
#' @param lambda.optimizer Algorithm used to maximize p(lambda|data). Possible choices are "nlminb","ucminf","nlm","LevMarq" (Default: "nlminb").
#' @param lprior.lambda Log of the prior density for a <lambda> component if \code{lambda.family} set to "myprior".
#' @param theta0 (Optional) Vector containing starting values for the regression parameters.
#' @param lambda0 Vector of penalty parameters for the additive terms (Default: 10 for each additive term).
#' @param ci.level Confidence levels of the computed credible intervals for the regression parameters.
#' @param verbose Verbose mode (logical)
#'
#' @return an object of type \code{\link{ordgam.object}}.
#'
#' @references
#' Lambert, P. and Gressani, O. (2022) Penalty parameter selection and
#' asymmetry corrections to Laplace approximations in Bayesian P-splines models.
#' arXiv:2210.01668.
#'
#' @seealso \code{\link{ordregr}}, \code{\link{ordgam.object}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
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
                  descending=TRUE,
                  select.lambda=TRUE, lambda.family="dgamma",
                  lambda.optimizer="nlminb",
                  lprior.lambda = function(x) dgamma(x,1,1e-4,log=TRUE),
                  theta0=NULL, lambda0=NULL, ci.level=.95, verbose=FALSE){
    ptm <- proc.time() ## Start timer
    cl <- match.call()
    match.arg(lambda.family, c("none","dgamma","BetaPrime","myprior"))
    match.arg(lambda.optimizer, c("nlminb","ucminf","nlm","LevMarq"))
    ## Design matrix, response & other objects for the given <formula> and <data>
    ## --------------------------------------------------------------------------
    regr = DesignFormula(formula=formula, data=data,
                         K=K, pen.order=pen.order, nointercept=TRUE)
    ##
    y = regr$y
    n = nrow(regr$Xcal) ## Number of data entries
    Xcal = regr$Xcal ## Design matrix (omitting intercept columns)
    nbeta = ncol(Xcal) ## Number of regression parameters (omitting intercepts)
    if (is.null(nbeta)) nbeta = 0
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
        if (is.null(lambda0)) lambda0 = rep(100,JJ)
        lambda = lambda0
        names(lambda0) = lambda.lab
        beta0 = rep(0,nbeta)
        Prec.beta = Pcal.fun(nfixed=nfixed, lambda=lambda, Pd.x)
        ## diag(Prec.beta)[1:nfixed,1:nfixed] = diag(Prec.beta)[1:nfixed,1:nfixed] + 1e-6
        prior = list(mean=beta0,Prec=Prec.beta)
    } else {
        prior = list(mean=NULL,Prec=NULL)
    }
    ## Fit Proportional Odds model
    ## ... with penalty for additive terms
    ## -----------------------------------
    ##
    select.lambda = select.lambda & (JJ > 0) ## Only make sense if additive terms...
    if (select.lambda){
        ## log p(lambda) (i.e. log prior for a penalty parameter
        lprior.lambda = switch(lambda.family,
                            "none" = function(lambda) 0,
                            "dgamma" = function(lambda) dgamma(lambda,1,1e-4,log=TRUE),
                            "BetaPrime" = function(lambda) -log(pi)-.5*log(lambda)-log(1+lambda),
                            "myprior" = lprior.lambda)
        ##
        ## Straightforward optimization to select <lambda>
        ## -----------------------------------------------
        nlambda = JJ
        ## Loss function to select <lambda>:
        ##     loss.fn(nu) = -log p(lambda=exp(nu) | data)
        loglambda.loss = function(loglambda,data){
            lambda = exp(loglambda)
            Prec.beta = Pcal.fun(nfixed=nfixed, lambda=lambda, Pd.x)
            prior = list(mean=beta0,Prec=Prec.beta)
            temp = ordregr(y=y, nc=NULL, Xcal=Xcal, descending=descending,
                           prior=prior, theta0=theta0, ci.level=ci.level)
            ##
            lmargllik = temp$levidence ## Log marginal likelihood
            lprior = sum(sapply(lambda, lprior.lambda))
            ans = -lmargllik - lprior
            return(ans)
        }
        ##
        switch(lambda.optimizer,
               "LevMarq" = {
                   ## Levenberg-Marquardt algorithm to select <lambda>
                   obj.ml = marqLevAlg::mla(b=log(lambda0),fn=loglambda.loss,data=data)
                                             ## nproc=2, clustertype="FORK")
                   lambda = exp(obj.ml$b) ## Selected <lambda> --> lambda.hat
               },
               "nlminb" = {
                 ## <nlminb> to select <lambda>
                 obj.ml = nlminb(start=log(lambda0),objective=loglambda.loss,
                                 lower=rep(0,JJ),upper=(rep(10,JJ)))
                 lambda = exp(obj.ml$par) ## Selected <lambda> --> lambda.hat
               },
               "ucminf" = {
                 ## <ucminf> to select <lambda>
                 obj.ml = ucminf::ucminf(par=log(lambda0),fn=loglambda.loss)
                 lambda = exp(obj.ml$par) ## Selected <lambda> --> lambda.hat
               },
               "nlm" = {
                 ## <nlm> to select <lambda>
                 obj.ml = nlm(f=loglambda.loss,p=log(lambda0))
                 lambda = exp(obj.ml$est) ## Selected <lambda> --> lambda.hat
               }
        )
        ##
        names(lambda) = lambda.lab
        if (verbose) cat("lambda:",lambda,"\n")
        ## Update smoothness prior accordingly
        Prec.beta = Pcal.fun(nfixed=nfixed, lambda=lambda, Pd.x)
        prior = list(mean=beta0,Prec=Prec.beta)
        ##
        if (lambda.family != "none"){
            ## nu = log(lambda) posterior
            ##  (CAREFUL: log(Jacobian) added as compared to loglambda.loss)
            ## --------------------------
            nu.lpost = function(nu,data){
                lambda = exp(nu)
                ans = -loglambda.loss(loglambda=nu,data=data) + sum(nu)
                return(ans)
            }
            switch(lambda.optimizer,
              "LevMarq" = {
                  ## Levenberg-Marquardt algorithm to get <nu.hat> from p(nu=loglambda) | data)
                  obj.nu = marqLevAlg::mla(b=log(lambda0),fn=nu.lpost,
                                           minimize=FALSE, data=data)
                                           ## nproc=2, clustertype="FORK")
                  nu = obj.nu$b ## Posterior mode of <nu> = log(lambda) (--> differs from log(lambda.hat) !!)
                  if (verbose) cat("nu:",nu,"\n")
                  ## Reconstruct Var-Cov from its upper triangular output from maqLevAlg
                  nlambda = length(lambda0)
                  mat = matrix(0,nrow=nlambda,ncol=nlambda)
                  mat[!lower.tri(mat)] = obj.nu$v ## Upper triangular + diag of the Hessian
                  mat[lower.tri(mat)] = mat[upper.tri(mat)] ## Lower triang. w/o diagonal
                  rownames(mat) = colnames(mat) = lambda.lab
                  V.nu = mat
              },
             "nlminb" = {
               ## nlminb to get nu.hat from p(nu=loglambda) | data)
                 obj.nu = nlminb(start=log(lambda0),
                                 function(x) -nu.lpost(x,data=data),
                                 lower=rep(0,JJ),upper=rep(10,JJ))
               nu = obj.nu$par
               if (verbose) cat("nu:",nu,"\n")
               V.nu = numDeriv::hessian(function(x) -nu.lpost(x,data=data), nu)
             },
            "ucminf" = {
              ## <ucminf> to get nu.hat from p(nu=loglambda) | data)
                obj.nu = ucminf::ucminf(par=log(lambda0),
                                        fn=function(x) -nu.lpost(x,data=data),
                                        hessian=TRUE)
              nu = obj.nu$par
              if (verbose) cat("nu:",nu,"\n")
              V.nu = obj.nu$hessian
            },
            "nlm" = {
              ## nlm to get nu.hat from p(nu=loglambda) | data)
              obj.nu = nlm(function(x) -nu.lpost(x), p=log(lambda0), hessian=TRUE)
              nu = obj.nu$est
              if (verbose) cat("nu:",nu,"\n")
              V.nu = obj.nu$hessian
            }
            )
            ##
            se.nu = sqrt(diag(V.nu))
            names(se.nu) = lambda.lab
            ## Approximate marginal posterior for nu = log(lambda[k])
            nu.dp = list()
            for (k in 1:nlambda){
                gg = function(x){
                    arg = nu
                    arg[k] = x
                    ans = nu.lpost(arg,data=data)
                    return(ans)
                }
                # ## Method 1
                # xr = nu[k] + c(-1,1) * 5 * se.nu[k]
                # xg = seq(xr[1],xr[2],length=100)
                # ggrid = sapply(xg, gg)
                # yg = exp((ggrid-max(ggrid)))
                # dp = stmatch(xg,yg)$dp ## snmatch(xg,yg)$dp
                ##
                ## Method 2
                xr = nu[k] + c(-1,1) * 3 * se.nu[k]
                xg = seq(xr[1],xr[2],length=10)
                ggrid = sapply(xg, gg)
                lyg = ggrid-max(ggrid)
                dp = STapprox(xg,lyg)$dp ## Match a skew-t to log p(nu[k]|Data)
                ##
                nu.dp[[k]] = dp
            }
            names(nu.dp) = lambda.lab
        }
    } ## Endif (select.lambda)
    ##
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
        ED = Chi2 = Tr = Pval.Chi2 = Pval.Tr = numeric(JJ)
        for (j in 1:JJ){
            idx = (nalpha + nfixed) + (j-1)*KK + (1:KK)
            beta.j = obj$theta[idx]
            ED[j] = sum(obj$ED.full[idx])
            ## Wood Chi2 test
            ngrid = 200
            knots.x = regr$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
            pen.order = regr$pen.order
            x.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
            ## Centered B-spline basis
            cB = centeredBasis.gen(x.grid,knots=knots.x,cm=NULL,pen.order)$B
            bele = Tr.test(beta.j,Xt=cB,V=obj$Sigma[idx,idx],ED[j])
            Tr[j] = bele$stat
            Pval.Tr[j] = bele$pval
            ## cat("Wood:",Chi2[j],Pval[j],ED1[j],"\n")
            ## End Wood
            ##
            ## Straight Chi2 test
            Chi2[j] = sum(beta.j * c(solve(obj$Sigma[idx,idx]) %*% beta.j))
            Pval.Chi2[j] = 1 - pchisq(Chi2[j],ED[j])
            ## cat("Chi2:",Chi2[j],Pval[j],ED1[j],"\n")
            ## End Chi2
        }
        ED.Chi2 = cbind(ED,Chi2,Pval.Chi2)
        ED.Tr = cbind(ED,Tr,Pval.Tr)
        rownames(ED.Chi2) = rownames(ED.Tr) = paste("f(",regr$additive.lab,")",sep="")
        colnames(ED.Chi2) = c("edf","Chi2","Pval")
        colnames(ED.Tr) = c("edf","Tr","Pval")
        obj$ED.Chi2 = ED.Chi2
        obj$ED.Tr = ED.Tr
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
    if (JJ > 0){
        obj$lambda0 = lambda0
        obj$lambda = lambda
        obj$select.lambda = select.lambda
        obj$lambda.family = lambda.family
        if (select.lambda){
            obj$lprior.lambda = lprior.lambda
            obj$loglambda.loss = loglambda.loss
            if (lambda.family != "none"){
                obj$nu.lpost = nu.lpost
                obj$nu.hat= nu
                obj$V.nu = V.nu
                obj$se.nu = se.nu
                obj$nu.dp = nu.dp
            } else {
                obj$nu.dp = NA
            }
        }
    }
    ## ## Nbr of non-penalized parameters in the additive part
    ## nparms = with(obj, sum(ED.full))
    ##
    ## Nbr of non-penalized spline parameters (in the additive part)
    nparms = with(obj$regr, (pen.order-1) * J)
    obj$levidence = obj$levidence + .5 * nparms * log(2*pi)
    ##
    obj$call = cl
    obj$formula = formula
    obj$elapsed.time <- (proc.time()-ptm)[1] ## Elapsed time
    ans = obj
    class(ans) = c("ordregr","ordgam")
    return(ans)
}
## End ordgam
