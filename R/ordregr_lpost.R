#' Log-posterior function for a proportional odds model
#'
#' @param y Vector containing the ordinal response (coded using integers in 1:nc).
#' @param nc (optional) Maximum value of \code{y}.
#' @param Xcal Design matrix.
#' @param theta Vector c(alpha,beta) with intercepts <alpha> and regression parameters <beta>.
#' @param descending Logical indicating if the odds of the response taking a value in the upper scale should be preferred over values in the lower scale.
#' @param prior (optional) List given the mean and Prec(ision) of the regression parameters.
#' @param gradient Logical indicating if the gradient of the log-posterior should be computed.
#' @param Hessian Logical indicating if the Hessian of the log-posterior should be computed.
#'
#' @return The log-posterior with the following attributes:
#' \itemize{
#' \item{\code{Salpha} : \verb{ }}{gradient wrt intercepts 'alpha'.}
#' \item{\code{Sbeta} : \verb{ }}{gradient wrt regression parameters 'beta'.}
#' \item{\code{grad} : \verb{ }}{gradient wrt  c(alpha,beta).}
#' \item{\code{Halpha} : \verb{ }}{Hessian wrt intercepts 'alpha'.}
#' \item{\code{Hbeta} : \verb{ }}{Hessian wrt regression parameters 'beta'.}
#' \item{\code{Hba} : \verb{ }}{cross-derivatives (Hessian) submatrix wrt 'alpha' & 'beta'.}
#' \item{\code{Hessian} : \verb{ }}{Hessian wrt c(alpha,beta).}
#' \item{\code{dtheta} : \verb{ }}{step in a Newton-Raphson iteration: solve(-Hessian,grad).}
#' }
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling. doi:10.1177/1471082X231181173. Preprint: arXiv:2210.01668.
#'
#' @export
#'
ordregr_lpost = function(y,nc,Xcal,theta, descending=FALSE,
                        prior=list(mean=NULL,Prec=NULL),
                        gradient=TRUE,Hessian=TRUE){
    ## Check if prior$mean & prior$Prec are both specified for <beta>
    if (is.list(prior)){
        use.prior = !any(unlist(lapply(prior,is.null)))
    } else {
        use.prior = FALSE
    }
    ##
    n = length(y) ## Sample size
    alpha = sort(theta[1:(nc-1)]) ## Intercepts
    if (!is.null(Xcal)){
        beta = tail(theta,-(nc-1))    ## Regression parameters
        nbeta = length(beta)
        mu = c(Xcal %*% beta) ## Part of the linear predictor related to covariates
    } else { ## No regressors, just intercepts
        beta = NULL
        nbeta = 0
        mu = 0
    }
    ## Descending: TRUE if want to model odds of being in the upper response scale
    desc = ifelse(descending,-1,1)
    ## Fij = P(Yi <= j) ; pij = P(Yi=j)
    Fij = pij = matrix(nrow=n,ncol=nc)
    for (j in 1:(nc-1)){
        eta = alpha[j] + desc * mu ## logit(P(Yi<= j))
        Fij[,j] = 1 / (1 + exp(-eta)) ## P(Y <= j)
    }
    Fij[,nc] = 1.0
    pij[,1] = Fij[,1] ## P(Y=1) = P(Y <= 1)
    for (j in 2:(nc-1)) pij[,j] = Fij[,j] - Fij[,j-1] ## P(Y=j) = P(Y <= j) - P(y <= j-1)
    pij[,nc] = 1 - Fij[,nc-1] ## P(Y=nc) = 1 - P(Y <= nc-1)
    ##
    ## Log-lik
    idx = cbind(1:n, y)
    llik = sum(log(pij[idx]))
    ##
    if (use.prior){
        ## lprior.fun:  compute log(p(beta)) with its gradient and precision matrix
        lprior.fun = function(beta,prior){
            ev = svd(prior$Prec)$d
            delta = beta - prior$mean
            Pdelta = c(prior$Prec %*% delta)
            quad = sum(delta * Pdelta)
            ##
            lprior = .5*sum(log(ev[ev>1e-6])) - .5*quad
            grad = -Pdelta
            Prec = prior$Prec
            ##
            ans = list(lprior=lprior, grad=grad, Prec=Prec)
            return(ans)
        }
        ## lpost
        lprior = lprior.fun(beta=beta, prior=prior)
        lpost = llik + lprior$lprior
    }
    ##
    ## Gradient
    if (gradient){
        ## Score (alpha)
        vij = Fij * (1-Fij)
        idx0 = cbind(1:n, y-1)
        temp = matrix(0.,nrow=n,ncol=nc)
        temp[idx] = vij[idx] ; temp[idx0] = -vij[idx0]
        Salpha.i = (1/pij[idx]) * temp[,-nc] ## n x nalpha
        ## Salpha = apply(Salpha.i,2,sum) ## vector(nalpha)
        Salpha = colSums(Salpha.i) ## vector(nalpha)
        attr(llik,"Salpha") = Salpha
        if (nbeta > 0){
            ## Score (beta)
            w.i = (1 + pij - 2*Fij)[idx]
            Sbeta.i = desc * Xcal *  w.i ## n x nbeta
            ## Sbeta = apply(Sbeta.i,2,sum) ## vector(nbeta)
            Sbeta = colSums(Sbeta.i) ## vector(nbeta)
        } else {
            Sbeta = NULL
        }
        attr(llik,"Sbeta") = Sbeta
        ## Score (theta = (alpha,beta))
        grad = c(Salpha,Sbeta)
        attr(llik,"grad") = grad
        ## Score for <lpost>
        if (use.prior){
            attr(lpost,"Salpha") = attr(llik,"Salpha")
            if (nbeta > 0) {
                attr(lpost,"Sbeta") = attr(llik,"Sbeta") + lprior$grad
            } else {
                attr(lpost,"Sbeta") = NULL
            }
            attr(lpost,"grad") = c(attr(lpost,"Salpha"),attr(lpost,"Sbeta"))
        }
        ## Hessian
        if (Hessian){
            ## Hessian (alpha)
            ## ---------------
            ##
            Halpha = -crossprod(Salpha.i)
            ## Correction to diagonal matrix
            zij = (1-2*Fij)*vij
            temp2 = matrix(0.,nrow=n,ncol=nc)
            temp2[idx] = zij[idx] ; temp2[idx0] = -zij[idx0]
            temp2 = (1/pij[idx]) * temp2[,-nc]
            ## diag(Halpha) = diag(Halpha) + apply(temp2,2,sum)
            diag(Halpha) = diag(Halpha) + colSums(temp2)
            attr(llik,"Halpha") = Halpha
            ##
            if (nbeta >0) {
                ## Hessian (beta)
                ## --------------
                term2 = t(apply(pij*(1+pij-2*Fij),1,cumsum))
                w2.i = pij[idx] * w.i - 2*term2[idx]
                Hbeta = t(Xcal) %*% (w2.i * Xcal) ## nbeta x nbeta
                ## ## Approximate Hessian (beta)
                ## Hbeta = crossprod(Salpha.i)
                attr(llik,"Hbeta") = Hbeta
                ##
                ## Hessian (cross derivatives beta-alpha)
                temp3 = matrix(0.,nrow=n,ncol=nc)
                temp3[idx] = vij[idx] ; temp3[idx0] = vij[idx0]
                Hba = - desc * (t(Xcal) %*% temp3)[,-nc]
                attr(llik,"Hba") = Hba ## nalpha x nbeta
                ##
                ## Hessian (theta = (alpha,beta))
                Hes = Matrix::bdiag(Halpha,Hbeta)
                ida = 1:(nc-1) ; idb = nc:ncol(Hes)
                Hes[idb,ida] = Hba
                Hes[ida,idb] = t(Hba)
            } else {
                Hes = Halpha
            }
            attr(llik,"Hessian") = Hes ## ntheta x ntheta
            attr(llik,"Hessian") = as.matrix(Hes) ## ntheta x ntheta
            ##
            dtheta = as.vector(solve(-Hes+1e-6*diag(ncol(Hes)),grad))
            attr(llik,"dtheta") = dtheta
            ## Hessian for <lpost>
            if (use.prior){
                attr(lpost,"Halpha") = attr(llik,"Halpha")
                attr(lpost,"Hbeta") = attr(llik,"Hbeta") - lprior$Prec
                attr(lpost,"Hba") = attr(llik,"Hba")
                ##
                Hes.lpost = attr(llik,"Hessian")
                Hes.lpost[idb,idb] = Hes.lpost[idb,idb] - lprior$Prec
                attr(lpost,"Hessian") = Hes.lpost   ## Hessian with prior
                attr(lpost,"Sigma") = MASS::ginv(-Hes.lpost) ## solve(-Hes.lpost)
                ##
                ## dtheta = as.vector(solve(-Hes.lpost, attr(lpost,"grad")))
                dtheta = as.vector(MASS::ginv(-Hes.lpost) %*% attr(lpost,"grad"))
                attr(lpost,"dtheta") = dtheta
            }
        }
    }
    ##
    if (!use.prior) lpost = llik
    ans = list(llik=llik, lpost=lpost, use.prior=use.prior)
    if (use.prior){
        ans$lpost = lpost
        ans$lprior = lprior
    }
    ##
    return(ans)
} ## End lpost.ordregr
