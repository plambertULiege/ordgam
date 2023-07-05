#' Compute the additive terms estimated using an 'ordgam' model
#'
#' @param obj.ordgam An object of class 'ordgam'.
#' @param ngrid Number of grid points where the additive terms are computed.
#' @param ci.level Credibility level for the pointwise credible region for the additive terms
#'
#' @return a list containing:
#' \itemize{
#' \item{\code{nalpha} : \verb{ }}{number of intercepts in the proportional odds model.}
#' \item{\code{nfixed} : \verb{ }}{number of non-penalized regression parameters in 'beta'.}
#' \item{\code{J} : \verb{ }}{number of additive terms.}
#' \item{\code{additive.lab} : \verb{ }}{labels of the additive terms.}
#' \item{\code{K} : \verb{ }}{number of spline parameters to specify an additive term.}
#' \item{\code{knots} : \verb{ }}{list of length J containing the knots for the B-spline basis associated to a given additive term.}
#' \item{\code{f.grid} : \verb{ }}{list of length J with, for each additive term, a list of length 2 with 'x': a vector of grid values for the covariate ; 'y.mat': a matrix with 3 columns (est,low,up) giving the additive term and its pointwise credible region}
#' \item{\code{f} : \verb{ }}{a list of length J with, for each additive term <x>, a list with f$x: a function computing the additive term f(x) for a given covariate value 'x' ; attributes(f$x): support, label, range.}
#' \item{\code{f.se} : \verb{ }}{a list of length J with, for each additive term <x>, a list with f.se$x: a function computing the s.e. of f(x) for a given covariate value 'x' ; attributes(f.se$x): support, label, range}
#' }
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling (in press). Preprint: arXiv:2210.01668.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @examples
#' library(ordgam)
#' data(freehmsData)
#' mod = ordgam(freehms ~ gndr + s(eduyrs) + s(age),
#'              data=freehmsData, descending=TRUE)
#' obj = ordgam_additive(mod)
#' names(obj)
#' with(obj$f.grid$age,
#'       matplot(x, y.mat, lty=c(1,2,2),type="l",col=1,
#'               xlab="Age", ylab="f(Age)"))
#'
#' @export
ordgam_additive <- function(obj.ordgam,ngrid=300, ci.level=.95){
    obj = obj.ordgam
    nfixed = obj$regr$nfixed ## Number of non-additive regression parms in <beta>
    nalpha = obj$nalpha ## Number of intercepts
    J = obj$regr$J ## Number of additive terms
    K = obj$regr$K ## Number of centered B-splines per additive term
    f = f.se = list() ## Fitted additive terms
    alpha = 1 - ci.level ; z.alpha = qnorm(1-.5*alpha)
    ans = NULL
    ##
    rangetol = function(x,tol=.01){
        temp = range(x) ; df = diff(temp)
        ans = c(temp[1]-.5*tol*df,temp[2]+.5*tol*df)
        return(ans)
    }
    ## Additive part
    ans$nalpha = obj$nalpha ; ans$nfixed=nfixed ; ans$J = J
    if (J > 0){
        add.lab = obj$regr$additive.lab
        ans$additive.lab = add.lab
        ans$K=K ; ans$knots = obj$regr$knots.x
        mat = matrix(nrow=ngrid,ncol=J) ; colnames(mat) = add.lab
        Sigma = solve(-obj$Hessian) ## Sigma(theta)
        ## x.grid = f.grid = f.grid.se = f.grid.low = f.grid.up = mat
        f.grid = list()
        ##
        for (j in 1:J){
            idx = nalpha + nfixed + (j-1)*K + (1:K)
            beta.j = obj$theta.mat[idx,"est"] ## Centered B-splines coefs for jth additive term
            knots.x = obj$regr$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
            pen.order = obj$regr$pen.order
            x.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
            cB = centeredBasis.gen(x.grid,knots=knots.x,cm=NULL,pen.order)$B ## Centered B-spline basis
            y.grid = c(cB %*% beta.j)
            y.grid.se = sqrt(diag(cB %*% (Sigma[idx,idx]%*%t(cB))))
            ylow = y.grid - z.alpha*y.grid.se
            yup  = y.grid + z.alpha*y.grid.se
            ##
            f.grid[[add.lab[j]]]$x = x.grid
            f.grid[[add.lab[j]]]$y.mat = cbind(est=y.grid,low=ylow,up=yup)
            ##
            f[[add.lab[j]]]    = splinefun(x.grid, y.grid)
            f.se[[add.lab[j]]] = splinefun(x.grid, y.grid.se)
            attr(f[[add.lab[j]]],"support") = attr(f.se[[add.lab[j]]],"support") = c(xL,xU)
            attr(f[[add.lab[j]]],"label") = attr(f.se[[add.lab[j]]],"label") = add.lab[j]
            attr(f[[add.lab[j]]],"range") = rangetol(y.grid)
        }
        ans$f.grid = f.grid
        ans$f = f ; ans$f.se = f.se ## ; ans$ED1 = obj$fit$ED1
    }
    ##
    return(ans)
}
