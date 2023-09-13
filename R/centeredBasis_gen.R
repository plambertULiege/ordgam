## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018)
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' Generation of a recentered B-spline basis matrix in additive models
#' @description Generation of a B-spline basis matrix with recentered columns
#'  to handle the identifiability constraint in additive models. See Wood (CRC Press 2017, pp. 175-176) for more details.
#' @param x vector of values where to compute the "recentered" B-spline basis
#' @param knots vector of knots (that should cover the values in <x>)
#' @param cm (Optional) values subtracted from each column of the original B-spline matrix
#' @param pen.order penalty order for the B-spline parameters (Default: 2)
#' @param verbose verbose indicator (Default: TRUE)
#'
#' @return List containing
#' \itemize{
#' \item{\code{B} : \verb{ }}{centered B-spline matrix (with columns recentered to have mean 0 over equi-spaced x values on the range of the knots).}
#' \item{\code{Dd} : \verb{ }}{difference matrix (of order <pen.order>) for the associated centered B-spline matrix.}
#' \item{\code{Pd} : \verb{ }}{penalty matrix (of order <pen.order>) for the associated centered B-spline matrix.}
#' \item{\code{K} : \verb{ }}{number of centered B-splines in the basis.}
#' \item{\code{cm} : \verb{ }}{values subtracted from each column of the original B-spline matrix. By default, this is a vector containing the mean of each column in the original B-spline matrix.}
#'}
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling. <doi:10.1177/1471082X231181173>. Preprint: <arXiv:2210.01668>.
#'
#' @export
#'
#' @examples
#' x = seq(0,1,by=.01)
#' knots = seq(0,1,length=5)
#' obj = centeredBasis.gen(x,knots)
#' matplot(x,obj$B,type="l",ylab="Centered B-splines")
#' colMeans(obj$B)
#'
centeredBasis.gen = function(x,knots,cm=NULL,pen.order=2,verbose=TRUE){
  if ((max(x)>max(knots))|(min(x)<min(knots))){
    if (verbose) cat("The knots do no cover the values of x !!\n")
    return(NULL)
  }
  ##
  knots.x = knots ; pen.order.x = pen.order
  ##
  temp = cubicBsplines::Bsplines(x, knots.x)
  idx = 1:(ncol(temp)-1)
  B <- temp[,idx,drop=FALSE] ## Reduce the number of B-splines by 1 !!
  if (is.null(cm)){
    ## Consider a B-spline matrix over a dense equi-spaced set of values on range(knots)
    ##  to compute a recentering value for each column of the matrix
    B.ref = cubicBsplines::Bsplines(seq(min(knots),max(knots),length=100), knots.x)[,idx]
    cm = colMeans(B.ref)
  }
  B = sweep(B,2,cm) ## Substract column mean from each column
  Dd = diff(diag(ncol(temp)),dif=pen.order.x)[,idx] ## Penalty order <---- Important to remove last column !!
  Pd = t(Dd) %*% Dd
  return(list(B=B,Dd=Dd,Pd=Pd,K=ncol(B),cm=cm))
}
