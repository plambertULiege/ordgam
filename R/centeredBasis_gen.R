## Inspired by Wood (2017, pp. 175-176)
##
## GOAL:
##  Produce a B-spline basis matrix with columns recentered to handle the
##   identifiability constraint in additive models
##
## INPUT:
##  x: vector of values where to compute the "centered" B-spline basis
##  knots: vector of knots (that should cover the values in x)
##
## OUTPUT:
##  list with
##    B: centered B-spline matrix (with columns recentered to have mean 0
##         over equi-spaced x values on (0,1)) with the last column removed
##         omitted from the original B-spline matrix
##    Dd: difference matrix for the associated centered B-spline matrix
##    Dd: penalty matrix for the associated centered B-spline matrix
##    K: number of "free" spline parameters
##    cm: mean values substracted from each column of the original B-spline matrix
##
## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
centeredBasis.gen = function(x,knots,cm=NULL,pen.order=2){
  if ((max(x)>max(knots))|(min(x)<min(knots))){
    cat("The knots do no cover the values of x !!\n")
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
    ##  to compute a recentering value per column of the matrix
    B.ref = cubicBsplines::Bsplines(seq(min(knots),max(knots),length=100), knots.x)[,idx]
    cm = colMeans(B.ref)
  }
  B = sweep(B,2,cm) ## Substract column mean from each column
  Dd = diff(diag(ncol(temp)),dif=pen.order.x)[,idx] ## Penalty order <---- Important to remove last column !!
  Pd = t(Dd) %*% Dd
  return(list(B=B,Dd=Dd,Pd=Pd,K=ncol(B),cm=cm))
}
