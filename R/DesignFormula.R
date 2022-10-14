## GOAL
##   Produce a list with the objects required to describe and fit the desired additive model
##    from a list of fixed effect and additive effect covariates
##   --> for use with <LocationScaleFit>
## INPUT:
##   formula: a formula describing the fixed effects and the additive terms in the desired model
##   data: a dataframe containing the data
##   K: number of B-splines to describe an additive term
##   pen.order: desired penalty order for additive terms
##   n: number of units (by default: obtained from the design matrix constructed from the formula and the data frame)
##
## OUTPUT:
##   list with
##     Z: (n x nfixed) design matrix with fixed effects (including a first column of 1)
##     X: (n x J) design matrix with covariates on (0,1) having additive effects
##     nfixed: number of fixed effect regression parameters
##     J: number of additive terms
##     K: number of B-splines in a basis used to estimate an additive term
##     Bx: list with J objects (one per additive term) including (B,Dd,Pd,K,cm)
##     Pd.x, Dd.x: penalty and difference penalty matrices applied on spline parameters of an additive term
##     knots.x: list of length J with the knots associated to each of the J additive terms
##     Bcal: column-stacked matrix with the J centered B-spline bases
##     Xcal: Z column-stacked with the J centered B-spline bases to yield the full design matrix (with column labels)
##     pen.order: reminded
##     lambda.lab: labels for the penalty parameters
##
## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' Internal function extracting design matrices from formulas in the DALSM function and computing penalty related matrices.
#' @keywords internal

DesignFormula = function(formula, data, K=10, pen.order=2, knots.x=NULL, n=NULL, nointercept=FALSE){
  if (!inherits(formula, "formula"))
    stop("Incorrect model formula")
  if ((formula=="~1")&(missing(data))){
    if (is.null(n)){
      cat("Model with only the intercept: the sample size <n> or a data frame should be provided !\n")
      return(NULL)
    }
    XX = model.matrix(~ 1, data = data.frame(rep(1,n)))
  } else {
    ## Extract design matrix
    if (missing(data)) {
      mf <- stats::model.frame(formula)
      XX <- stats::model.matrix(mf)
      if (nointercept){
          if (colnames(XX)[1] == "(Intercept)") XX = XX[,-1]
      }
    }
    else {
      mf <- stats::model.frame(formula, data = data)
      XX <- stats::model.matrix(mf, data = data)
      if (nointercept){
          if (colnames(XX)[1] == "(Intercept)") XX = XX[,-1,drop="FALSE"]
      }
    }
  }
  ## Extract response
  y = stats::model.extract(mf,"response")
  ## Identify additive terms
  smterms <- grepl("s(", colnames(XX), fixed = TRUE)
  X = subset(XX, select=smterms) ## Covariates with additive effects in the design matrix
  ## Reorder the design matrix to start with 'fixed effect' covariates
  idx = c(which(!smterms),which(smterms))
  XX <- subset(XX,select=idx) ## Reordering
  if (any(is.infinite(XX)))
    stop("Covariates contain Inf, NA or NaN values")
  J <- sum(smterms) ## Nbr of additive terms
  nfixed <- ncol(XX) - J
  n <- nrow(XX) ## Nbr of units
  ##
  if (nfixed == 0) Z = NULL
  else Z <- subset(XX, select=1:nfixed) ## 'Fixed effects' part of the design matrix
  ## if (ncol(Z) == 1)
  ##   colnames(Z) <- "(Intercept)"
  ## Some labels
  fixed.lab = colnames(Z) ## Labels of the fixed effects
  additive.lab = lambda.lab = Bcal = NULL
  ## Additives terms
  Bx = NULL ## B-spline for the additive terms
  if (J > 0){
    K = floor(K) ## Make sure that this is an integer
    ## Number of knots
    nknots = K-1 ## Number of knots for the cubic B-splines in the basis
    if (!is.vector(nknots, mode = "numeric") || length(nknots) > 1 || is.na(nknots))
      stop("nknots must be numeric of length 1")
    if (K < 5 || K > 60)
      stop("K must be between 5 and 60")
    ## Penalty order
    pen.order <- floor(pen.order)
    if (pen.order < 1 || pen.order > 4)
      stop("Penalty order must be between 1 and 4")
    ## knots.x = NULL
    set.knots = is.null(knots.x)
    for (j in 1:J){  ## Loop over functional components in location
      xj = XX[,nfixed+j]
      if (set.knots) knots.x[[j]] = seq(min(xj),max(xj),length=nknots)
      ## B-spline matrix for the jth additive term
      Bx[[j]] = centeredBasis.gen(xj,knots=knots.x[[j]],cm=NULL,pen.order)
      colnames(Bx[[j]]$B) = paste(colnames(XX)[nfixed+j],".",1:K,sep="")
    }
    ## Global design matrix
    Bcal = NULL ## Global design matrix for the B-splines of the additive terms
    for (j in 1:J) Bcal = cbind(Bcal,Bx[[j]]$B)
    ## Extra labels
    # additive.lab = colnames(XX)[-(1:nfixed)]
    if (nfixed > 0) additive.lab = unname(sapply(colnames(XX)[-(1:nfixed)], function(x) substring(x,3,nchar(x)-1)))
    else additive.lab = unname(sapply(colnames(XX), function(x) substring(x,3,nchar(x)-1)))
    lambda.lab = paste("lambda.",additive.lab,sep="")
  }
  if (J > 0) {
    Xcal = cbind(Z,Bcal) ## Global design matrix
  } else Xcal = Z
  ##
  Pd.x = Dd.x = NULL
  if (J > 0){
    Pd.x = Bx[[1]]$Pd ## Same penalty matrix for all functional components
    Dd.x = Bx[[1]]$Dd ## Same difference matrix for all functional components
  }
  ##
  ans = list(y=y,Z=Z,X=X,Xcal=Xcal,
             nfixed=nfixed,J=J,K=K)
  if (J > 0){
    ans = c(ans,
            list(
                 Bcal=Bcal,
                 Bx=Bx,Pd.x=Pd.x,Dd.x=Dd.x,knots.x=knots.x,
                 pen.order=pen.order,additive.lab=additive.lab,lambda.lab=lambda.lab))
  }
  return(ans)
}
