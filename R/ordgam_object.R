#' Object resulting from the fit of an additive proportional odds model using 'ordgam'
#'
#' An object returned by the \code{\link{ordgam}} function: this is a list
#' with various components related to the fit of such a model.
#'
#' @return An \code{ordgam} object is a list with  following elements:
#' \itemize{
#' \item{\code{val} : \verb{ }}{Value of the log-posterior at convergence.}
#' \item{\code{val.start} : \verb{ }}{Value of the log-posterior at the start of the Newton-Raphson (N-R) algorithm.}
#' \item{\code{theta} : \verb{ }}{(Penalized) MLE or MAP of the regression coefficients.}
#' \item{\code{grad} : \verb{ }}{Gradient of the log-posterior at \code{theta}.}
#' \item{\code{Hessian} : \verb{ }}{Hessian of the log-posterior at \code{theta}.}
#' \item{\code{iter} : \verb{ }}{Number of iterations of the N-R algorithm.}
#' \item{\code{Hessian0} : \verb{ }}{Hessian of the (non-penalized) log-likelihood at \code{theta}.}
#' \item{\code{Sigma.theta} : \verb{ }}{Variance-covariance of 'theta'.}
#' \item{\code{ED.full} : \verb{ }}{Effective degrees of freedom associated to each regression parameter, penalized parameters included.}
#' \item{\code{se.theta} : \verb{ }}{Standard errors of the regression coefficents.}
#' \item{\code{theta.mat} : \verb{ }}{Matrix containing the point estimate, standard error, credible interval, Z-score and P-value for \code{theta}.}
#' \item{\code{nc} : \verb{ }}{Number of categories for the ordinal response.}
#' \item{\code{nalpha} : \verb{ }}{Number of intercepts in the proportional odds model (=\code{nc}-1) .}
#' \item{\code{nbeta} : \verb{ }}{Number of regression parameters (intercepts excluded).}
#' \item{\code{nfixed} : \verb{ }}{Number of non-penalized regression parameters.}
#' \item{\code{ci.level} : \verb{ }}{Nominal coverage of the credible intervals (Default: .95).}
#' \item{\code{n} : \verb{ }}{Sample size.}
#' \item{\code{call} : \verb{ }}{Function call.}
#' \item{\code{descending} : \verb{ }}{Logical indicating if the odds of the response taking a value in the upper scale should be preferred over values in the lower scale.}
#' \item{\code{use.prior} : \verb{ }}{Logical indicating if a prior (such as a penalty) is assumed for the regression parameters.}
#' \item{\code{lpost} : \verb{ }}{Value of the log-posterior at convergence.}
#' \item{\code{levidence} : \verb{ }}{Log of the marginal likelihood (also named 'evidence').}
#' \item{\code{y} : \verb{ }}{Vector containing the values of the ordinal response.}
#' \item{\code{regr} : \verb{ }}{List created by the internal function \code{DesignFormula} and containing diverse objects associated to the model specification, including the part of the design matrix 'X' associated to regressors and its extended version 'Xcal' with B-spline bases for additive term.}
#' \item{\code{ED} : \verb{ }}{Matrix containing the Effective Degrees of Freedom associated to the additive terrms with their respective significance test and P-value.}
#' \item{\code{Wood.test} : \verb{ }}{Logical indicating if the significance test described by S. Wood (Biometrika 2013)  should be used to test for additive term significance. When FALSE, an approximate Chi2 test based on the effective degrees of freedom for each additive term is reported.}
#' \item{\code{lpost.fun} : \verb{ }}{Function with arguments (theta,lambda,gradient=TRUE,Hessian=TRUE) computing the log-posterior for given regression (and possibly spline) parameters \code{theta} and vector of penalty parameters \code{lambda} associated to the additive terms. Gradient and Hessian are also computed if requested.}
#' \item{\code{lambda} : \verb{ }}{Vector of penalty parameters. Its length corresponds to the number of additive terms.}
#' \item{\code{formula} : \verb{ }}{Formula used during the model specification.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @seealso \code{\link{ordgam}}, \code{\link{print.ordregr}}, \code{\link{plot.ordgam}}
#'
#' @name ordgam.object
NULL
