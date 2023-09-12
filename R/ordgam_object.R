#' Object resulting from the fit of an additive proportional odds model using 'ordgam'
#'
#' An object returned by the \code{\link{ordgam}} function: this is a list
#' with various components related to the fit of such a model.
#'
#' @return An \code{ordgam} object is a list with following elements:
#' \itemize{
#' \item{\code{val} : \verb{ }}{Value of the log-posterior at convergence.}
#' \item{\code{val.start} : \verb{ }}{Value of the log-posterior at the start of the Newton-Raphson (N-R) algorithm.}
#' \item{\code{theta} : \verb{ }}{(Penalized) MLE or MAP of the regression coefficients.}
#' \item{\code{grad} : \verb{ }}{Gradient of the log-posterior at \code{theta}.}
#' \item{\code{Hessian} : \verb{ }}{Hessian of the log-posterior at \code{theta}.}
#' \item{\code{iter} : \verb{ }}{Number of iterations of the N-R algorithm.}
#' \item{\code{llik} : \verb{ }}{Multinomial log likelihood.}
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
#' \item{\code{AIC} : \verb{ }}{Aikake information criterion: AIC = -2 logLik + 2 x edf where edf stands for the effective degrees of freedom.}
#' \item{\code{BIC} : \verb{ }}{Schwarz information criterion: BIC = -2 logLik + n x log(edf) where edf stands for the effective degrees of freedom.}
#' \item{\code{y} : \verb{ }}{Vector containing the values of the ordinal response.}
#' \item{\code{regr} : \verb{ }}{List created by the internal function \code{DesignFormula} and containing diverse objects associated to the model specification, including the part of the design matrix 'X' associated to regressors and its extended version 'Xcal' with B-spline bases for additive term.}
#' \item{\code{ED.Chi2} : \verb{ }}{Matrix containing the Effective Degrees of Freedom associated to the additive terms with their respective significance Chi2 test and P-value.}
#' \item{\code{ED.Tr} : \verb{ }}{Matrix containing the Effective Degrees of Freedom associated to the additive terrms with their respective significance <Tr> test (described by S. Wood, Biometrika 2013) and P-value.}
#' \item{\code{lpost.fun} : \verb{ }}{Function with arguments (theta,lambda,gradient=TRUE,Hessian=TRUE) computing the log-posterior for given regression (and possibly spline) parameters \code{theta} and vector of penalty parameters \code{lambda} associated to the additive terms. Gradient and Hessian are also computed if requested.}
#' \item{\code{lambda0} : \verb{ }}{Initial values for the vector of penalty parameters. Its length corresponds to the number of additive terms.}
#' \item{\code{lambda} : \verb{ }}{(Selected) vector of penalty parameters. Its length corresponds to the number of additive terms.}
#' \item{\code{select.lambda} : \verb{ }}{Logical indicating if \code{lambda} should be selected by maximizing the marginal likelihood or its marginal posterior.}
#' \item{\code{lambda.family} : \verb{ }}{Chosen prior for \code{lambda}: possible choices are "none", "dgamma" (i.e. dgamma(1,1e-4)), "BetaPrime" (BetaPrime(.5,.5)) or "myprior" (with log of the prior density function in \code{myprior}). When "none" is selected, the marginal likelihood is directly maximized.}
#' \item{\code{lprior.lambda} : \verb{ }}{Log of the prior density for the penalty parameters \code{lambda} when \code{select.lambda} is TRUE.}
#' \item{\code{loglambda.loss} : \verb{ }}{The function of \code{log(lambda)} that is minimized to select \code{lambda}. It is minus the log marginal likelihood (when \code{lambda.family} is "none") or minus the log of the marginal posterior for \code{lambda} otherwise.}
#' \item{\code{nu.lpost} : \verb{ }}{Function giving the log of the marginal posterior density of \code{nu=log(lambda)}.}
#' \item{\code{nu.hat} : \verb{ }}{The mode of the marginal posterior density \code{nu.lpost} for \code{nu=log(lambda)}.}
#' \item{\code{V.nu} : \verb{ }}{Variance of the marginal posterior for \code{nu=log(lambda)}.}
#' \item{\code{se.nu} : \verb{ }}{Standard error of \code{nu=log(lambda)}, i.e. the square-root of the diagonal elements of \code{V.nu}.}
#' \item{\code{nu.dp} : \verb{ }}{List containing the parameters of the skew-t approximation to the marginal posterior of \code{nu[j]=loglambda[j]} associated to each of the \code{J} additive terms.}
#' \item{\code{formula} : \verb{ }}{Formula used during the model specification.}
#' \item{\code{elapsed.time} : \verb{ }}{Elapsed time.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling. <doi:10.1177/1471082X231181173>. Preprint: <arXiv:2210.01668>.
#'
#' @seealso \code{\link{ordgam}}, \code{\link{print.ordregr}}, \code{\link{plot.ordgam}}
#'
#' @name ordgam.object
NULL
