#' Perception of gay men and lesbians in Belgium.
#'
#' @docType data
#'
#' @description Data extracted from the European Social Survey (2018) conducted in
#' a series of European countries including Belgium.
#' Each of the 1737 participants (aged at least 15) was asked to react to the following statement,
#' "Gay men and lesbians should be free to live their own life as they wish",
#' with a positioning on a Likert scale going from 1 (=Agree strongly) to 5 (=Disagree strongly),
#' with 3 labelled as Neither agree nor disagree.
#'#'
#' @usage data(freehmsDataBE)
#'
#' @format A data frame with 1737 rows and 5 columns.
#' \describe{
#'  \item{\code{freehms}}{Ordinal response (1: Agree strongly to 5: Disagree strongly).}
#'  \item{\code{gndr}}{Gender of the respondent.}
#'  \item{\code{age}}{Age of the respondent.}
#'  \item{\code{eduyrs}}{Number of years of education completed.}
#'  \item{\code{region}}{Belgian region of residence: WAL (Wallonia), FL (Flanders) or BXL (Brussels).}
#' }
#'
#' @references
#' Lambert, P. and Gressani, 0. (2023) Penalty parameter selection and asymmetry corrections
#' to Laplace approximations in Bayesian P-splines models. Statistical Modelling (in press). Preprint: arXiv:2210.01668.
#'
#' @references European Social Survey Round 9 Data (2018)
#'  Data file edition 3.1. NSD - Norwegian Centre for Research Data, Norway.
#'
"freehmsDataBE"
