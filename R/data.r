#' Example data set: Game addiction
#'
#' @description
#' Dataset containing variables on the gaming addiction of young people.
#' The data consists of gaming addiction, coded as either 'no' or 'yes',
#' measured over four consecutive years for 500 individuals, three covariates
#' and one time-dependent covariate. The yearly states
#' are recorded in columns 1 (\code{T1_abuse}) to 4 (\code{T4_abuse}).
#'
#' The three covariates are
#'
#' \itemize{
#'   \item \code{Gender} (female or male),
#'   \item \code{Age} (measured at time 1),
#'   \item \code{Track} (school or apprenticeship).
#' }
#'
#' The time-varying covariate consists of the individual's relationship to
#' gambling at each of the four time points, appearing in columns
#' \code{T1_gambling}, \code{T2_gambling},
#' \code{T3_gambling}, and \code{T4_gambling}. The states are either
#' no, gambler or problematic gambler
#'
#'
#' @docType data
#' @keywords datasets
#' @name gameadd
#' @usage data(gameadd)
#' @format A data frame containing 500 rows, 4 states variable, 3 covariates
#' and a time-dependent covariate.
"gameadd"
