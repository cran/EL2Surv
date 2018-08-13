
#' Simulated Survival with Crossing Hazard Functions
#'
#' @description The data frame \code{hazardcross} is simulated from two groups of piecewise exponential
#' lifetime distributions with crossing hazard functions. The estimated survival functions remain ordered even when the estimated hazard 
#' functions are crossed.
#' See \code{\link{supELtest}} for the application.
#' @format The \code{hazardcross} is a data frame with 100 simulated observations of 3 variables,
#' and has the following columns:
#' \itemize{
#'    \item \code{time} the survival time
#'    \item \code{censor} the censoring indicator
#'    \item \code{group} the grouping variable
#' }
#' @seealso \code{\link{supELtest}}

"hazardcross"