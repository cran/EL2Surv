#' The pointwise likelihood ratio test
#' 
#' @description \code{ptwiseELtest} gives pointwise EL statistic values at uncensored time span.
#' The pointwise statistic considers only the decision on each single time point;
#' thus, it is different from the \code{\link[=intELtest]{integral type}} and
#' \code{\link[=supELtest]{sup type}} statistics.
#' @param data a data frame/matrix with 3 columns. The first column is
#' the survival time. The second is the censoring indicator. The last is
#' the grouping variable. An example as the input to \code{data} provided is
#' \code{\link{hepatitis}}.
#' @param g1 the group with the longer survival that
#' should take a value from the third column of \code{data}
#' @param t1 pre-specified \eqn{t_1} based on domain knowledge
#' with the default value of \eqn{0}
#' @param t2 pre-specified \eqn{t_2} based on domain knowledge
#' with the default value of \eqn{\infty}
#' @param sided 2 if two-sided test, and 1 if one-sided test.
#' It assumes the default value of 2.
#' @param nboot number of bootstrap replications in calculating critical values
#' @param alpha pre-specified significance level of the test
#' @param seed the parameter to \code{\link[=Random]{set.seed}} for the random number generator in \R.
#' The \code{set.seed} is used implicitly in \code{supELtest}.
#' @param nlimit the splitting unit. To deal with large data problems, the bootstrap algorithm is
#' to split the number of bootstrap replicates into \code{nsplit} parts. The number \code{nsplit}
#' is the smallest integer not less than \eqn{\left\| U\right\|/}\code{nlimit}.
#' @return \code{ptwiseELtest} returns a list with four elements:
#' \itemize{
#'    \item \code{time_pts} the values of statistics at each uncensored time point
#'    \item \code{decision} logical values. See \code{stat_ptwise}.
#'    \item \code{stat_ptwise} the decision of the test in which the null hypothesis os rejected at a
#'    specific day if the decision exhibits 1 and not rejected if otherwise
#'    \item \code{critval_ptwise} the critical values of the statistic at each uncensored time point
#' }
#' @references H. W. Chang, "Empirical likelihood tests for stochastic
#' ordering based on censored and biased data," \emph{Columbia University Academic Commons} (2014).
#' \url{http://academiccommons.columbia.edu/catalog/ac\%3A177230}
#' @seealso \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{supELtest}}
#' @examples
#' library(EL2Surv)
#' ptwiseELtest(hepatitis)
#' ## It produces the estimates on 47 distinct uncensored days
#' ## out of 57 possibly repeated uncensored days.
#' 
#' ptwiseELtest(hepatitis, t1 = 30, t2 = 60)
#' ## It produces the estimates on 12 distinct uncensored days
#' ## on the restricted time interval [30, 60].
#' 
#' @export
#' @importFrom stats quantile

ptwiseELtest <- function(data, g1 = 1, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, compo = FALSE, seed = 1011, nlimit = 200) {
  at_ts <- neg2ELratio(data, g1, t1, t2, sided, nboot, details.return = TRUE, seed, nlimit)
  if (compo == "TRUE") {
    boot_ptw <- apply(as.matrix(at_ts$Up2_boot_compo), 2, quantile, 1 - alpha)
    return(list(time_pts       = at_ts$bootstrap_ts,
                decision       = as.numeric(at_ts$stat_compo > boot_ptw),
                stat_ptwise    = at_ts$stat_compo,
                critval_ptwise = boot_ptw))
  } else {
    boot_ptw <- apply(as.matrix(at_ts$neg2ELratio_bootstrap_at_ts[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 2, quantile, 1 - alpha)
    return(list(time_pts       = at_ts$ts,
                decision       = as.numeric(at_ts$neg2ELratio_at_ts > boot_ptw),
                stat_ptwise    = at_ts$neg2ELratio_at_ts,
                critval_ptwise = boot_ptw))
  }#END compo if-else
}#END ptwiseELtest