#' The maximally selected likelihood ratio test
#' 
#' @description \code{supELtest} provides a maximal deviation type statistics that
#' is better adapted at detecting local differences:
#' \deqn{\sup_{t\in U}\{-2\log R(t)\},}
#' where \eqn{R(t)} is an empirical likelihood
#' (EL) ratio that compares two survival functions at each time point \eqn{t} in the set of
#' observed uncensored lifetimes, \eqn{U}.
#' @param data a data frame/matrix with 3 columns. The first column is
#' the survival time. The second is the censoring indicator. The last is
#' the grouping variable. An example as the input to \code{data} provided is
#' \code{\link{hepatitis}}.
#' @param g1 the group with longer survival in one-sided testing with the default value of \eqn{1}.
#' @param t1 pre-specified \eqn{t_1} based on domain knowledge
#' with the default value of \eqn{0}
#' @param t2 pre-specified \eqn{t_2} based on domain knowledge
#' with the default value of \eqn{\infty}
#' @param sided 2 if two-sided test, and 1 if one-sided test.
#' It assumes the default value of \eqn{2}.
#' @param nboot number of bootstrap replications in calculating critical values
#' with the defualt value of \eqn{1000}.
#' @param alpha pre-specified significance level of the test with the default value of \eqn{0.05}
#' @param compo FALSE if taking the standardized square of the difference as the local statisic
#' for two-sided testing, and TRUE if constructing for one-sided testing, but only the positive
#' part of the difference included. It assumes the default value of \eqn{FALSE}.
#' @param seed the parameter with the default value of \eqn{1011} to \code{\link[=Random]{set.seed}} for 
#' generating bootstrap-based critical values in \R.
#' The \code{set.seed} is used implicitly in \code{intELtest}.
#' @param nlimit the splitting unit with the default value of \eqn{200}. To deal with large data problems, the bootstrap algorithm is
#' to split the number of bootstrap replicates into \code{nsplit} parts. The number \code{nsplit}
#' is the smallest integer not less than \eqn{\left\| U\right\|/}\code{nlimit}.
#' @return \code{supELtest} returns a list with three elements:
#' \itemize{
#'    \item \code{teststat} the resulting integrated test statistic
#'    \item \code{critval} the critical value
#'    \item \code{pvalue} the p-value based on the integrated statistic
#' }
#' @references H.-w. Chang and I. W. McKeague, "Empirical likelihood based tests 
#' for stochastic ordering under right censorship," \emph{Electronic Journal of Statistics},
#' Vol. 10, No. 2, pp. 2511-2536 (2016).
#' @seealso \code{\link{hazardcross}}, \code{\link{intELtest}}, \code{\link{ptwiseELtest}}
#' @examples
#' library(EL2Surv)
#' supELtest(hazardcross)
#' 
#' ## OUTPUT:
#' ## $teststat
#' ## [1] 8.945539
#' ## 
#' ## $critval
#' ## [1] 8.738189
#' ## 
#' ## $pvalue
#' ## [1] 0.045
#' 
#' @export
#' @importFrom stats quantile

supELtest <- function(data, g1 = 1, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, compo = FALSE, seed = 1011, nlimit = 200) {
  at_ts <- neg2ELratio(data, g1, t1, t2, sided, nboot, details.return = TRUE, seed, nlimit)
  if (is.null(at_ts)) return (NULL)
  if (compo == TRUE){
    suptest     <- max(at_ts$stat_compo)
    sup_boot_H1 <- apply(as.matrix(at_ts$Up2_boot_compo), 1, max)
    EL_SOcrit   <- as.vector(quantile(sup_boot_H1, 1 - alpha))
  } else {
    suptest     <- max(at_ts$neg2ELratio_at_ts)
    sup_boot_H1 <- apply(as.matrix(at_ts$neg2ELratio_bootstrap_at_ts[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 1, max)
    EL_SOcrit   <- as.vector(quantile(sup_boot_H1, 1 - alpha))
  }
  return(list(teststat = suptest, critval = EL_SOcrit, pvalue = mean(sup_boot_H1 > suptest)))
}
