#' The maximally selected likelihood ratio test
#' 
#' @description \code{supELtest} provides a maximal deviation type statistics that
#' is better adapted at detecting local differences:
#' \deqn{\sup_{t\in U}\{-2\log R(t)\},}
#' where \eqn{w(t)} is an objective weight function, and \eqn{R(t)} is an empirical likelihood
#' (EL) ratio that compares two survival functions at each time point \eqn{t} in the set of
#' observed uncensored lifetimes, \eqn{U}.
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
#' @return \code{supELtest} returns a list with three elements:
#' \itemize{
#'    \item \code{critval} the critical value
#'    \item \code{teststat} the resulting integrated test statistic
#'    \item \code{pvalue} the p-value based on the integrated statistic
#' }
#' @references H. W. Chang, "Empirical likelihood tests for stochastic
#' ordering based on censored and biased data," \emph{Columbia University Academic Commons} (2014).
#' \url{http://academiccommons.columbia.edu/catalog/ac\%3A177230}
#' @seealso \code{\link{hepatitis}}, \code{\link{intELtest}}, \code{\link{ptwiseELtest}}
#' @examples
#' library(EL2Surv)
#' supELtest(hepatitis, 1, sided = 2)
#' 
#' ## OUTPUT:
#' ## $critval
#' ## [1] 7.247709
#' ## 
#' ## $teststat
#' ## [1] 10.3581
#' ## 
#' ## $pvalue
#' ## [1] 0.013
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
  return(list(critval = EL_SOcrit, teststat = suptest, pvalue = mean(sup_boot_H1 > suptest)))
}
