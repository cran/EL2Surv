#' The integrated likelihood ratio test
#' 
#' @description \code{intELtest} gives a class of the weighted likelihood ratio statistics:
#' \deqn{\sum_{t\in U}w(t)\{-2\log R(t)\},}
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
#' @param wt a string for the integral statistic with a specific weight function.
#' There are four types of integral statistics provided: \code{"p.event"}, \code{"dt"},
#' \code{"db"}, and \code{"dF"}. See 'Details' for more about the integral statistics.
#' @param alpha pre-specified significance level of the test
#' @param seed the parameter to \code{\link[=Random]{set.seed}} for the random number generator in \R.
#' The \code{set.seed} is used implicitly in \code{intELtest}.
#' @param nlimit the splitting unit. To deal with large data problems, the bootstrap algorithm is
#' to split the number of bootstrap replicates into \code{nsplit} parts. The number \code{nsplit}
#' is the smallest integer not less than \eqn{\left\| U\right\|/}\code{nlimit}.
#' @return \code{intELtest} returns a list with three elements:
#' \itemize{
#'    \item \code{critval} the critical value
#'    \item \code{teststat} the resulting integrated test statistic
#'    \item \code{pvalue} the p-value based on the integrated statistic
#' }
#' @details \code{intELtest} calculates the weighted likelihood ratio statistics:
#' \deqn{\sum_{i=1}^{h}w_i\cdot \{-2\log R(t_i)\},}
#' where \eqn{w_1,...,w_h} are the values of the weight function evaluated at
#' the distinct ordered uncensored times \eqn{t_1,...,t_h} in \eqn{U}.
#' There are four types of weight functions considered.
#' \itemize{
#'    \item (\code{wt = "dt"}) \cr
#'      By means of an extension of the integral statistic derived by Pepe and Fleming (1989),
#'      \deqn{w_i=\left\{\begin{array}{ll}t_{i+1}-t_i & \textrm{if }i\neq h\\ 0 & \textrm{if }i=h\end{array}\right.}
#'    \item (\code{wt = "p.event"}) \cr
#'      According to the integral statistic derived by Uno et al. (2013):
#'      \deqn{w_i=\frac{1}{n_1+n_2},}
#'      where \eqn{n_1} and \eqn{n_2} are the sample sizes of each group.
#'      The role of \eqn{w_i} assigns equal weight \eqn{1/n} to each observation when no
#'      tie is involved (\eqn{n = n_1 + n_2}); otherwise, it assigns heavy weight to
#'      the observations with multiplicity.
#'    \item (\code{wt = "dF"}) \cr
#'      Based on the integral statistic built by Barmi and McKeague (2013),
#'      the weight function is the derivative of the empirical distribution function \eqn{\hat{F}(t)}. 
#'      \deqn{w_i=\left.{\frac{d\hat{F}(t)}{dt}}\right|_{t=t_i}}
#'      This is an empirical version of taking expectation.
#'    \item (\code{wt = "db"}) \cr
#'      The weight function is of the form:
#'      \deqn{w_i=\left.\frac{d\hat{b}(t)}{dt}\right|_{t=t_i}}
#'      with \eqn{\hat{b}(t)=\hat{\sigma}^2(t)/(1+\hat{\sigma}^2(t))}.
#'      The \eqn{\hat{b}(t)} is chosen so that the limiting distribution
#'      \deqn{\int^{x_2}_{x_1}\frac{B^2_+(x)}{x(1 - x)}dx}
#'      is the same as the asymptotic null distribution in Barmi and McKeague (2013) with the
#'      \eqn{[x_1, x_2]} restriction.
#' }
#' @references
#' \itemize{
#'    \item H. W. Chang, "Empirical likelihood tests for stochastic
#'      ordering based on censored and biased data," \emph{Columbia University Academic Commons} (2014).
#'      \url{http://academiccommons.columbia.edu/catalog/ac\%3A177230}
#'    \item M. S. Pepe and T. R. Fleming, "Weighted Kaplan-Meier
#'      Statistics: A Class of Distance Tests for Censored Survival Data," \emph{Biometrics},
#'      Vol. 45, No. 2, pp. 497-507 (1989).
#'      \url{https://www.jstor.org/stable/2531492?seq=1#page_scan_tab_contents}
#'    \item H. Uno, L. Tian, B. Claggett, and L. J. Wei, "A versatile test for equality of
#'      two survival functions based on weighted differences of Kaplan-Meier curves,"
#'      \emph{Statistics in Medicine}, Vol. 34, No. 28, pp. 3680-3695 (2015).
#'      \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.6591/abstract}
#'    \item H. E. Barmi and I. W. McKeague, "Empirical likelihood-based tests
#'      for stochastic ordering," \emph{Bernoulli}, Vol. 19, No. 1, pp. 295-307 (2013).
#'      \url{https://projecteuclid.org/euclid.bj/1358531751}
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{supELtest}}, \code{\link{ptwiseELtest}}
#' @examples
#' library(EL2Surv)
#' intELtest(hepatitis, 1, sided = 2, wt = "p.event")
#' 
#' ## OUTPUT:
#' ## $critval
#' ## [1] 0.8993514
#' ## 
#' ## $teststat
#' ## [1] 1.406029
#' ## 
#' ## $pvalue
#' ## [1] 0.012
#' 
#' @export
#' @importFrom stats quantile

intELtest <- function(data, g1 = 1, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, wt = "p.event", alpha = 0.05, compo = FALSE, seed = 1011, nlimit = 200) {
  at_ts <- neg2ELratio(data, g1, t1, t2, sided, nboot, details.return = TRUE, seed, nlimit)
  if (is.null(at_ts)) return (NULL)
  if (compo == TRUE) {
    if (wt == "p.event") {
      inttest_dbarNt                  <- sum(at_ts$stat_compo * at_ts$wt_dbarNt)
      Up2_boot_H1_1sided_times_dbarNt <- at_ts$Up2_boot_compo * at_ts$wt_dbarNt_boot
      int_dbarNt_boot_H1              <- apply(Up2_boot_H1_1sided_times_dbarNt, 1, sum)
      int_dbarNtEL_SOcrit             <- as.vector(quantile(int_dbarNt_boot_H1, 1 - alpha))
      
      critval  <- int_dbarNtEL_SOcrit
      teststat <- inttest_dbarNt
      pvalue   <- mean(int_dbarNt_boot_H1 > inttest_dbarNt)
    } else if (wt == "dt") {
      inttest_dt                  <- sum((at_ts$stat_compo_dt * at_ts$wt_dt)[-length(at_ts$wt_dt)])
      Up2_boot_H1_1sided_times_dt <- as.matrix((at_ts$Up2_boot_compo_dt * at_ts$wt_dt_boot)[, -length(at_ts$wt_dt)])
      int_dt_boot_H1              <- apply(Up2_boot_H1_1sided_times_dt, 1, sum)
      int_dtEL_SOcrit             <- as.vector(quantile(int_dt_boot_H1, 1 - alpha))
      
      critval  <- int_dtEL_SOcrit
      teststat <- inttest_dt
      pvalue   <- mean(int_dt_boot_H1 > inttest_dt)
    } else if (wt == "db") {
      inttest_db                     <- sum(at_ts$stat_compo * at_ts$wt_db)
      Up2_boot_H1_1sided_times_ds1ps <- at_ts$Up2_boot_compo * at_ts$wt_db_boot
      int_boot_H1                    <- apply(Up2_boot_H1_1sided_times_ds1ps, 1, sum)
      intEL_SOcrit                   <- as.vector(quantile(int_boot_H1, 1 - alpha))
      
      critval  <- intEL_SOcrit
      teststat <- inttest_db
      pvalue   <- mean(int_boot_H1 > inttest_db)
    } else if (wt == "dF") {
      inttest_dF                  <- sum(at_ts$stat_compo * at_ts$wt_dF)
      Up2_boot_H1_1sided_times_dF <- at_ts$Up2_boot_compo * at_ts$wt_dF_boot
      int_dF_boot_H1              <- apply(Up2_boot_H1_1sided_times_dF, 1, sum)
      int_dFEL_SOcrit             <- as.vector(quantile(int_dF_boot_H1, 1 - alpha))
      
      critval  <- int_dFEL_SOcrit
      teststat <- inttest_dF
      pvalue   <- mean(int_dF_boot_H1 > inttest_dF)
    }
  } else {
    if (wt == "p.event") {
      inttest_dbarNt                  <- sum(at_ts$neg2ELratio_at_ts * at_ts$wt_dbarNt[at_ts$lowerbindx_boot:at_ts$upperbindx_boot])
      Up2_boot_H1_1sided_times_dbarNt <- at_ts$neg2ELratio_bootstrap_at_ts * at_ts$wt_dbarNt_boot
      int_dbarNt_boot_H1              <- apply(as.matrix(Up2_boot_H1_1sided_times_dbarNt[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 1, sum)
      int_dbarNtEL_SOcrit             <- as.vector(quantile(int_dbarNt_boot_H1, 1 - alpha))
      
      critval  <- int_dbarNtEL_SOcrit
      teststat <- inttest_dbarNt
      pvalue   <- mean(int_dbarNt_boot_H1 > inttest_dbarNt)
    } else if (wt == "dt") {
      inttest_dt                  <- sum(at_ts$neg2ELratio_at_ts[-(at_ts$upperbindx_boot - at_ts$lowerbindx_boot + 1)] * at_ts$wt_dt[at_ts$lowerbindx_boot:(at_ts$upperbindx_boot - 1)])
      Up2_boot_H1_1sided_times_dt <- at_ts$neg2ELratio_bootstrap_at_ts * at_ts$wt_dt_boot
      int_dt_boot_H1              <- apply(as.matrix(Up2_boot_H1_1sided_times_dt[, at_ts$lowerbindx_boot:(at_ts$upperbindx_boot - 1)]), 1, sum)
      int_dtEL_SOcrit             <- as.vector(quantile(int_dt_boot_H1, 1 - alpha))
      
      critval  <- int_dtEL_SOcrit
      teststat <- inttest_dt
      pvalue   <- mean(int_dt_boot_H1 > inttest_dt)
    } else if (wt == "db") {
      inttest_db                     <- sum(at_ts$neg2ELratio_at_ts * at_ts$wt_db[at_ts$lowerbindx_boot:at_ts$upperbindx_boot])
      Up2_boot_H1_1sided_times_ds1ps <- at_ts$neg2ELratio_bootstrap_at_ts * at_ts$wt_db_boot
      int_boot_H1                    <- apply(as.matrix(Up2_boot_H1_1sided_times_ds1ps[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 1, sum)
      intEL_SOcrit                   <- as.vector(quantile(int_boot_H1, 1 - alpha))
      
      critval  <- intEL_SOcrit
      teststat <- inttest_db
      pvalue   <- mean(int_boot_H1 > inttest_db)
    } else if (wt == "dF") {
      inttest_dF                  <- sum(at_ts$neg2ELratio_at_ts * at_ts$wt_dF[at_ts$lowerbindx_boot:at_ts$upperbindx_boot])
      Up2_boot_H1_1sided_times_dF <- at_ts$neg2ELratio_bootstrap_at_ts * at_ts$wt_dF_boot
      int_dF_boot_H1              <- apply(as.matrix(Up2_boot_H1_1sided_times_dF[, at_ts$lowerbindx_boot:at_ts$upperbindx_boot]), 1, sum)
      int_dFEL_SOcrit             <- as.vector(quantile(int_dF_boot_H1, 1 - alpha))
      
      critval  <- int_dFEL_SOcrit
      teststat <- inttest_dF
      pvalue   <- mean(int_dF_boot_H1 > inttest_dF)
    }
  }
  return (list(critval = critval, teststat = teststat, pvalue = pvalue))
}
