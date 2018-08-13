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
#' @param g1 the group with longer survival in one-sided testing with the default value of \eqn{1}.
#' @param t1 pre-specified \eqn{t_1} based on domain knowledge
#' with the default value of \eqn{0}
#' @param t2 pre-specified \eqn{t_2} based on domain knowledge
#' with the default value of \eqn{\infty}
#' @param sided 2 if two-sided test, and 1 if one-sided test.
#' It assumes the default value of \eqn{2}.
#' @param nboot number of bootstrap replications in calculating critical values
#' with the defualt value of \eqn{1000}.
#' @param wt a string for the integral statistic with a specific weight function.
#' There are four types of integral statistics provided: \code{"p.event"}, \code{"dF"},
#' \code{"dt"}, and \code{"db"}. It assumes the default value of \code{"p.event"}. See 'Details' for more about the integral statistics.
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
#' @return \code{intELtest} returns a list with three elements:
#' \itemize{
#'    \item \code{teststat} the resulting integrated test statistic
#'    \item \code{critval} the critical value
#'    \item \code{pvalue} the p-value based on the integrated statistic
#' }
#' @details \code{intELtest} calculates the weighted likelihood ratio statistics:
#' \deqn{\sum_{i=1}^{h}w_i\cdot \{-2\log R(t_i)\},}
#' where \eqn{w_1,...,w_h} are the values of the weight function evaluated at
#' the distinct ordered uncensored times \eqn{t_1,...,t_h} in \eqn{U}.
#' There are four types of weight functions considered.
#' \itemize{
#'     \item (\code{wt = "p.event"}) \cr
#'      This default option is an objective weight,
#'      \deqn{w_i=\frac{d_i}{n}}
#'      In other words, this \eqn{w_i} assigns weight proportional to the number of events
#'      at each observed uncensored time \eqn{t_i}.
#'    \item (\code{wt = "dF"}) \cr
#'      Based on the integral statistic built by Barmi and McKeague (2013), another weigth function is 
#'      \deqn{w_i= \hat{F}(t_i)-\hat{F}(t_{i-1})}
#'      for \eqn{i=1,\ldots,m},where \eqn{\hat{F}(t)=1-\hat{S}(t)}, \eqn{\hat{S}(t)} is the pooled KM estimator, and \eqn{t_0 \equiv 0}. 
#'      This reduces to the objective weight when there is no censoring.The resulting \eqn{I_n} can be seen as an empirical 
#'      version of \eqn{E(-2\log\mathcal{R}(T))}, where \eqn{T} denotes the lifetime random variable of interest distributed 
#'      as the common distribution under \eqn{H_0}.
#'    \item (\code{wt = "dt"}) \cr
#'      By means of an extension of the integral statistic derived by Pepe and Fleming (1989), another weight function is 
#'      \deqn{w_i= t_{i+1}-t_i} 
#'      for \eqn{i=1,\ldots,m}, where \eqn{t_{m+1} \equiv t_{m}}. This gives more weight to the time intervals where 
#'      there are fewer observed uncensored times, but may be affected by extreme observations. 
#'    \item (\code{wt = "db"}) \cr
#'      According to a weigthing method mentioned in Chang and McKeague (2016), the other weight function is
#'      \deqn{w_i= \hat{b}(t_i)-\hat{b}(t_{i-1})}
#'      where \eqn{\hat{b}(t)=\hat{\sigma}^2(t)/(1+\hat{\sigma}^2(t))}, and \eqn{\hat{\sigma}^2(t)} is given.
#'      The \eqn{\hat{b}(t)} is chosen so that the limiting distribution is the same as the asymptotic null 
#'      distribution in EL Barmi and McKeague (2013).
#' }
#' @references
#' \itemize{
#'    \item H.-w. Chang and I. W. McKeague, "Empirical likelihood based tests 
#'      for stochastic ordering under right censorship," \emph{Electronic Journal of Statistics},
#'      Vol. 10, No. 2, pp. 2511-2536 (2016).
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
#' intELtest(hepatitis)
#' 
#' ## OUTPUT:
#' ## $teststat
#' ## [1] 1.406016
#' ## 
#' ## $critval
#' ## [1] 0.8993514
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
  return (list(teststat = teststat, critval = critval, pvalue = pvalue))
}
