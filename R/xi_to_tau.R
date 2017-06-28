
xi_to_tau <- function(tau, fit, fit1, fit2, xi) {
  return (division00(1/(1 + sigma2_hat(tau, fit, fit1, fit2)),
                     1/sigma2_hat(tau, fit, fit1, fit2)) - xi)
}
