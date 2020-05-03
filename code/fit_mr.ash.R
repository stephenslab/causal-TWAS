library(mr.ash.alpha)
#'
#' MR.ASH
#' from mr.ash.alpha github repo
fit.mr.ash = function(X, y, X.test, y.test, seed = 1, sa2 = NULL) {
  
  # run mr.ash
  t.mr.ash           = system.time(
    fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2,
                                max.iter = 2000,
                                standardize = standardize,
                                tol = list(epstol = 1e-12, convtol = 1e-8)))
  beta               = fit.mr.ash$beta
  pip                = 1 - mr.ash.alpha:::get_phi(fit.mr.ash)$phi[,1]
  
  return (list(fit = fit.mr.ash, t = t.mr.ash[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict(fit.mr.ash, X.test), '2')))
}