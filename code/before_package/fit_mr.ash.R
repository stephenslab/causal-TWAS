library(mr.ash.alpha)
#'
#' MR.ASH
#' from mr.ash.alpha github repo
# fit.mr.ash = function(X, y, X.test, y.test, seed = 1, sa2 = NULL) {
#
#   # run mr.ash
#   t.mr.ash           = system.time(
#     fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2,
#                                 max.iter = 2000,
#                                 standardize = standardize,
#                                 tol = list(epstol = 1e-12, convtol = 1e-8)))
#   beta               = fit.mr.ash$beta
#   pip                = 1 - mr.ash.alpha:::get_phi(fit.mr.ash)$phi[,1]
#
#   return (list(fit = fit.mr.ash, t = t.mr.ash[3], beta = beta, pip = pip,
#                rsse = norm(y.test - predict(fit.mr.ash, X.test), '2')))
# }

get_phi <- function(fit) {
  X <- fit$data$X

  # compute residual
  r            = fit$data$y - X %*% fit$beta

  if (class(X) == "matrix") {
    Xt <- t(X)
  } else if (class(X) == "FBM") {
    if (file.exists(paste0(drop_ext(X$backingfile), ".t.rds"))){
      Xt <- readRDS(paste0(drop_ext(X$backingfile), ".t.rds"))
    } else{
      Xt <- big_transpose(X, backingfile = paste0(drop_ext(X$backingfile), ".t"))
      Xt$save()
    }
  }
  # compute bw and S2inv

  bw           = as.vector((Xt %*% r) + fit$data$w * fit$beta)
  S2inv        = 1 / outer(fit$data$w, 1/fit$data$sa2, '+');

  # compute mu, phi
  mu           = bw * S2inv;
  phi          = -log(1 + outer(fit$data$w, fit$data$sa2))/2 + mu * (bw / 2 /
                                                                       fit$sigma2);
  phi          = c(fit$pi) * t(exp(phi - apply(phi,1,max)));
  phi          = t(phi) / colSums(phi);
  return (list(phi = phi, mu = mu, r = r))
}

get_pve <- function(fit) {
  y <- data.frame(y = fit$data$y)
  n <- dim(y)[1]
  v <- fit$data$sa2 %*% fit$pi
  v1 <- sum(fit$data$w)/n * v
  pve <- v1/(v1 + fit$sigma2)
  return (pve)
}

get_pip <- function(fit) {
  pip  <- 1 - get_phi(fit)$phi[,1]
  return (pip)
}

