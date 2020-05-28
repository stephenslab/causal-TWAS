#' add up two mr.ash models.
#' y = X1 + X2 + e
#' @param iter iterations

mr.ash2 <- function(X1, X2, y, iter = 20){
  r1 <- y
  fit1 <- mr.ash(X1, r1, method = "caisa")
  for (i in 1:iter) {
    r2 <- y - X1 %*% fit1$beta
    fit2 <- mr.ash(X2, r2, method = "caisa")

    r1 <- y - X2 %*% fit2$beta
    fit1 <- mr.ash(X1, r1, method = "caisa")
    gc()
  }
  fit <- list("fit1" = fit1, "fit2" = fit2)
  return(fit)
}

#' add up two mr.ash models.
#' update of sigma, init with mr.ash, one iter each for X1 and X2.
mr.ash2s <- function(X1, X2, y, iter = 100){
  n <- dim(X1)[1]
  p1 <- dim(X1)[2]
  p2 <- dim(X2)[2]

  # initiate using mr.ash for X1, X2 combined
  X <- cbind(X1, X2)
  fit.init <- mr.ash(X, y)

  beta1 <- fit.init$beta[1: p1, , drop = F]
  beta2 <- fit.init$beta[(p1+1): (p1 + p2), , drop = F]

  pi1 <- fit.init$pi
  pi2 <- fit.init$pi
  sa2.1 <- fit.init$sa2
  sa2.2 <- fit.init$sa2

  r <- y - X1 %*% beta1 - X2 %*% beta2
  sigma2 <- fit.init$sigma2

  K <- dim(fit.init$pi)[1]

  rm(fit.init)
  rm(X);gc()

  # start iteration (from X1)
  piK1_rec <- rep(0, iter)
  piK2_rec <- rep(0, iter)
  sigma2_rec <- rep(0, iter)

  for (i in 1:iter) {
    y1 <- y - X2 %*% beta2
    fit1 <- mr.ash(X1, y1,
                   beta.init = beta1,
                   sigma2 = sigma2,
                   sa2 = sa2.1,
                   pi = pi1,
                   max.iter = 1)

    y2 <- y - X1 %*% beta1

    fit2 <- mr.ash(X2, y2,
                   beta.init = beta2,
                   sigma2 = sigma2,
                   sa2 = sa2.2,
                   pi = pi2,
                   max.iter = 1)

    # update beta
    beta1 <- fit1$beta
    beta2 <- fit2$beta

    # update pi
    pi1 <- fit1$pi
    pi2 <- fit2$pi
    sa2.1 <- fit1$sa2
    sa2.2 <- fit2$sa2

    # update sigma
    r1 <- y1 - X1 %*% beta1
    s1 <- n * fit1$sigma2 - sum(r1^2) # s1 is d^2 (btilde -bhat)bhat for X1 in algorithm 4 in mr.ash paper.

    r2 <- y2 - X2 %*% beta2
    s2 <- n * fit2$sigma2 - sum(r2^2) # s2 is d^2 (btilde -bhat)bhat for X2 in algorithm 4 in mr.ash paper.

    sigma2 <- 1/n * (sum(r^2) + s1 + s2)

    # update r
    r <- y - X1 %*% beta1 - X2 %*% beta2

    piK1_rec[i] <- pi1[K,]
    piK2_rec[i] <- pi2[K,]
    sigma2_rec[i] <- sigma2

    gc()
  }
  fit <- list("fit1" = fit1, "fit2" = fit2, "piK1" = piK1_rec, "piK2" = piK2_rec, "sigma2" = sigma2_rec)
  return(fit)
}
