library(mr.ash.alpha)
library(glmnet)

#' this is a simplified version of mr.ash
#' Difference:
#' - doesn't have the remove covariates function. no `Z` will be accepted.
#' - preprocessing (center/scale) for X and y need to be done before running
#'   this function (no `standardize`, `intercept` arguments).
#' - precalculation for w can be provided as an argument `w`.

mr.ashs                      = function(X, y, sa2 = NULL,
                                       K = 20, method = c("caisa","accelerate","block","sigma_scaled","sigma_indep"),
                                       w = NULL,
                                       max.iter = 1000, min.iter = 1,
                                       beta.init = NULL,
                                       update.pi = TRUE, pi = NULL,
                                       update.sigma = TRUE, sigma2 = NULL,
                                       update.order = NULL,
                                       mixsqpiter = 5, mode = FALSE,
                                       tol = set_default_tolerance(),
                                       verbose = TRUE){

  # get sizes
  n            = dim(X)[1];
  p            = dim(X)[2];

  # match method
  method      <- match.arg(method)

  # set default tolerances unless specified
  tol0         = set_default_tolerance()
  tol          = modifyList(tol0,tol,keep.null = TRUE)

  data         = list(y = y, alpha = mean(y))

  # initialize beta
  if ( is.null(beta.init) ){
    data$beta  = as.vector(double(p));
  } else {
    data$beta  = as.vector(beta.init);
  }

  data$beta[1] = data$beta[1] + 0;      # to make sure beta.init is not modified

  # initialize r
  r            = data$y - X %*% data$beta;

  # sigma2
  if ( is.null(sigma2) ) {
    sigma2 = c(var(r));
  }

  # set sa2 if missing
  if ( is.null(sa2) ) {
    sa2             = (2^((0:19) / 20) - 1)^2
  }
  K                 = length(sa2)
  data$sa2          = sa2

  # precalculate --- > this step will double memory usage of X.
  if ( is.null(w)) {
      w        = colSums(X^2);
  }
  data$w       = w

  # initialize other parameters
  if ( is.null(pi) ) {
    if ( is.null(beta.init) ){

      Phi          = matrix(1,p,K)/K;
      pi           = rep(1,K)/K;

    } else {

      S            = outer(1/w, sa2, '+') * sigma2;
      Phi          = -data$beta^2/S/2 - log(S)/2;
      Phi          = exp(Phi - apply(Phi,1,max));
      Phi          = Phi / rowSums(Phi);
      pi           = colMeans(Phi);

    }
  } else {
    Phi          = matrix(rep(pi, each = p), nrow = p)
  }

  # run algorithm
  if ( is.null(update.order) ) {
    update.order   = 1:p
    if (method == "caisa") {
      if (update.pi) {
        out          = caisa_em   (X, w, sa2, pi, data$beta, r, sigma2,
                                   max.iter, min.iter, tol$convtol, tol$epstol,
                                   update.sigma, verbose)
      } else {
        out          = caisa_fix_pi(X, w, sa2, pi, data$beta, r, sigma2,
                                    max.iter, min.iter, tol$convtol, tol$epstol,
                                    update.sigma, verbose)
      }
    } else if (method == "accelerate") {
      out           = caisa_acc (X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, mixsqpiter, tol$convtol, tol$epstol,
                                 update.sigma, verbose)
    } else if (method == "block") {
      stepsize      = 1
      out           = caisa_g  (X, w, sa2, Phi, pi, data$beta, r, sigma2,
                                max.iter, min.iter, tol$convtol, tol$epstol,
                                stepsize, update.sigma, mode, verbose)
    } else if (method == "sigma_scaled") {
      out          = caisa_em2  (data$y, X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma, verbose)
      out$beta     = out$beta * sqrt(out$sigma2)
    } else if (method == "sigma_indep") {
      out          = caisa_em3  (X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma, verbose)
    }
  } else if (is.numeric(update.order)) {
    o   = rep(update.order - 1, max.iter)
    out = caisa_order(X, w, sa2, pi, data$beta, r, sigma2,
                      o, max.iter, min.iter, tol$convtol, tol$epstol,
                      update.sigma, verbose)
  } else if (update.order == "random") {
    update.order = random_order(p, max.iter)
    out = caisa_order(X, w, sa2, pi, data$beta, r, sigma2,
                      update.order, max.iter, min.iter, tol$convtol, tol$epstol,
                      update.sigma, verbose)
  }

  out$data = data

  out$update.order = update.order

  class(out) <- c("mr.ash", "list")

  out$data$X <- X

  return(out)
}

environment(mr.ashs) <- asNamespace('mr.ash.alpha')


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
#' start iteration from X1 (fix X2).
#' X1, X2, y need to be preprocessed consistently (all either centered/standardized or not)
#' @param mr.ash.init beta.init for mr.ash
#' - "lasso", init using lasso
#' - "lassoX1", init with lasso for X1, X2 with NULL
#' - "lassoX2", init with lasso for X2, X1 with NULL
#' - NULL, init with NULL for both X1 and X2.

mr.ash2s <- function(X1, X2, y, iter = 30, mr.ash.init = NULL){
  n <- dim(X1)[1]
  p1 <- dim(X1)[2]
  p2 <- dim(X2)[2]

  w1 <- colSums(X1^2)
  w2 <- colSums(X2^2)
  w <- c(w1, w2)
  gc()

  # init

  if (is.null(mr.ash.init)){

    X <- cbind(X1, X2)
    fit.init <- mr.ashs(X, y, w = w)

  } else if (mr.ash.init == "lasso"){

    X <- cbind(X1, X2) # this step needs optimization for memory usage, will double X.
    cvfit <- cv.glmnet(X, y) # this step needs optimization for memory usage, 4 times of X.
    b <- coef(cvfit, s = "lambda.min")[2: (p1 + p2 + 1), 1] # omit intercept
    rm(cvfit); gc()
    fit.init <- mr.ashs(X, y, w = w, beta.init = b)

  } else if (mr.ash.init == "lassoX1"){

    cvfit <- cv.glmnet(X1, y, nfolds = 5)
    b <- coef(cvfit, s = "lambda.min")[2: (p1 + 1), 1]
    rm(cvfit); gc()
    X <- cbind(X1, X2)
    fit.init <- mr.ashs(X, y, w = w, beta.init = c(b, rep(0, p2)))

  } else if (mr.ash.init == "lassoX2"){

    cvfit <- cv.glmnet(X2, y)
    b <- coef(cvfit, s = "lambda.min")[2: (p2 + 1), 1]
    rm(cvfit); gc()
    X <- cbind(X1, X2)
    fit.init <- mr.ashs(X, y, w = w, beta.init = c(rep(0, p1), b))

  }

  beta1 <- fit.init$beta[1: p1, , drop = F]
  beta2 <- fit.init$beta[(p1+1): (p1 + p2), , drop = F]

  pi1 <- fit.init$pi
  pi2 <- fit.init$pi

  sa2 <- fit.init$data$sa2

  sigma2 <- fit.init$sigma2

  rm(fit.init)
  if (!exists(X1)) {
    X1 <- X[, 1:p1]
    X2 <- X[, (p1 + 1): (p1 + p2)]
  }
  rm(X);gc()

  # start iteration (from X1)
  pi01_rec <- rep(0, iter)
  pi02_rec <- rep(0, iter)
  sigma2_rec <- rep(0, iter)

  for (i in 1:iter) {
    message("iter ", i)
    y1 <- y - X2 %*% beta2
    fit1 <- mr.ashs(X1, y1,
                   beta.init = beta1,
                   sigma2 = sigma2,
                   sa2 = sa2,
                   pi = pi1,
                   w = w1,
                   max.iter = 1)

    beta1 <- fit1$beta

    y2 <- y - X1 %*% beta1

    fit2 <- mr.ashs(X2, y2,
                   beta.init = beta2,
                   sigma2 = sigma2,
                   sa2 = sa2,
                   pi = pi2,
                   w = w2,
                   max.iter = 1)

    beta2 <- fit2$beta

    # update pi
    pi1 <- fit1$pi
    pi2 <- fit2$pi

    # update sigma
    r1 <- y1 - X1 %*% beta1
    s1 <- n * fit1$sigma2 - sum(r1^2) # s1 is d^2 (btilde - bhat) bhat for X1 in algorithm 4 in mr.ash paper.

    r2 <- y2 - X2 %*% beta2
    s2 <- n * fit2$sigma2 - sum(r2^2) # s2 is d^2 (btilde - bhat) bhat for X2 in algorithm 4 in mr.ash paper.

    r <- y - X1 %*% beta1 - X2 %*% beta2

    sigma2 <- 1/n * (sum(r^2) + s1 + s2)

    pi01_rec[i] <- pi1[1,]
    pi02_rec[i] <- pi2[1,]
    sigma2_rec[i] <- sigma2

    gc()
  }
  fit <- list("fit1" = fit1, "fit2" = fit2, "pi01" = pi01_rec, "pi02" = pi02_rec, "sigma2" = sigma2_rec)
  return(fit)
}
