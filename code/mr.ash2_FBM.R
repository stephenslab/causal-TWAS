library(mr.ash.alpha)
library(biglasso)

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
#' @param X1, matrix
#' @param X2, matrix

mr.ash2s_iter <- function(X1, X2,
                          y,
                          w1 = NULL , w2 = NULL,
                          beta1 = NULL, beta2 = NULL,
                          pi1 = NULL, pi2 = NULL,
                          sigma2 = NULL,
                          iter = 30){
  n <- dim(X1)[1]
  p1 <- dim(X1)[2]
  p2 <- dim(X2)[2]

  if (is.null(w1)){
    w1 <- colSums(X1^2)
  }

  if (is.null(w2)){
    w2 <- colSums(X2^2)
  }

  init.par <- tibble::lst(beta1, beta2, pi1, pi2, sigma2, iter)

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
                    pi = pi1,
                    w = w1,
                    max.iter = 1)

    beta1 <- fit1$beta

    y2 <- y - X1 %*% beta1

    fit2 <- mr.ashs(X2, y2,
                    beta.init = beta2,
                    sigma2 = sigma2,
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

  fit <- tibble::lst(fit1, fit2, pi01_rec, pi02_rec, sigma2_rec, init.par)
  return(fit)
}


#' Provide initiation parameters for mr.ash2s.
#' Run mr.ash for all data (snp and expr) and use the
#' parameters as initiating parameters for mr.ash2s
#' @param snp FBM
#' @param expr matrix
#' @param mr.ash.init beta.init for mr.ash
#' - "NULL", init mr.ash using default
#' - "lasso", init mr.ash using lasso
#' - "lassoSNP", init with lasso for SNP and NULL for gene
#' @param init.order "es", expr-snp; "se", "snp-expr"
#' @param iter.order "es", expr-snp; "se", "snp-expr"
#'
mr.ash2s <- function(expr,
                     snp,
                     y,
                     mr.ash.init = c(NULL, "lasso", "lassoSNP"),
                     init.order = c("es", "se"),
                     iter.order =c("es", "se"),
                     outname = "mr.ash2s",
                     backingfile = outname,
                     ncores = 1, iter = 30
                     ){

  mr.ash.init <- match.arg(mr.ash.init)
  init.order <- match.arg(init.order)
  iter.order <- match.arg(iter.order)

  n <- snp$nrow
  p.expr <- dim(expr)[2]
  p.snp <- snp$ncol

  # prepare X
  if (file.exists(paste0(backingfile,".bk"))) {
    X <- readRDS(paste0(backingfile, ".rds"))
  } else {
    X <- as_FBM(expr, backingfile = backingfile)$save()
    X$add_columns(p.snp)

    a <- big_apply(snp, a.FUN = function(m, ind) {
      X[, ind + p.expr] <<- m[, ind, drop = FALSE]; return(0)
    }, a.combine = 'c', block.size = 2000, ncores = ncores)
  }

  # prepare w
  print("calculating w")
  w.snp <- big_apply(snp, a.FUN = function(m, ind) {
    colSums(m[, ind, drop = FALSE]^2)
  }, a.combine = 'c', block.size = 2000, ncores = ncores)

  w.expr <- colSums(expr^2)
  w = c(w.expr, w.snp)

  # prepare mr.ash beta.init
  print("prepare mr.ash.init")
  if (is.null(mr.ash.init)){
    b <- NULL
  } else if (mr.ash.init == "lasso"){
    X.bm <- X$bm()
    time.cvfit <- system.time(
      cvfit <- cv.biglasso(X.bm, y, screen = 'SSR-BEDPP', seed = 1234, ncores = ncores, nfolds = 5)
    )
    save(cvfit, time.cvfit, file = paste0(outname,".lasso-beta.Rd"))
    b <- coef(cvfit)
    b <- b[-1] # discard intercept
  } else if (mr.ash.init == "lassoSNP"){
    snp.bm <- snp$bm()
    time.cvfit <- system.time(
      cvfit <- cv.biglasso(snp.bm, y, screen = 'SSR-BEDPP', seed = 1234, ncores = ncores, nfolds = 5)
    )
    save(cvfit, time.cvfit, file = paste0(outname,".lassoSNP-beta.Rd"))
    b.snp <- coef(cvfit)
    b <- c(rep(0, p.expr), b.snp[-1])
  }

  # prepare order
  if (init.order == "es"){
    o <- 1: (p.expr + p.snp)
  } else if ( init.order == "se"){
    o <- c(1: p.expr + p.snp, 1: p.snp)
  }

  # mr.ash for all data as initiation
  print("run mr.ash to initiate")
  X.m <- as.matrix(X$bm())
  fit.init <- mr.ashs(X.m, y, w = w, beta.init = b, update.order = o)

  # start interation
  beta.expr <- fit.init$beta[1: p.expr, , drop = F]
  beta.snp <- fit.init$beta[(p.expr + 1): (p.expr + p.snp), , drop = F]

  pi.expr <- pi.snp <- fit.init$pi

  sigma2 <- fit.init$sigma2

  rm(fit.init, X.m); gc()

  if (iter.order == "es"){

    mr.ash2s.fit <- mr.ash2s_iter(
      X1 = expr, X2 = as.matrix(snp$bm()),
      y = y,
      w1 = w.expr, w2 = w.snp,
      beta1 = beta.expr, beta2 = beta.snp,
      pi1 = pi.expr, pi2 = pi.snp,
      sigma2 = sigma2,
      iter = iter)

  } else if (iter.order == "se"){

    mr.ash2s.fit <- mr.ash2s_iter(
      X2 = expr, X1 = as.matrix(snp$bm()),
      y = y ,
      w2 = w.expr, w1 = w.snp,
      beta2 = beta.expr, beta1 = beta.snp,
      pi2 = pi.expr, pi1 = pi.snp,
      sigma2 = sigma2,
      iter = iter)
  }

  mr.ash2s.fit$init.order <- init.order
  mr.ash2s.fit$mr.ash.init <- mr.ash.init
  mr.ash2s.fit$iter.order <- iter.order

  if (iter.order == "es"){
    names(mr.ash2s.fit)[1] <- "g.fit"
    names(mr.ash2s.fit)[2] <- "s.fit"
  } else if (iter.order == "se"){
    names(mr.ash2s.fit)[1] <- "s.fit"
    names(mr.ash2s.fit)[2] <- "g.fit"
  }

  mr.ash2s.fit$s.fit$data$X <- snp

  return(mr.ash2s.fit)
}

