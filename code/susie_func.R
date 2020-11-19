# @title update each effect once
# @param X an n by p matrix of regressor variables
# @param Y an n vector of response variable
# @param s a SuSiE fit
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
update_each_effect = function (X, Y, s, estimate_prior_variance = FALSE,
                               estimate_prior_method = "optim",
                               check_null_threshold) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  L = nrow(s$alpha)
  if (L > 0)
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$Xr = s$Xr - susieR:::compute_Xb(X,s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      R = Y - s$Xr

      res = susieR:::single_effect_regression(R,X,s$V[l,],s$sigma2,s$pi,
                                     estimate_prior_method,
                                     check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l,]      = res$V  # different from susieR
      s$lbf[l]    = res$lbf_model
      s$KL[l]     = -res$loglik +
        susieR:::SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
                               res$alpha * res$mu2)

      s$Xr = s$Xr + susieR:::compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }
  return(s)
}

assignInNamespace("update_each_effect", update_each_effect, "susieR")

# Update a susie fit object in order to initialize susie model.
init_finalize = function (s, X = NULL, Xr = NULL) {
  # different form susieR
  # if(length(s$V) == 1)
  #   s$V = rep(s$V, nrow(s$alpha))

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")

  # different from susieR
  # if (nrow(s$alpha) != length(s$V))
  #   stop("Input prior variance V must have length of nrow of alpha in ",
  #        "input object")

  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = susieR:::compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}

assignInNamespace("init_finalize", init_finalize ,"susieR")


susie_get_cs <- function (res, X = NULL, Xcorr = NULL, coverage = 0.95, min_abs_corr = 0.5,
                          dedup = TRUE, squared = FALSE)
{
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
    stop("Xcorr matrix must be symmetric")
  }
  if (inherits(res, "susie")) {
    null_index = res$null_index
    if (is.numeric(res$V))
      include_idx = rep(TRUE, nrow(res$alpha)) # different from susieR
    else include_idx = rep(TRUE, nrow(res$alpha))
  }
  else null_index = 0
  status = susieR:::in_CS(res$alpha, coverage)
  cs = lapply(1:nrow(status), function(i) which(status[i, ] !=
                                                  0))
  include_idx = include_idx * (lapply(cs, length) > 0)
  if (dedup)
    include_idx = include_idx * (!duplicated(cs))
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL, coverage = coverage))
  cs = cs[include_idx]
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) = paste0("L", which(include_idx))
    return(list(cs = cs, coverage = coverage))
  }
  else {
    purity = data.frame(do.call(rbind, lapply(1:length(cs),
                                              function(i) {
                                                if (null_index > 0 && null_index %in% cs[[i]])
                                                  c(-9, -9, -9)
                                                else susieR:::get_purity(cs[[i]], X, Xcorr, squared)
                                              })))
    if (squared)
      colnames(purity) = c("min.sq.corr", "mean.sq.corr",
                           "median.sq.corr")
    else colnames(purity) = c("min.abs.corr", "mean.abs.corr",
                              "median.abs.corr")
    threshold = ifelse(squared, min_abs_corr^2, min_abs_corr)
    is_pure = which(purity[, 1] >= threshold)
    if (length(is_pure) > 0) {
      cs = cs[is_pure]
      purity = purity[is_pure, ]
      row_names = paste0("L", which(include_idx)[is_pure])
      names(cs) = row_names
      rownames(purity) = row_names
      ordering = order(purity[, 1], decreasing = T)
      return(list(cs = cs[ordering], purity = purity[ordering,
                                                     ], cs_index = which(include_idx)[is_pure[ordering]],
                  coverage = coverage))
    }
    else {
      return(list(cs = NULL, coverage = coverage))
    }
  }
}

assignInNamespace("susie_get_cs", susie_get_cs ,"susieR")



susie_get_pip <- function (res, prune_by_cs = FALSE, prior_tol = 1e-09)
{
  if (inherits(res, "susie")) {
    if (res$null_index > 0)
      res$alpha = res$alpha[, -res$null_index, drop = FALSE]
    if (is.numeric(res$V))
      include_idx = 1:nrow(res$alpha) # different from susieR
    else include_idx = 1:nrow(res$alpha)
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx, res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)
    if (length(include_idx) > 0) {
      res = res$alpha[include_idx, , drop = FALSE]
    }
    else {
      res = matrix(0, 1, ncol(res$alpha))
    }
  }
  return(as.vector(1 - apply(1 - res, 2, prod)))
}
assignInNamespace("susie_get_pip", susie_get_pip ,"susieR")
