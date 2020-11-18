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
      s$Xr = s$Xr - compute_Xb(X,s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      R = Y - s$Xr

      res = single_effect_regression(R,X,s$V[l],s$sigma2,s$pi,
                                     estimate_prior_method,
                                     check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l,]      = res$V  # different from susieR
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$loglik +
        SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
                               res$alpha * res$mu2)

      s$Xr = s$Xr + compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }
  return(s)
}

assignInNamespace("update_each_effect", update_each_effect,"susieR")

# Update a susie fit object in order to initialize susie model.
init_finalize = function (s, X = NULL, Xr = NULL) {
  if(length(s$V) == 1)
    s$V = rep(s$V,nrow(s$alpha))

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
    s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}

assignInNamespace("init_finalize", init_finalize ,"susieR")
