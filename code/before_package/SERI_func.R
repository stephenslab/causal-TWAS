#' @title Bayesian single-effect linear regression of y on X,
#' adapted from susieR (single_effect_regression.R)
#'
#' @param Y An n vector.
#'
#' @param X An n by p matrix of covariates. centered and scaled.
#'
#' @param V The prior variance.
#'
#' @param residual_variance The residual variance.
#'
#' @param prior_weights A p vector of prior weights.
#'
#' @return a list
#'
#' @importFrom Matrix colSums
#'
#'
get_lbf =
  function (Y, X, V, residual_variance = 1) {

    ytX = crossprod(Y, X)
    Xty = t(ytX)
    xTx = colSums(X * X)
    betahat = (1/xTx) * Xty
    shat2 = residual_variance/xTx


    # log(bf) for each SNP
    lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
      dnorm(betahat,0,sqrt(shat2),log = TRUE)

    # Deal with special case of infinite shat2 (e.g., happens if X does
    # not vary).
    lbf[is.infinite(shat2)] = 0

    z = betahat/sqrt(shat2)
    z[is.infinite(z)] = 0

    pvalue = 2*pnorm(-abs(z))

    lbfout <- cbind(lbf, betahat, shat2, z, pvalue)
    colnames(lbfout) <- c("lbf", "betahat", "shat2", "z", "pvalue")

    return(lbfout)
  }


#' @title Bayesian single-effect linear regression of y on X,
#' adapted from susieR (single_effect_regression.R)
#'
#' @param lbf a vector
#'
#' @return A list with the following elements:
#'
#' \item{alpha}{Vector of posterior inclusion probabilities;
#'   \code{alpha[i]} is posterior probability that that the ith
#'   coefficient is non-zero.}
#'
#' \item{lbf_model}{Log-Bayes factor for the single effect regression.}
#'
#' \item{loglik}{The logarithm of the likelihood, \eqn{p(y | X, V)}.}
#'
#' @importFrom Matrix colSums
#'
#'
SER =
  function (lbf, prior_weights = NULL, null_weight = NULL) {
    if (is.null(prior_weights) & is.null(null_weight))
      prior_weights = c(rep(1/length(lbf), length(lbf)))

    if (is.null(prior_weights) & !is.null(null_weight))
      prior_weights = c(rep(1/length(lbf)*(1-null_weight), length(lbf)))

    if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
    if (!is.null(null_weight)) {
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
      lbf <- c(lbf,0)
    }

    # w is proportional to BF, but subtract max for numerical stability.
    maxlbf = max(lbf)
    w = exp(lbf - maxlbf)

    # Posterior prob for each SNP.
    w_weighted = w * prior_weights
    weighted_sum_w = sum(w_weighted)
    alpha = w_weighted/weighted_sum_w

    if (!is.null(null_weight)) {
      alpha <- alpha[-length(lbf)]
    }

    # BF for single effect model.
    lbf_model = maxlbf + log(weighted_sum_w)
    loglik = lbf_model

    return(list(alpha = alpha,
                lbf_model = lbf_model,loglik = loglik))
  }



