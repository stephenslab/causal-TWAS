# fit gene first, then with SNP
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
