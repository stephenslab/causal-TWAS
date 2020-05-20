# simulate data with two groups of effects
# Two groups of variables: k1 and k2 non-zero effects with effect s.d. sigma1/sigma2.
# Group 1: SNP effect, Group 2: gene effect.
# Variables in group 2 are in “complete LD” with some of the varaibles in group 1.
n <- 3000
p <- 1000
p2 <- 100
p1 <- p - p2
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
#X2 <- matrix(rnorm(n*p2), nrow=n, ncol=p2)
#X1 <- cbind(matrix(rnorm(n*(p1-p2)), nrow=n, ncol=p1-p2), X2)
#X <- cbind(X1, X2)
k1 <- 20
sigma1 <- 0.1
k2 <- 20
sigma2 <- 0.2
beta1 <- c(rnorm(k1, 0, sigma1), rep(0, p1-k1))
beta2 <- c(rnorm(k2, 0, sigma2), rep(0, p2-k2))
beta <- c(beta1, beta2)
sigma.e <- 1
y <- X %*% beta + rnorm(n, 0, sigma.e)
# fit a single group
fit <- mr.ash(X, y, method="caisa")
beta.single <- fit$beta
which(abs(beta) > 0.05)
which(abs(beta.single) > 0.05)
# iteratively fit two groups using mr.ash
X1 <- X[, 1:p1]
X2 <- X[, (p1+1):p]
r1 <- y
r2 <- y
fit2 <- mr.ash(X2, r2, method = "caisa")
beta2.pm <- fit2$beta
r1 <- y - X2 %*% beta2.pm
fit1 <- mr.ash(X1, r1, method="caisa")
beta1.pm <- fit1$beta
r2 <- y - X1 %*% beta1.pm
beta.pm <- c(beta1.pm, beta2.pm)
which(abs(beta.pm) > 0.05)