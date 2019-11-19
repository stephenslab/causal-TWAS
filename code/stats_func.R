# Short script to center and scale the columns (i.e., SNPs) of a
# matrix containing genotypes. Credit to Peter
library(matrixStats)
library(Rcpp)

sourceCpp(paste0(codedir,"scale.cpp"))

scaleRcpp <- function(x) {
  cat("Scaling genotype data.\n")
  # Load the genotype data as a matrix of floating-point numbers.
  storage.mode(x) <- "double"
  # Remove all SNPs that do not vary.
  # cat("Filtering SNPs.\n")
  # s    <- colSds(x)
  # x <- x[,s > 0]
  # gc()
  # center and scale
  timing <- system.time({
    x.scaled <- x
    mu <- colMeans(x)
    s  <- colSds(x)
    scale_rcpp(x.scaled,mu,s)
  })
  print(timing)
  # verify
  cat("Get the largest & smallest column mean:\n")
  print(range(colMeans(x.scaled)))
  cat("Get the largest & smallest column s.d.:\n")
  print(range(colSds(x.scaled)))
  return(x.scaled)
}
