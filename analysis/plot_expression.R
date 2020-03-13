exprlist <- list()
for (gname in names(weightslist)){
  weights <- weightslist[[gname]]
  exprlist[[gname]] <- as.matrix(X.scaled[ , match(weights$labels, labels)]) %*% weights$alpha
  if (is.na(exprlist[[gname]])) exprlist[[gname]] <- 0 #TODO:eQTL not genotyped should be omited, instead of simply put gene expression to 0.
}

for (i in c(1:20)){
  print(i)
  plot(expr[,i],implist[[i]], main=names(implist)[i], xlab = "true expression" , ylab =" imputed expression")
  abline(coef = c(0,1), col="red")
}
