library(ctwas)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("5 arguments must be supplied:
       * genotype file name (multiple file names in a .txt file)
       * weight_1 file
       * weight_2 file
       * outfile name
       * outputdir", call.=FALSE)
}

pgenfs <- read.table(args[1], header = F, stringsAsFactors = F)[,1]

weight_1 <- args[2]
weight_2 <- args[3]

outname <- args[4]
outputdir <- args[5]

weight <- c(weight_1, weight_2)

exprfs <- impute_expr(pgenfs = pgenfs, weight = weight,
                           method = "lasso", outputdir = outputdir,
                           outname = outname)

out <- data.frame("file" = exprfs, stringsAsFactors = F)
write.table(out, file = paste0(outputdir, "/", outname, ".expr.txt"),
             row.names=F, col.names=F, sep="\t", quote = F)
