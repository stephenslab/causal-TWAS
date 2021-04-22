setwd("/home/simingz/causalTWAS/GTEx_v7_eQTL")
lookupfile <- "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz"
egenefile <- "Adipose_Subcutaneous.v7.egenes.txt.gz"
eqtlfile <- "Adipose_Subcutaneous.v7.signif_variant_gene_pairs.txt.gz"

library(data.table)
eqtl <- fread(eqtlfile)
ref <- fread(lookupfile)
a <- merge(eqtl, ref, all.x = T, by = "variant_id")
rm(ref, eqtl)
egene <- fread(egenefile)
egene <- egene[, 1:11, drop =F]
out <- merge(a, egene, all.x = T, by = "gene_id")
write.table(out , file= "Adipose_Subcutaneous.v7.eQTL.txt" , row.names=F, col.names=T, sep="\t", quote = F)
