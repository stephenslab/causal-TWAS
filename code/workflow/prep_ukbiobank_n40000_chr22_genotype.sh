cd /home/szhao1/uk-biobank/uk-biobank-genotypes; /home/szhao1/uk-biobank/uk-biobank-phenotypes/tools/plink2 --sample /home/szhao1/uk-biobank/uk-biobank-phenotypes/data/genotypes/ukb27386_imp_v3_s487324.sample --keep /home/szhao1/causal-TWAS/data/ukbiobank_samples40000.txt --bgen ukb_imp_chr22_v3.bgen --chr 22 --threads 19 --make-pgen --out ~/causal-TWAS/data/ukb_chr22_s40000 --maf 0.05 --geno 0.05 --memory 30000