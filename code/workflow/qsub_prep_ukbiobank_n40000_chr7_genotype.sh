
#!/bin/bash

#PBS -N bgen2plink
#PBS -S /bin/bash
#PBS -l mem=64gb
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=20

cd /home/szhao1/causal-TWAS/data; bash prep_ukbiobank_n40000_chr7_genotype.sh
