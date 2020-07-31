# causal-TWAS

New method to identify causal genes in TWAS studies. 


## Required libraries

* plink2R
```
wget https://github.com/gabraham/plink2R/archive/master.zip
unzip master.zip
install.packages('plink2R-master/plink2R/',repos=NULL)
```

## Next steps

1. Use UK biobank genotype data or [hapgen2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html), may need n~10,000.

2. Simulation process:
To simulate phenotype, we use fusion's weights ($\alpha$), these weights were learned from GTEx or other dataset with both expression and genotype data. Given a large cohort with genotype data (like ukbionank or hapgen2 ), we then simulate $Y$ following 
$$Y= G\alpha\beta + G\theta+ \epsilon$$ ,where $G$ is genotype, $\beta$ and $\theta$ follows distribution described in the [sparse model](sparse_model.html). 

    To simulate expression, we choose a separate cohort with genotype data and $X =  G\alpha + \xi$. 

    To run our program, we train expression model to obtain $\alpha$ (note this will be different from the true $alpha$, as it is estimated), then impute expression in the large genotype-phenotype cohort.

3. Try mr.ash
Basic idea of mr.ash:You want to add up 2 linear regressions, which we might call “learners”.
A learner is, roughly, a fucntion, that learns a relationshipo between Y and some predictors X
(here the Xs are different for your 2 different learners) Let’s write this $Y=L1(X1) + L2(X2) + E$
. The idea is that you can fit this model iteratively, by first applying $L1$, and then applying $L2$, each time applying it to the appropriate residuals.So you apply $L1$ to learn a relationship between $R=Y-L2$ (residuals) and $X1$, and then $L2$ to learn a relationship between $R=Y-L1$ and $X2$ and you just iterate….
probably best to start by trying the mr.ash.alpha package. It would be great if you can follow the instructions in https://stephenslab.github.io/mr-ash-workflow/ to start reproducing some of the results in the paper.

Try mr.ash using real genotype data. We can try 0.5 heritability. If using snps on one chromosome, we may need a few thousand samples, if using whole genome, we may need 10,000 samples. 
    
4. About memory usage

R `matrix`, type `double`

ukb chr22, 20000 samples, Rd file (0.9G)

* simulate data

  - `load` need 13.8G (`dat` is 13.8G by `object.size()`).
  - `scale` 25.8G.
  - `cis_expr` 15G
  - `simulate_phenotype` 15G

* mr.ash2s
mr.ash2s without lasso uses 15G, cv.glmnet uses 59G.

`bigstatr`, type `FBM`, `unsigned short`
backup file and in R: 3G. 

5. data class

Change from `matrix` type in R to `FBM` type in `bigstatr` package. Advantage: 1. can apply `biglasso` on it: save memory for lasso step. 2. Call by reference, so this will not copy the matrix. 3. easy connect to Rcpp. 4. file backed, no need to load everything to memory.

After switch to FBM, all parts should not limit by memory except for mr.ash. However, when using multiple cores, and core can't read simutaneously from the same back up file. 
1. For `bigstatsr::big_apply`, only use one core. when using multiple cores, even if they are from the same node, they will interfere with each other 
2. For `biglasso`, only one node can access the .bk file at one time. can use multiple core from the same node. 

When writing code, pay attention to 1. and when running, use multiple copies of gentype file to avoid reading from the same file.

simulate phenotype: need to change for FBM class type.

Some other options: `snpnet` for lasso and `pgen` to read data. problem: need to start from pgen file, not sure how to add or delete features/columns. 

# packages not managed by conda:
mr.ash.alpha: `devtools::install_github("stephenslab/mr.ash.alpha@1acabf8a5032be3b13a1707b7fe514016d3ae0de")`
bigstatsr
bigreadr
biglasso


# Known bugs:
Currently, snp ref and alt can be different between Fusion weights and ukbiobank, the definition of ref and alt comes from plink bim files in GTEx and pvar files in UKBB, they are different. The snp name, chr and position matches in general.

