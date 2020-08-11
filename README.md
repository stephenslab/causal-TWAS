# causal-TWAS

New method to identify causal genes in TWAS studies. 


## Required libraries

* plink2R
```
wget https://github.com/gabraham/plink2R/archive/master.zip
unzip master.zip
install.packages('plink2R-master/plink2R/',repos=NULL)
```

* package not managed by conda
mr.ash.alpha: `devtools::install_github("stephenslab/mr.ash.alpha@1acabf8a5032be3b13a1707b7fe514016d3ae0de")`
bigstatsr, bigreadr, biglasso: installed from CRAN.

## About memory usage

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

## Data class

Change from `matrix` type in R to `FBM` type in `bigstatr` package. Advantage: 1. can apply `biglasso` on it: save memory for lasso step. 2. Call by reference, so this will not copy the matrix. 3. easy connect to Rcpp. 4. file backed, no need to load everything to memory.

After switch to FBM, all parts should not limit by memory except for mr.ash. However, when using multiple cores, and core can't read simutaneously from the same back up file. 

1. For `bigstatsr::big_apply`, only use one core. when using multiple cores, even if they are from the same node, they will interfere with each other.

2. For `biglasso`, only one node can access the .bk file at one time. can use multiple core from the same node. 

When writing code, pay attention to 1. and when running, use multiple copies of genotype file to avoid reading from the same file.

simulate phenotype: need to change for FBM class type. 

Some other options: `snpnet` for lasso and `pgen` to read data. Problem: need to start from pgen file, not sure how to add or delete features/columns. 


## Next steps

1. Simulation process:

To simulate phenotype, we use fusion's weights ($\alpha$), these weights were learned from GTEx or other dataset with both expression and genotype data. Given a large cohort with genotype data (like ukbionank or hapgen2 ), we then simulate $Y$ following 
$$Y= G\alpha\beta + G\theta+ \epsilon$$ ,where $G$ is genotype, $\beta$ and $\theta$ follows distribution described in the [sparse model](sparse_model.html). 

To simulate expression, we choose a separate cohort with genotype data and $X =  G\alpha + \xi$. 

To run our program, we train expression model to obtain $\alpha$ (note this will be different from the true $alpha$, as it is estimated), then impute expression in the large genotype-phenotype cohort.
    
2. MR.ash full version workflow

initialize: use lasso
grid: not that important
rPIP: be strict, try local false sign rate
susie version: seems w0 is the best

sa2 --> make it to an argument?
rPIP --> make it to an arugment?
s80k: sigma2 infer still has some problem.


3. Try mr.ash + iterative version of susie
3.1 mr.ash for all genes





# Known bugs:
Currently, snp ref and alt can be different between Fusion weights and ukbiobank, the definition of ref and alt comes from plink bim files in GTEx and pvar files in UKBB, they are different. The snp name, chr and position matches in general.

