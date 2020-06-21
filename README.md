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

ukb chr22, 20000 samples, Rd file (0.9G)

* simulate data

  - `load` need 13.8G (`dat` is 13.8G by `object.size()`).
  - `scale` 25.8G.
  - `cis_expr` 15G
  - `simulate_phenotype` 15G

* mr.ash2s
27


