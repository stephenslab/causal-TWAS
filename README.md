# causal-TWAS analysis

New method to identify causal genes in TWAS studies. 

# Run a workflow

Dry run first:
```
snakemake -s ~/causalTWAS/causal-TWAS/code/workflow/Snakefile-test -n --dag | dot -Tsvg > dag.svg
```

Run on node:
```
snakemake -s ~/causalTWAS/causal-TWAS/code/workflow/Snakefile-test -j 1
```

Run on cluster:

```
snakemake -s ~/causalTWAS/causal-TWAS/code/workflow/Snakefile-test -j 999 --cluster-config cluster-broadwl.json --cluster "sbatch -n {cluster.n}  -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.out} --error {cluster.err} --partition {cluster.partition}"
```