A Quick Start of cola Package
=============================================================

**Author**: Zuguang Gu ( z.gu@dkfz.de )

**Date**: 2021-07-14

**Package version**: 1.8.0

-------------------------------------------------------------






## Consensus partitioning

Assume your matrix is stored in an object called `mat`, to perform consensus
partitioning with *cola*, you only need to run following code:


```r
# code only for demonstration
mat = adjust_matrix(mat)  # optional
rl = run_all_consensus_partition_methods(mat, cores = ...)
cola_report(rl, output_dir = ..., cores = ...)
```

In above code, there are three steps:

1. Adjust the matrix. In this step, rows with too many `NA`s are removed. Rows
   with very low variance are removed. `NA` values are imputed if there are
   less than 50% in each row. Outliers are adjusted in each row.
2. Run consensus partitioning with several methods.
   Partitioning methods are `hclust` (hierarchical clustering with cutree),
   `kmeans` (k-means clustering), `skmeans::skmeans` (spherical k-means
   clustering), `cluster::pam` (partitioning around medoids) and
   `Mclust::mclust` (model-based clustering). The default methods to extract
   top n rows are `SD` (standard deviation), `CV` (coefficient of variation),
   `MAD` (median absolute deviation) and `ATC` (ability to correlate to other
   rows). 
3. Generate a detailed HTML report for the complete analysis.


`run_all_consensus_partition_methods()` runs multiple methods in sequence, which might
take long time for big datasets. Users can also run consensus partitioining with
a specific top-value methods (e.g. SD) and partitioning methods (e.g. skmeans) by 
`consensus_partition()` function:


```r
res = consensus_partition(mat, top_value_method = ..., partition_method = ...)
cola_report(res, output_dir = ..., cores = ...)
```

You can refer to [the main vignette](cola.html) for more details.

For extremely large datasets, users can run `consensus_partition_by_down_sampling()` by randomly 
sampling a subset of samples for classification, later the classes of the remaining
samples are predicted by the signatures of the _cola_ classification. More details
can be found in the vignette ["Work with Big Datasets"](work_with_big_datasets.html).


```r
res = consensus_partition_by_down_sampling(mat, subset = ...,
    top_value_method = ..., partition_method = ...)
cola_report(res, output_dir = ..., cores = ...)
```

There are examples on real datasets for _cola_ analysis that can be found at https://jokergoo.github.io/cola_collection/.

## Hierarchical consensus partitioning

From version 2.0.0, there is a new function `hierarchical_partition()` that applies consensus partitioning in
a hierarchical way. Simply use `hierarchical_partition()` with the matrix:


```r
rh = hierarchical_partition(mat, cores = ...)
cola_report(rh, output_dir = ..., cores = ...)
```

With big matrix, argument `subset` can be set so that down sampling consensus partitioning will be internally used. E.g.


```r
rh = hierarchical_partition(mat, subset = 500, cores = ...)
```

Please refer to the vignette ["Hierarchical Consensus Partitioning"](hierarchical.html) for more details on this method.

