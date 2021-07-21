Predict Classes for New Samples
=============================================================

**Author**: Zuguang Gu ( z.gu@dkfz.de )

**Date**: 2021-07-14

**Package version**: 1.8.0

-------------------------------------------------------------




`predict_classes()` function predicts sample classes based on cola
classification. It is mainly used in the following two scenarios: 

1. Predict sample classes in a new study.
2. For a dataset with huge number of samples, users can apply *cola* on a randomly 
   sampled subset of samples and predict classes for the remaining samples.

To use a *cola* classification, users need to select the result from a specific
top-value method and partitioning method, i.e., the input should be a
`ConsensusPartition` object. In the following example, we use the analysis
from [Golub dataset](https://jokergoo.github.io/cola_examples/Golub_leukemia/)
and we take the result from the method `ATC:skmeans`.


```r
data(golub_cola)
res = golub_cola["ATC:skmeans"]
res
```

```
## A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
##   On a matrix with 4116 rows and 72 columns.
##   Top rows (412, 824, 1235) are extracted by 'ATC' method.
##   Subgroups are detected by 'skmeans' method.
##   Performed in total 750 partitions by row resampling.
##   Best k for subgroups seems to be 3.
## 
## Following methods can be applied to this 'ConsensusPartition' object:
##  [1] "cola_report"             "collect_classes"         "collect_plots"          
##  [4] "collect_stats"           "colnames"                "compare_partitions"     
##  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
## [10] "functional_enrichment"   "get_anno"                "get_anno_col"           
## [13] "get_classes"             "get_consensus"           "get_matrix"             
## [16] "get_membership"          "get_param"               "get_signatures"         
## [19] "get_stats"               "is_best_k"               "is_stable_k"            
## [22] "membership_heatmap"      "ncol"                    "nrow"                   
## [25] "plot_ecdf"               "predict_classes"         "rownames"               
## [28] "select_partition_number" "show"                    "suggest_best_k"         
## [31] "test_to_known_factors"   "top_rows_heatmap"
```

`predict_classes()` needs at least three arguments: a `ConsensusPartition`
object, the number of subgroups and the new matrix. The first two are for
extracting the classification as well as the signatures that best separate
subgroups. The new matrix should have the same number of rows as the matrix
used for *cola* analysis, also the row orders should be the same. **Be careful
that the scaling of the new matrix should also be the same as the one applied
in cola analysis.**

The prediction is based on **the signature centroid matrix**. The processes
are as follows:

1. For the provided `ConsensusPartition` object and a selected k, the
   signatures that discriminate classes are extracted by `get_signatures()`.
   If number of signatures is more than 2000, only 2000 signatures are randomly
   sampled.
2. The signature centroid matrix is a k-column matrix where each column is the
   centroid of samples in the corresponding class, i.e., the mean across
   samples. If rows were scaled in *cola* analysis, the signature centroid
   matrix is the mean of scaled values and vise versa. Please note the samples
   with silhouette score less than `silhouette_cutoff` (0.5 as the default)
   are removed for calculating the centroids.
3. With the signature centroid matrix and the new matrix, it perform the
   prediction.

The class prediction is applied as follows. For each sample in the new matrix,
the task is basically to test which signature centroid the current sample is
the closest to. There are three methods: the Euclidean distance, cosine distance (it is `1-cos(x, y)`) and the
correlation (Spearman) distance. 

For the Euclidean distance and cosine distance method, for the vector denoted as $x$ which
corresponds to sample $i$ in the new matrix, to test which class should be
assigned to sample $i$, the distance between sample $i$ and all $k$ signature
centroids are calculated and denoted as $d_1$, $d_2$, ..., $d_k$. The class
with the smallest distance is assigned to sample $i$.

To test whether the class assignment is statistically significant, or to test
whether the class that is assigned is significantly closest to sample $i$, we
design a statistic named "_difference ratio_", denoted as $r$ and calculated as
follows. First, The distances for $k$ centroids are sorted increasingly, and we 
calculate $r$ as: 

$$r = \frac{|d_{(1)} - d_{(2)}|}{\bar{d}}$$

which is the difference between the smallest distance and the second smallest
distance, normalized by the mean distance. To test the statistical
significance of $r$, we randomly permute rows of the signature centroid matrix
and calculate $r_{rand}$. The random permutation is performed 1000 times and the p-value is calculated as the proportion of
$r_{rand}$ being larger than $r$.

For the correlation method, the distance is calculated as the Spearman
correlation between sample $i$ and signature centroid $k$. The label for the
class with the maximal correlation value is assigned to sample $i$. The
p-value is simply calculated by `stats::cor.test()` between sample $i$ and
centroid $k$.

If a sample is tested with a p-value higher than 0.05, the
corresponding class label is set to `NA`.

To demonstrate the use of `predict_classes()`, we use the same matrix as for *cola*
analysis.


```r
mat = get_matrix(res)
```

Note the matrix was row-scaled in *cola* analysis, thus, `mat` should also be scaled
with the same method (z-score scaling).


```r
mat2 = t(scale(t(mat)))
```

And we predict the class of `mat2` with the 3-group classification from `res`.



```r
cl = predict_classes(res, k = 3, mat2)
```

```
## The matrix has been scaled in cola analysis, thus the new matrix should also be scaled
## with the same method ('z-score'). Please double check.
## Set `help = FALSE` to suppress this message. 
## 
## * 70/72 samples (in 3 classes) remain after filtering by silhouette (>= 0.5).
## * cache hash: 5e81712ade334cd09013d9ac870cab7e (seed 888).
## * calculating row difference between subgroups by Ftest.
## * split rows into 4 groups by k-means clustering.
## * 2076 signatures (50.4%) under fdr < 0.05, group_diff > 0.
## * Predict classes based on 3-group classification (euclidean method) on a 72-column matrix.
```

<img src="figure/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="600" style="display: block; margin: auto;" />

```r
cl
```

```
##           class     p
## sample_39     3 0.000
## sample_40     3 0.000
## sample_42     2 0.000
## sample_47     2 0.000
## sample_48     2 0.000
## sample_49     3 0.000
## sample_41     2 0.000
## sample_43     2 0.000
## sample_44     2 0.000
## sample_45     2 0.000
## sample_46     2 0.000
## sample_70     2 0.000
## sample_71     2 0.000
## sample_72     2 0.000
## sample_68     2 0.000
## sample_69     2 0.000
## sample_67     2 0.000
## sample_55     3 0.000
## sample_56     3 0.000
## sample_59     2 0.000
## sample_52     1 0.000
## sample_53     1 0.000
## sample_51     1 0.000
## sample_50     1 0.000
## sample_54     2 0.000
## sample_57     1 0.000
## sample_58     1 0.000
## sample_60     1 0.058
## sample_61     1 0.000
## sample_65     1 0.000
## sample_66     2 0.000
## sample_63     1 0.000
## sample_64     1 0.000
## sample_62     1 0.000
## sample_1      3 0.000
## sample_2      2 0.000
## sample_3      3 0.000
## sample_4      3 0.000
## sample_5      2 0.000
## sample_6      3 0.000
## sample_7      3 0.000
## sample_8      3 0.000
## sample_9      2 0.000
## sample_10     3 0.000
## sample_11     2 0.000
## sample_12     1 0.000
## sample_13     2 0.000
## sample_14     2 0.000
## sample_15     2 0.000
## sample_16     2 0.000
## sample_17     2 0.000
## sample_18     3 0.000
## sample_19     2 0.000
## sample_20     2 0.000
## sample_21     2 0.000
## sample_22     1 0.032
## sample_23     3 0.000
## sample_24     2 0.000
## sample_25     1 0.837
## sample_26     2 0.001
## sample_27     3 0.000
## sample_34     1 0.000
## sample_35     1 0.000
## sample_36     1 0.000
## sample_37     1 0.000
## sample_38     1 0.000
## sample_28     1 0.000
## sample_29     2 0.000
## sample_30     1 0.000
## sample_31     1 0.000
## sample_32     1 0.000
## sample_33     1 0.000
```

We compare to the original classification:


```r
data.frame(cola_class = get_classes(res, k = 3)[, "class"],
           predicted = cl[, "class"])
```

```
##    cola_class predicted
## 1           3         3
## 2           3         3
## 3           2         2
## 4           2         2
## 5           2         2
## 6           3         3
## 7           2         2
## 8           2         2
## 9           2         2
## 10          2         2
## 11          2         2
## 12          2         2
## 13          2         2
## 14          2         2
## 15          2         2
## 16          2         2
## 17          2         2
## 18          3         3
## 19          3         3
## 20          2         2
## 21          1         1
## 22          1         1
## 23          1         1
## 24          1         1
## 25          2         2
## 26          1         1
## 27          1         1
## 28          1         1
## 29          1         1
## 30          1         1
## 31          2         2
## 32          1         1
## 33          1         1
## 34          1         1
## 35          3         3
## 36          2         2
## 37          3         3
## 38          3         3
## 39          2         2
## 40          3         3
## 41          3         3
## 42          3         3
## 43          2         2
## 44          3         3
## 45          2         2
## 46          1         1
## 47          2         2
## 48          2         2
## 49          2         2
## 50          2         2
## 51          2         2
## 52          3         3
## 53          2         2
## 54          2         2
## 55          2         2
## 56          1         1
## 57          3         3
## 58          2         2
## 59          1         1
## 60          2         2
## 61          3         3
## 62          1         1
## 63          1         1
## 64          1         1
## 65          1         1
## 66          1         1
## 67          1         1
## 68          2         2
## 69          1         1
## 70          1         1
## 71          1         1
## 72          1         1
```

`predict_classes()` generates a plot which shows the process of the prediction. 
The left heatmap corresponds to the new matrix and the right small heatmap corresponds
to the signature centroid matrix. The purple annotation on top of the first heatmap
illustrates the distance of each sample to the k signatures.


And if we change to correlation method:


```r
cl = predict_classes(res, k = 3, mat2, dist_method = "correlation")
```

```
## The matrix has been scaled in cola analysis, thus the new matrix should also be scaled
## with the same method ('z-score'). Please double check.
## Set `help = FALSE` to suppress this message. 
## 
## * 70/72 samples (in 3 classes) remain after filtering by silhouette (>= 0.5).
## * cache hash: 5e81712ade334cd09013d9ac870cab7e (seed 888).
## * calculating row difference between subgroups by Ftest.
##   - row difference is extracted from cache.
## * use k-means partition that are already calculated in previous runs.
## * 2076 signatures (50.4%) under fdr < 0.05, group_diff > 0.
## * Predict classes based on 3-group classification (correlation method) on a 72-column matrix.
```

<img src="figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="600" style="display: block; margin: auto;" />

```r
cl
```

```
##           class             p
## sample_39     3  1.879136e-72
## sample_40     3  9.818580e-50
## sample_42     2  4.726915e-28
## sample_47     2  4.368519e-52
## sample_48     2 9.595020e-164
## sample_49     3 8.966020e-145
## sample_41     2 1.355841e-155
## sample_43     2  2.124976e-28
## sample_44     2  1.129990e-63
## sample_45     2  2.862219e-43
## sample_46     2  1.463918e-39
## sample_70     2  9.555535e-09
## sample_71     2  1.883125e-08
## sample_72     2  1.890317e-21
## sample_68     2 7.434706e-169
## sample_69     2 8.727398e-164
## sample_67     2  2.342032e-26
## sample_55     3  1.167792e-50
## sample_56     3  2.129338e-84
## sample_59     2  8.745044e-09
## sample_52     1  4.990070e-35
## sample_53     1  2.225819e-79
## sample_51     1 3.855936e-128
## sample_50     1 2.690689e-133
## sample_54     2  2.396841e-06
## sample_57     1  9.747602e-38
## sample_58     1  8.959459e-91
## sample_60     1  3.846989e-04
## sample_61     1  6.257010e-47
## sample_65     1  6.262522e-59
## sample_66     2  1.996473e-58
## sample_63     1  4.024436e-66
## sample_64     1  3.250437e-59
## sample_62     1  6.065638e-47
## sample_1      3  1.265338e-45
## sample_2      2  1.177492e-23
## sample_3      3  1.683713e-65
## sample_4      3  1.531594e-54
## sample_5      2 6.191909e-170
## sample_6      3  4.257850e-73
## sample_7      3  7.828966e-99
## sample_8      3  4.445058e-97
## sample_9      2  4.291423e-26
## sample_10     3  1.040207e-13
## sample_11     2  1.188730e-33
## sample_12     1  6.107449e-21
## sample_13     2 3.186049e-177
## sample_14     2  4.878387e-33
## sample_15     2 2.349793e-198
## sample_16     2  1.502211e-36
## sample_17     2  1.656820e-16
## sample_18     3  2.807438e-17
## sample_19     2  9.592796e-40
## sample_20     2 2.289534e-176
## sample_21     2  1.806297e-81
## sample_22     1  3.539478e-18
## sample_23     3  4.171210e-76
## sample_24     2  2.445660e-76
## sample_25     1  1.157628e-01
## sample_26     2  1.750686e-07
## sample_27     3 1.362532e-135
## sample_34     1 4.331738e-101
## sample_35     1  6.367166e-78
## sample_36     1 8.758828e-114
## sample_37     1  6.245463e-98
## sample_38     1  5.815944e-57
## sample_28     1  1.616753e-51
## sample_29     2  2.970516e-19
## sample_30     1  2.921410e-53
## sample_31     1  1.693478e-77
## sample_32     1  3.601003e-67
## sample_33     1 5.194996e-163
```

As we can see from the above two heatmaps, correlation method is less strict than
the Euclidean method that the two samples that cannot be assigned to any class with
Euclidean method are assigned with certain classes under correlation method.

`predict_classes()` can also be directly applied to a signature centroid matrix.
Following is how we manually generate the signature centroid matrix for 3-group
classification from Golub dataset:


```r
tb = get_signatures(res, k = 3, plot = FALSE)
```

```
## * 70/72 samples (in 3 classes) remain after filtering by silhouette (>= 0.5).
## * cache hash: 5e81712ade334cd09013d9ac870cab7e (seed 888).
## * calculating row difference between subgroups by Ftest.
##   - row difference is extracted from cache.
## * use k-means partition that are already calculated in previous runs.
## * 2076 signatures (50.4%) under fdr < 0.05, group_diff > 0.
```

```r
# the centroids are already in `tb`, both scaled and unscaled, we just simply extract it
sig_mat = tb[, grepl("scaled_mean", colnames(tb))]
sig_mat = as.matrix(sig_mat)
colnames(sig_mat) = paste0("class", seq_len(ncol(sig_mat)))
head(sig_mat)
```

```
##           class1     class2     class3
## [1,]  0.52395127 -0.5519215  0.3740407
## [2,]  0.06789322 -0.4025209  0.7546084
## [3,]  0.58853214 -0.3077090 -0.2459700
## [4,]  0.47963596 -0.1548451 -0.4051056
## [5,] -0.15871813 -0.2048139  0.6803041
## [6,]  0.40763038 -0.3427458  0.1061578
```

And `sig_mat` can be used in `predict_classes()`:


```r
cl = predict_classes(sig_mat, mat2)
cl = predict_classes(sig_mat, mat2, dist_method = "correlation")
```



```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## attached base packages:
##  [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods  
## [10] base     
## 
## other attached packages:
##  [1] hu6800.db_3.13.0     org.Hs.eg.db_3.13.0  AnnotationDbi_1.54.1 IRanges_2.26.0      
##  [5] S4Vectors_0.30.0     Biobase_2.52.0       BiocGenerics_0.38.0  GetoptLong_1.0.5    
##  [9] mvtnorm_1.1-2        matrixStats_0.59.0   circlize_0.4.13      ComplexHeatmap_2.8.0
## [13] cola_1.8.0           markdown_1.1         knitr_1.33           BiocManager_1.30.16 
## [17] colorout_1.2-2      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_2.0-2       rjson_0.2.20           ellipsis_0.3.2         mclust_5.4.7          
##  [5] XVector_0.32.0         GlobalOptions_0.1.2    clue_0.3-59            bit64_4.0.5           
##  [9] RSpectra_0.16-0        fansi_0.5.0            xml2_1.3.2             codetools_0.2-18      
## [13] splines_4.1.0          doParallel_1.0.16      cachem_1.0.5           impute_1.66.0         
## [17] polyclip_1.10-0        jsonlite_1.7.2         Cairo_1.5-12.2         umap_0.2.7.0          
## [21] annotate_1.70.0        cluster_2.1.2          png_0.1-7              data.tree_1.0.0       
## [25] compiler_4.1.0         httr_1.4.2             assertthat_0.2.1       Matrix_1.3-4          
## [29] fastmap_1.1.0          tools_4.1.0            gtable_0.3.0           glue_1.4.2            
## [33] GenomeInfoDbData_1.2.6 dplyr_1.0.7            Rcpp_1.0.6             slam_0.1-48           
## [37] eulerr_6.1.0           vctrs_0.3.8            Biostrings_2.60.1      iterators_1.0.13      
## [41] polylabelr_0.2.0       xfun_0.24              stringr_1.4.0          mime_0.11             
## [45] lifecycle_1.0.0        irlba_2.3.3            XML_3.99-0.6           dendextend_1.15.1     
## [49] zlibbioc_1.38.0        scales_1.1.1           microbenchmark_1.4-7   RColorBrewer_1.1-2    
## [53] gridExtra_2.3          memoise_2.0.0          reticulate_1.20        ggplot2_3.3.5         
## [57] stringi_1.6.2          RSQLite_2.2.7          highr_0.9              genefilter_1.74.0     
## [61] foreach_1.5.1          shape_1.4.6            GenomeInfoDb_1.28.0    rlang_0.4.11          
## [65] pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.14          lattice_0.20-44       
## [69] purrr_0.3.4            bit_4.0.4              tidyselect_1.1.1       magrittr_2.0.1        
## [73] R6_2.5.0               magick_2.7.2           generics_0.1.0         DBI_1.1.1             
## [77] pillar_1.6.1           survival_3.2-11        KEGGREST_1.32.0        RCurl_1.98-1.3        
## [81] tibble_3.1.2           crayon_1.4.1           utf8_1.2.1             viridis_0.6.1         
## [85] skmeans_0.2-13         blob_1.2.1             digest_0.6.27          xtable_1.8-4          
## [89] brew_1.0-6             openssl_1.4.4          munsell_0.5.0          viridisLite_0.4.0     
## [93] askpass_1.1
```

