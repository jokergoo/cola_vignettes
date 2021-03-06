Compare Two Partitioning Results
=============================================================

**Author**: Zuguang Gu ( z.gu@dkfz.de )

**Date**: 2021-07-14

**Package version**: 1.9.3

-------------------------------------------------------------




**cola** allows to perform multiple partitioning methods in paralle. Different partitioning methods might give different
partitions. Here the function `compare_partitions()` can compare two partitioning results by a
HTML report.

We take `golub_cola` as an example which contains consensus partitioning results of the Golub dataset.


```r
library(cola)
data(golub_cola)
```

Basically we can directly provide two `ConsensusPartition` objects to this function:


```r
compare_partitions(golub_cola["ATC:skmeans"], golub_cola["SD:kmeans"])
```

`compare_partitions()` generates a HTML file, and it automatically opens it in the web browser. If the matrix rows are genes
with IDs as Ensemble IDs, RefSeq IDs or gene symbols, functional enrichment will also be applied to the signature genes under
the two sets of partitions to see which is more biological meaningful.

Genes in `golub_cola` have microarray probe IDs, so here we need to provide a gene ID mapping vector which maps probe IDs
to Entrez IDs, as what is done in the following code:


```r
require(hu6800.db)
x = hu6800ENTREZID
mapped_probes = mappedkeys(x)
id_mapping = unlist(as.list(x[mapped_probes]))
head(id_mapping)
```

```
##     A28102_at   AB000114_at   AB000115_at   AB000220_at AB000381_s_at   AB000409_at 
##        "2556"        "4958"       "10964"       "10512"        "2765"        "8569"
```

Note if the ID mapping is not provided, the functional enrichment will not be applied.


```r
compare_partitions(golub_cola["ATC:skmeans"], golub_cola["SD:kmeans"], 
    id_mapping = id_mapping, output_file = "compare_partitions_example.html")
```

An example output can be found [in this link](compare_partitions_example.html).
