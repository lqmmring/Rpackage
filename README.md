# scESI

scESI performs data imputation for single cell RNA sequencing data.
We propose an evolutionary sparse imputation algorithm for single-cell transcriptomes, which constructs a sparse representation model based on gene regulation relationships between cells.

This framework takes into account the topological relationship between cells and the variety of gene expression to iteratively search the global optimal solution, thereby learning the Pareto optimal cell-cell affinity matrix. Finally, we use the learned sparse relationship model between cells to improve data quality and reduce data noise. 

Developer
------------
Qiaoming Liu (cslqm@hit.edu.cn)

Installation
----------------------
Download scESI_0.1.0.tar.gz
```R
install.packages("scESI_0.1.0.tar.gz",repos=NULL,type="source")
```
or install through GitHub
```R
library(devtools)
install_github("lqmmring/Rpackage/scESI")
```


Usage
----------------------

```R
library(scESI)
sce<-readRDS(file = "/data/covid-gene-0.001-0.999-cell-0.003-0.997.rds")

imputation.sce<-sparse_imputation_with_selected_genes(data = sce@assays[["RNA"]]@counts,
                                                      processing = TRUE,
                                                      num.pop=20,
                                                      num.Iteration=30,
                                                      crossover.p=0.7,
                                                      set_num_list=c(50,250,500,800,1000),
                                                      paralle = TRUE,
                                                      cores = 10)
head(imputation.sce[["predictCount"]])
```
An example to show how SCMarker [improve identification of NK cell in GBM data.](https://github.com/KChen-lab/SCMarker/blob/master/test/NK%20cell%20identification%20from%20GBM%20data.pdf)

Acknowledge
-----------------------
The authors would like to appreciate the support and guidance from Dr. G.H. Wang and Dr. J. Li.
