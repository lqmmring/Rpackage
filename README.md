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
library(scRNAseq)

sce<-scRNAseq::UsoskinBrainData(ensembl = FALSE, location = TRUE)
imputation.sce<-scESI::sparse_imputation_with_selected_genes(data = sce@assays[["RNA"]]@counts,
                                                      processing = TRUE,
                                                      num.pop=20,
                                                      num.Iteration=30,
                                                      crossover.p=0.7,
                                                      set_num_list=c(25,50,80,100,150)
                                                      )
head(imputation.sce[["predictCount"]])
```
An example to show how scESI [discover new cell types on CBMCs from newborns of mothers infected with SARS-CoV-2.](https://github.com/lqmmring/Rpackage)

Acknowledge
-----------------------
The authors would like to appreciate the support and guidance from Dr. G.H. Wang (ghWang@hit.edu.cn)
and Dr. J. Li.
