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
install_github("lqmmring/scESI/scESI")
```


Usage
----------------------

```R
library(SCMarker)
data(melanoma)
melanoma1=as.matrix(melanoma[,2:dim(melanoma)[2]])
row.names(melanoma1)=melanoma[,1]
res=ModalFilter(data=melanoma1,geneK=10,cellK=10,width=2)# default width = 1 for UMI data, width =2 for TPM data.
res=GeneFilter(obj=res)
res=getMarker(obj=res,k=300,n=30)
head(res$marker)

```

An example to show how SCMarker [improve identification of NK cell in GBM data.](https://github.com/KChen-lab/SCMarker/blob/master/test/NK%20cell%20identification%20from%20GBM%20data.pdf)

Publication
-----------------------
Wang, Fang, et al. "SCMarker: ab initio marker selection for single cell transcriptome profiling." PLoS computational biology 15.10 (2019): e1007445.
