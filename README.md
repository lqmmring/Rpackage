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

sce<-scRNAseq::BacherTCellData(ensembl = FALSE, location = TRUE)
imputation.sce<-scESI::sparse_imputation_with_selected_genes(data = sce@assays@data@listData[["rpm"]],
                                                      processing = TRUE,
                                                      num.pop=20,
                                                      num.Iteration=30,
                                                      crossover.p=0.7,
                                                      set_num_list=c(25,50,80,100,150)
                                                      )
head(row.names(imputation.sce[["predictCount"]]))
head(colnames(imputation.sce[["predictCount"]]))
saveRDS(imputation.sce, file='imputed-results.rds')
```

Clustering with [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and [Clustree](https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html)
----------------------

```R
library(Seurat)

imputation.sce<-readRDS(file='imputed-results.rds')
sce.impute <- CreateSeuratObject(counts = imputation.sce[["predictCount"]], 
                                 project = "imputed-covid")

sce.impute<- FindVariableFeatures(sce.impute, selection.method = 'vst',
                                  nfeatures = 1000)
sce.impute <- ScaleData(sce.impute)
sce.impute <- RunPCA(sce.impute, features = VariableFeatures(object = sce.impute)) 
sce.impute <- FindNeighbors(sce.impute, dims = 1:15)
sce.impute <- FindClusters(sce.impute, resolution = c(seq(.1,1.6,.2)))
sce.impute <- RunUMAP(sce.impute, dims = 1:15)
sce.impute <- RunTSNE(sce.impute, dims = 1:15)
clustree(sce.impute@meta.data, prefix = "RNA_snn_res.")

features = c("MS4A1", "CD79A","CD19", #B cell
             "MZB1", #Plasma cell
             "CD3E", "CD3D",'IL7R','CD4',"CD8A","CCR7", #T cell
             "CD14","LYZ","FCGR3A",'S100A4','S100A9', #Monocyte
             "CST3","CD1C", #Monocyte derived DC
             "IFITM3", "APOBEC3A","SERPINF1","ZEB2","CD1C", #pDC
             "PTPRC","CD68","CD163","C1QA","FPR1","ITGAM", #Macrophage
             "GNLY","NKG7",#NK
             "ELANE","LTF","MMP8",#Neutrophil
             "IL1B","NFKBIA","DUSP2",#Inflammatory cytokines
             "CXCL8"#chemokines
)

FeaturePlot(sce.impute, features = unique(features))
DotPlot(sce.impute, group.by = 'RNA_snn_res.0.9',features = unique(features)) + RotatedAxis()
DoHeatmap(subset(sce.impute),group.by = 'RNA_snn_res.0.9',features = features, size = 3)

saveRDS(sce.impute,file = "imputed-results.rds")
```

With other scRNA-seq data analysis tools, scESI can be extended to [infer single cell trajectories](http://cole-trapnell-lab.github.io/monocle-release/docs/) and [construct cell-cell communication network](https://scenic.aertslab.org/tutorials/).


Acknowledge
-----------------------
The authors would like to appreciate the support and guidance from Dr. G.H. Wang (ghWang@hit.edu.cn)
and Dr. J. Li.
