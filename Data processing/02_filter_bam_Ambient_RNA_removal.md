## filter bam and Ambient RNA removal:
### After cellranger 

1. Step1: Modify Bam
```
cd outs
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-antenan/outs/run_dedup.sh .
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-antenan/outs/parse_cellsorted_bam.py .
nohup bash run_dedup.sh gex_possorted_bam &
```

2. Step2: Manage peak bed file
```
gunzip atac_fragments.tsv.gz
sed -e '/#/d;' atac_fragments.tsv >fragments.tsv
cat fragments.tsv |sort -k1,1 -k2,2n | bgzip > fragments.tsv.gz
tabix -b 2 -e 3 -p bed fragments.tsv.gz
```

3. Step3: DecontX in R 

```
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(scater)
  #library(Signac)
  library(celda)
  library(DropletUtils)
  })

theme_set(theme_cowplot())
```
* Forager

```
1. Use SoupX for removing ambient RNAs by modified bam files 
Forager_data<-read.table("./counts_no_double_umi_001.tsv.gz")
Forager_counts = spread(Forager_data, V3, V2)
Forager_counts[is.na(Forager_counts)] <- 0
rownames(Forager_counts)<-Forager_counts$V1
Forager_counts<-Forager_counts[,-1]
Forager_barcode<-paste(colnames(Forager_counts),"-1",sep="")
colnames(Forager_counts)<-Forager_barcode
srat <-  CreateSeuratObject(counts = Forager_counts, min.cells=3,  project="Forager", assay = "RNA")
raw.matrix  <- Read10X_h5("./raw_feature_bc_matrix.h5",use.names = T)
srat.raw <- CreateSeuratObject(
        counts = raw.matrix$`Gene Expression`,
        assay = "RNA"
      )
    srat    <- SCTransform(srat, verbose = F)
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSNE(srat, dims = 1:50, verbose = F)
    srat    <- FindNeighbors(srat,reduction = 'pca', dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T,resolution =1.2)

    pdf('./decontX_outs/Original_UMAP_tSNE_F.pdf',width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsne", label = TRUE))
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_VlnPlot_F.pdf', width=15, height=8)
    print( VlnPlot(srat, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_FeaturePlot_F.pdf', width=16, height=16)
    print(FeaturePlot(srat, reduction = 'tsne', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
   # remove ambient RNAs
   library(Matrix)
   data_full<-as(as.matrix(Forager_counts),"dgCMatrix")
   #data_full<-Matrix(Forager_counts, sparse = TRUE) 
   data_empty<-raw.matrix$`Gene Expression`
   data_empty<-data_empty[rownames(data_full),]
   data_full_sce <- SingleCellExperiment(assays = list(counts = data_full))
   data_empty_sce <- SingleCellExperiment(assays = list(counts = data_empty))
   data_full_sce <- decontX(x = data_full_sce, background = data_empty_sce)
   sce.decontX<-data_full_sce
    # Cluster labels on UMAP
    umap <- reducedDim(sce.decontX, "decontX_UMAP")

    pdf('./decontX_outs/decontX_contamination_F.pdf',width=15, height=6)
    print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
            plotDecontXContamination(sce.decontX))
    dev.off()
    DropletUtils:::write10xCounts("./decontX_outs/F_decontXcounts", round(decontXcounts(sce.decontX)),
                                  barcodes = rownames(colData(data_full_sce)))  
    #########################################
    ### decontX
      sratDecontx  <- 
      Read10X("./decontX_outs/F_decontXcounts") %>%
      CreateSeuratObject(project = "Forager", min.cells = 3, min.features = 200)   
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-") 
    sratDecontx    <- SCTransform(sratDecontx, verbose = F)
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:30, reduction = 'pca',verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T,resolution =1.2)   
    pdf('./decontX_outs/decontX_UMAP_tSNE_F.pdf', width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsne", label = TRUE) )
    dev.off()

    pdf('./decontX_outs/decontX_coreceptor_VlnPlot_F.pdf',width=15, height=8)
    print( VlnPlot(sratDecontx, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/decontX_coreceptor_FeaturePlot_F.pdf', width=16, height=16)
    print(FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
```

* Nurse

``` 
Nurse_data<-read.table("./counts_no_double_umi_001.tsv.gz")
Nurse_counts = spread(Nurse_data, V3, V2)
Nurse_counts[is.na(Nurse_counts)] <- 0
rownames(Nurse_counts)<-Nurse_counts$V1
Nurse_counts<-Nurse_counts[,-1]
Nurse_barcode<-paste(colnames(Nurse_counts),"-1",sep="")
colnames(Nurse_counts)<-Nurse_barcode

srat <-  CreateSeuratObject(counts = Nurse_counts, min.cells=3,  project="Nurse", assay = "RNA")

raw.matrix  <- Read10X_h5("./raw_feature_bc_matrix.h5",use.names = T)

srat.raw <- CreateSeuratObject(
        counts = raw.matrix$`Gene Expression`,
        assay = "RNA"
      )

    srat    <- SCTransform(srat, verbose = F)
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSNE(srat, dims = 1:50, verbose = F)
    srat    <- FindNeighbors(srat,reduction = 'pca', dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T,resolution =1.2)

    pdf('./decontX_outs/Original_UMAP_tSNE_Nurse.pdf',width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsne", label = TRUE))
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_VlnPlot_Nurse.pdf', width=15, height=8)
    print( VlnPlot(srat, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_FeaturePlot_Nurse.pdf', width=16, height=16)
    print(FeaturePlot(srat, reduction = 'tsne', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()

   # remove ambient RNAs
   library(Matrix)
   data_full<-as(as.matrix(Nurse_counts),"dgCMatrix")
   #data_full<-Matrix(Nurse_counts, sparse = TRUE) 
   data_empty<-raw.matrix$`Gene Expression`
   data_empty<-data_empty[rownames(data_full),]
   data_full_sce <- SingleCellExperiment(assays = list(counts = data_full))
   data_empty_sce <- SingleCellExperiment(assays = list(counts = data_empty))
   data_full_sce <- decontX(x = data_full_sce, background = data_empty_sce)
   sce.decontX<-data_full_sce
   # Cluster labels on UMAP
   umap <- reducedDim(sce.decontX, "decontX_UMAP")
   pdf('./decontX_outs/decontX_contamination_Nurse.pdf',width=15, height=6)
   print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
           plotDecontXContamination(sce.decontX))
   dev.off()
   DropletUtils:::write10xCounts("./decontX_outs/Nurse_decontXcounts", round(decontXcounts(sce.decontX)),
                                 barcodes = rownames(colData(data_full_sce)))
       
    #########################################
    ### decontX
    sratDecontx  <- 
    Read10X("./decontX_outs/Nurse_decontXcounts") %>%
    CreateSeuratObject(project = "Nurse", min.cells = 3, min.features = 200)
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
    sratDecontx    <- SCTransform(sratDecontx, verbose = F)
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:30, reduction = 'pca',verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T,resolution =1.2)
    pdf('./decontX_outs/decontX_UMAP_tSNE_Nurse.pdf', width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsne", label = TRUE) )
    dev.off()
    pdf('./decontX_outs/decontX_coreceptor_VlnPlot_Nurse.pdf',width=15, height=8)
    print( VlnPlot(sratDecontx, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    pdf('./decontX_outs/decontX_coreceptor_FeaturePlot_Nurse.pdf', width=16, height=16)
    print(FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
```

* NE

``` 
NE_data<-read.table("./counts_no_double_umi_001.tsv.gz")
NE_counts = spread(NE_data, V3, V2)
NE_counts[is.na(NE_counts)] <- 0
rownames(NE_counts)<-NE_counts$V1
NE_counts<-NE_counts[,-1]
NE_barcode<-paste(colnames(NE_counts),"-1",sep="")
colnames(NE_counts)<-NE_barcode

srat <-  CreateSeuratObject(counts = NE_counts, min.cells=3,  project="NE", assay = "RNA")

raw.matrix  <- Read10X_h5("./raw_feature_bc_matrix.h5",use.names = T)

srat.raw <- CreateSeuratObject(
        counts = raw.matrix$`Gene Expression`,
        assay = "RNA"
      )

    srat    <- SCTransform(srat, verbose = F)
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSNE(srat, dims = 1:50, verbose = F)
    srat    <- FindNeighbors(srat,reduction = 'pca', dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T,resolution =1.2)

    pdf('./decontX_outs/Original_UMAP_tSNE_NE.pdf',width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsne", label = TRUE))
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_VlnPlot_NE.pdf', width=15, height=8)
    print( VlnPlot(srat, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/Original_coreceptor_FeaturePlot_NE.pdf', width=16, height=16)
    print(FeaturePlot(srat, reduction = 'tsne', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()

   # remove ambient RNAs
   library(Matrix)

   data_full<-as(as.matrix(NE_counts),"dgCMatrix")
   #data_full<-Matrix(NE_counts, sparse = TRUE) 

   data_empty<-raw.matrix$`Gene Expression`
   data_empty<-data_empty[rownames(data_full),]
   data_full_sce <- SingleCellExperiment(assays = list(counts = data_full))
   data_empty_sce <- SingleCellExperiment(assays = list(counts = data_empty))
   data_full_sce <- decontX(x = data_full_sce, background = data_empty_sce)
   sce.decontX<-data_full_sce
    # Cluster labels on UMAP
    umap <- reducedDim(sce.decontX, "decontX_UMAP")

    pdf('./decontX_outs/decontX_contamination_NE.pdf',width=15, height=6)
    print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
            plotDecontXContamination(sce.decontX))
    dev.off()

    DropletUtils:::write10xCounts("./decontX_outs/NE_decontXcounts", round(decontXcounts(sce.decontX)),
                                  barcodes = rownames(colData(data_full_sce)))
        
    #########################################
    ### decontX
      sratDecontx  <- 
      Read10X("./decontX_outs/NE_decontXcounts") %>%
      CreateSeuratObject(project = "NE", min.cells = 3, min.features = 200)
    
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
    

    sratDecontx    <- SCTransform(sratDecontx, verbose = F)
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:30, reduction = 'pca',verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T,resolution =1.2)
    
    
    pdf('./decontX_outs/decontX_UMAP_tSNE_NE.pdf', width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsne", label = TRUE) )
    dev.off()
    
    pdf('./decontX_outs/decontX_coreceptor_VlnPlot_NE.pdf',width=15, height=8)
    print( VlnPlot(sratDecontx, features = c('Or2',"Or51","Or52","Or4"), ncol = 1, pt.size = 0) )
    dev.off()
    
    pdf('./decontX_outs/decontX_coreceptor_FeaturePlot_NE.pdf', width=16, height=16)
    print(FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Or2',"Or51","Or52","Or4"), ncol = 2))
    dev.off()
```



