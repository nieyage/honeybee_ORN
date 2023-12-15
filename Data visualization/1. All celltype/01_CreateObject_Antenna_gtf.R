## create Seurat by gex matrix 
library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(DropletUtils)
library(Matrix)
library(DoubletFinder)
set.seed(1234)

######################################################################################
# 1_FILTERING
# Evaluate min.cells argument, create seurat object, normalize, cluster
######################################################################################

#....................................................................................
# EVALUATE MIN.CELLS ARGUMENT (CELLS PER GENE)
#....................................................................................

  # Create seurat object with no lowerbound on min.cells argument
  dataset.name="NE"
  sratDecontx  <- 
    Read10X("/md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/decontX_outs/NE_decontXcounts") %>%
    CreateSeuratObject(project = dataset.name, min.cells = 0, min.features = 200)
  #Create UMI count matrix and cell counts
  countMtx <- sratDecontx@assays[["RNA"]]@counts
  nonZeroCellNDf = data.frame(nonZeroCellN=rowSums(countMtx!=0))
  #PLOT: Cells per gene histogram (log)
  p1=nonZeroCellNDf %>%
    ggplot(aes(log(nonZeroCellN))) + geom_histogram(binwidth = 0.1)+
    xlab('Number of cells per gene (log)') + 
    ylab('Genes')+
    geom_vline(xintercept = c(log(12)), color='red')+
    ggtitle(paste(dataset.name,": Cells per Gene, vline = 12", sep=""))
  filename <- paste0("./01_QC/1a_filtering_CellsPerGene_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000, p1)

  dataset.name="Nurse"
  sratDecontx  <- 
    Read10X("/md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/decontX_outs/Nurse_decontXcounts") %>%
    CreateSeuratObject(project = dataset.name, min.cells = 0, min.features = 200)
  #Create UMI count matrix and cell counts
  countMtx <- sratDecontx@assays[["RNA"]]@counts
  nonZeroCellNDf = data.frame(nonZeroCellN=rowSums(countMtx!=0))
  #PLOT: Cells per gene histogram (log)
  p1=nonZeroCellNDf %>%
    ggplot(aes(log(nonZeroCellN))) + geom_histogram(binwidth = 0.1)+
    xlab('Number of cells per gene (log)') + 
    ylab('Genes')+
    geom_vline(xintercept = c(log(12)), color='red')+
    ggtitle(paste(dataset.name,": Cells per Gene, vline = 12", sep=""))
  filename <- paste0("./01_QC/1a_filtering_CellsPerGene_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000, p1)

  dataset.name="Forager"
  sratDecontx  <- 
    Read10X("/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/decontX_outs/F_decontXcounts") %>%
    CreateSeuratObject(project = dataset.name, min.cells = 0, min.features = 200)
  #Create UMI count matrix and cell counts
  countMtx <- sratDecontx@assays[["RNA"]]@counts
  nonZeroCellNDf = data.frame(nonZeroCellN=rowSums(countMtx!=0))
  #PLOT: Cells per gene histogram (log)
  p1=nonZeroCellNDf %>%
    ggplot(aes(log(nonZeroCellN))) + geom_histogram(binwidth = 0.1)+
    xlab('Number of cells per gene (log)') + 
    ylab('Genes')+
    geom_vline(xintercept = c(log(12)), color='red')+
    ggtitle(paste(dataset.name,": Cells per Gene, vline = 12", sep=""))
  filename <- paste0("./01_QC/1a_filtering_CellsPerGene_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000, p1)


#....................................................................................
# CREATE SEURAT OBJECT
#....................................................................................

gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']

# file path 
Forager_counts <- Read10X_h5("/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/filtered_feature_bc_matrix.h5")
Forager_fragpath <- "/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/fragments.tsv.gz"
NE_counts <- Read10X_h5("/md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/filtered_feature_bc_matrix.h5")
NE_fragpath <- "/md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/fragments.tsv.gz"
Nurse_counts <- Read10X_h5("/md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/filtered_feature_bc_matrix.h5")
Nurse_fragpath <- "/md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/fragments.tsv.gz"

# create a Seurat object containing the RNA adata;
F_decontXcounts<-Read10X("/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/decontX_outs/F_decontXcounts")
F_Peaks<-Forager_counts$Peaks;
F_barcode<-intersect(colnames(F_decontXcounts),colnames(F_Peaks))
F_decontXcounts <- F_decontXcounts[,F_barcode]
Forager  <- F_decontXcounts %>% CreateSeuratObject(project = "Forager",min.cells = 12, min.features = 200)
F_barcode<-intersect(colnames(Forager),colnames(F_Peaks))
F_Peaks<-F_Peaks[,F_barcode];
# create ATAC assay and add it to the object
Forager[["ATAC"]] <- CreateChromatinAssay(
  counts = F_Peaks,
  sep = c(":", "-"),
  fragments = Forager_fragpath,
  annotation = gene.coords
)
# Nurse 
Nurse_decontXcounts<-Read10X("/md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/decontX_outs/Nurse_decontXcounts")
Nurse_Peaks<-Nurse_counts$Peaks;
Nurse_barcode<-intersect(colnames(Nurse_decontXcounts),colnames(Nurse_Peaks))
Nurse_decontXcounts<-Nurse_decontXcounts[,Nurse_barcode];
Nurse <- Nurse_decontXcounts %>% CreateSeuratObject(project = "Nurse",min.cells = 12, min.features = 200)
Nurse_barcode<-intersect(colnames(Nurse),colnames(Nurse_Peaks))
Nurse_Peaks<-Nurse_Peaks[,Nurse_barcode];
Nurse[["ATAC"]] <- CreateChromatinAssay(
  counts = Nurse_Peaks,
  sep = c(":", "-"),
  fragments = Nurse_fragpath,
  annotation = gene.coords
)
# NE 
NE_decontXcounts<-Read10X("/md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/decontX_outs/NE_decontXcounts")
NE_Peaks<-NE_counts$Peaks;
NE_barcode<-intersect(colnames(NE_decontXcounts),colnames(NE_Peaks))
NE_decontXcounts<-NE_decontXcounts[,NE_barcode];
NE <- NE_decontXcounts %>% CreateSeuratObject(project = "NE",min.cells = 12, min.features = 200)
NE_barcode<-intersect(colnames(NE),colnames(NE_Peaks))
NE_Peaks<-NE_Peaks[,NE_barcode];
NE[["ATAC"]] <- CreateChromatinAssay(
  counts = NE_Peaks,
  sep = c(":", "-"),
  fragments = NE_fragpath,
  annotation = gene.coords
)

# create a new assay to store without decontX information
NE_data<-read.table("/md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/counts_no_double_umi_001.tsv.gz")
NE_data = spread(NE_data, V3, V2)
NE_data[is.na(NE_data)] <- 0
rownames(NE_data)<-NE_data$V1
NE_data<-NE_data[,-1]
NE_barcode<-paste(colnames(NE_data),"-1",sep="")
colnames(NE_data)<-NE_barcode
NE_data<-as(as.matrix(NE_data),"dgCMatrix")
NE_barcode<-intersect(colnames(NE_decontXcounts),colnames(NE_Peaks))
NE_data<-NE_data[,NE_barcode]
raw_assay <- CreateAssayObject(counts = NE_data)
# add this assay to the previously created Seurat object
NE[["raw_RNA"]] <- raw_assay

Nurse_data<-read.table("/md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/counts_no_double_umi_001.tsv.gz")
Nurse_data = spread(Nurse_data, V3, V2)
Nurse_data[is.na(Nurse_data)] <- 0
rownames(Nurse_data)<-Nurse_data$V1
Nurse_data<-Nurse_data[,-1]
Nurse_barcode<-paste(colnames(Nurse_data),"-1",sep="")
colnames(Nurse_data)<-Nurse_barcode
Nurse_data<-as(as.matrix(Nurse_data),"dgCMatrix")
Nurse_barcode<-intersect(colnames(Nurse_decontXcounts),colnames(Nurse_Peaks))
Nurse_data<-Nurse_data[,Nurse_barcode]
raw_assay <- CreateAssayObject(counts = Nurse_data)
# add this assay to the previously created Seurat object
Nurse[["raw_RNA"]] <- raw_assay
#Forager
Forager_data<-read.table("/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/counts_no_double_umi_001.tsv.gz")
Forager_data = spread(Forager_data, V3, V2)
Forager_data[is.na(Forager_data)] <- 0
rownames(Forager_data)<-Forager_data$V1
Forager_data<-Forager_data[,-1]
Forager_barcode<-paste(colnames(Forager_data),"-1",sep="")
colnames(Forager_data)<-Forager_barcode
Forager_data<-as(as.matrix(Forager_data),"dgCMatrix")
Forager_barcode<-intersect(colnames(F_decontXcounts),colnames(F_Peaks))
Forager_data<-Forager_data[,Forager_barcode]
raw_assay <- CreateAssayObject(counts = Forager_data)
# add this assay to the previously created Seurat object
Forager[["raw_RNA"]] <- raw_assay

# get MT gtf 
cat GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf | awk '{if($1 == "MT"){print $0}}' > GCF_003254395.2_Amel_HAv3.1_genomic_MitoFeatures.gtf
MT_gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_MitoFeatures.gtf')
# calculate the percentage of mtRNA 
objList <- list(NE,Nurse,Forager)
for (i in seq_len(length(objList))) {
  DefaultAssay(objList[[i]]) <- "RNA"
  MT_feature <- unique(MT_gtf$gene_name[which(MT_gtf$gene_name %in% rownames(objList[[i]]))])
  objList[[i]][["percent.mt"]] = colSums(x = GetAssayData(object = objList[[i]], 
  assay = "RNA", slot = "counts")[MT_feature,,drop = FALSE])/objList[[i]][[paste0("nCount_","RNA")]] * 100
    }

######################################################################################
# 2_DOUBLETFINDER
######################################################################################

# NORMALIZE for DoubletFinder
library(gridExtra)
sample <- c("NE","Nurse","Forager")
for (i in seq_len(length(objList))) {
  # Normalizing
  dataset.name <- sample[i]
  sratDecontx    <- SCTransform(objList[[i]], verbose = F)
  # Run cluster analysis
  sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
  sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindClusters(sratDecontx, verbose = T)
  # PLOT: umap/tsne
  p1 <- DimPlot(sratDecontx, reduction = "umap", label = TRUE)
  p2 <- DimPlot(sratDecontx, reduction = "tsne", label = TRUE)
  g = arrangeGrob(p1,p2, ncol = 2)
  filename <- paste0("./01_QC/1b_filtering_UMAP-tSNE_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(g)
  # PLOT: coreceptor tsne
  p4 <- FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Or2'), ncol = 1)
  filename <- paste0("./01_QC/1c_filtering_CoReceptor-tSNE_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p4)
  objList[[i]] <- sratDecontx
    }

objList2<-list()
for (i in seq_len(length(objList))) {
  dataset.name<-sample[i]
  print(dataset.name)
  sratDecontx <- objList[[i]]
  # Compute expected doublet rate
  cellN=nrow(sratDecontx@meta.data)
  expDoubletRate = (cellN*0.0008 + 0.0527)*0.01
  normalizationMethod='SCTransform'
  sweep.res.list_scData <- paramSweep_v3(sratDecontx, 
                                         PCs = 1:50, 
                                         sct = normalizationMethod == 'SCTransform', 
                                         num.cores = 4) #num.cores = 4
  sweep.stats_scData <- summarizeSweep(sweep.res.list_scData, GT = FALSE)
  bcmvn_scData <- find.pK(sweep.stats_scData)
  bcmvn_scData$pK <- as.numeric(as.character(bcmvn_scData$pK))
  pK1=bcmvn_scData$pK[bcmvn_scData$BCmetric==max(bcmvn_scData$BCmetric)]
  print(head(pK1))
  # PLOT: pK selection
  p1=ggplot(data=bcmvn_scData, 
            aes(x=pK, y=BCmetric, group=2)) +
    geom_line(color="blue")+
    geom_point()+
    geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
    labs(title="pK Selection",x="pK", y = "BCmvn")+
    theme_classic()
  filename <- paste0("./01_QC/2a_doubletfinder_pkselection_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)
  # More doublet finder
  pK1=as.numeric(as.character( pK1 ))
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- sratDecontx@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)   
  nExp_poi <- round(expDoubletRate*nrow(sratDecontx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sratDecontx <- doubletFinder_v3( sratDecontx, PCs = sratDecontx@commands$RunUMAP.SCT.pca$dims,
                                   pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  sratDecontx@meta.data$DoubletFinder =  sratDecontx@meta.data[,grep('DF.classifications', colnames( sratDecontx@meta.data))]
  # PLOT: Doublet Finder graphs
  p2 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'DoubletFinder')
  p3 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'DoubletFinder')
  g = arrangeGrob(p2,p3, ncol = 2)
  filename <- paste0("./01_QC/2b_doubletfinder_ScatterPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(filename, g)
  # PLOT: Violin Plots
  p4 <- VlnPlot(sratDecontx, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'DoubletFinder', pt.size = 0)
  filename <- paste0("./01_QC/2c_doubletfinder_DoubletFinder-ViolinPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p4)
  objList2[[i]] <- sratDecontx
}


######################################################################################
# 3_FEATURE FILTERING
######################################################################################

objList<- objList2
# Subset identified singlets
for (i in seq_len(length(objList))) {
    objList[[i]] <- subset(objList[[i]], subset = DoubletFinder == "Singlet")
}

# calculate the score of NS and TSS
for (i in seq_len(length(objList))) {
  DefaultAssay(objList[[i]]) <- "ATAC"
    objList[[i]] <- NucleosomeSignal(objList[[i]])
    objList[[i]] <- TSSEnrichment(objList[[i]],fast=FALSE)
    }
library(ggplot2)


# plot TSS and fragrament distribution plot 
# NE
  pdf("./01_QC/3a_TSS_distribution_NE.pdf")
  objList[[1]]$high.tss<-ifelse(objList[[1]]$TSS.enrichment > 1, 'High', 'Low')
  TSS<-TSSPlot(objList[[1]], group.by = 'high.tss') + NoLegend()+ labs(title = "NE")
  objList[[1]]$nucleosome_group <- ifelse(objList[[1]]$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
  #Frag<-FragmentHistogram(object = objList[[1]], group.by = 'nucleosome_group')+ labs(title = "NE")
  print(TSS);#print(Frag);
  dev.off();
# Nurse
  pdf("./01_QC/3a_TSS_distribution_Nurse.pdf")
  objList[[2]]$high.tss<-ifelse(objList[[2]]$TSS.enrichment > 1, 'High', 'Low')
  TSS<-TSSPlot(objList[[2]], group.by = 'high.tss') + NoLegend()+ labs(title = "Nurse")
  objList[[2]]$nucleosome_group <- ifelse(objList[[2]]$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
  #Frag<-FragmentHistogram(object = objList[[2]], group.by = 'nucleosome_group')+ labs(title = "Nurse")
  print(TSS);#print(Frag);
  dev.off();
# Forager
  pdf("./01_QC/3a_TSS_distribution_Forager.pdf")
  objList[[3]]$high.tss<-ifelse(objList[[3]]$TSS.enrichment > 1, 'High', 'Low')
  TSS<-TSSPlot(objList[[3]], group.by = 'high.tss') + NoLegend()+ labs(title = "Forager")
  objList[[3]]$nucleosome_group <- ifelse(objList[[3]]$nucleosome_signal > 2.5, 'NS > 2', 'NS < 2')
  #Frag<-FragmentHistogram(object = objList[[3]], group.by = 'nucleosome_group')+ labs(title = "Forager")
  print(TSS);#print(Frag);
  dev.off();


for (i in seq_len(length(objList))) {
  # plot QC plot 
  pdf(file = paste("./01_QC/3b_QC_before_",sample[i],".pdf", sep = ""),width=10,height=6)
  qc<-VlnPlot(object = objList[[i]],
            features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","percent.mt","nFeature_RNA"),
            ncol = 3,
            pt.size = 0.01
          )
  densityRNA<-plot(density(objList[[i]]@meta.data$nCount_RNA),xlim=c(0,5000))
  densityRNA_feature<-plot(density(objList[[i]]@meta.data$nFeature_RNA))
  densityATAC<-plot(density(objList[[i]]@meta.data$nCount_ATAC),xlim=c(0,10000))
  print(qc)
  print(densityRNA)
  print(densityRNA_feature)
  print(densityATAC)
  dev.off();
   }

# filter out low quality cells
# To remove doublets,select different cutoff#####
# NE 
  objList2<-c()
  objList2[[1]]<-subset(x=objList[[1]],
   subset = nCount_ATAC < 100000 & nCount_ATAC > 2000 &
    nCount_RNA < 30000 & nCount_RNA > 300 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 5
    )
  print(objList2[[1]])
# Nurse 
  objList2[[2]]<-subset(x=objList[[2]],
   subset = nCount_ATAC < 100000 &nCount_ATAC > 1900 &
    nCount_RNA < 30000 & nCount_RNA > 200 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1&
    percent.mt < 5
    )
  print(objList2[[2]])
# Forager 
  objList2[[3]]<-subset(x=objList[[3]],
   subset = nCount_ATAC < 100000 &
    nCount_RNA < 30000 &
    nCount_ATAC > 300 &
    nCount_RNA > 300 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1&
    percent.mt < 5
    )
  print(objList2[[3]])
objList<-objList2

  # Calculate feature filtering parameters

    objList2<-c()
for (i in seq_len(length(objList))) {
  srat <- objList[[i]]
  medianFeature<- median(srat$nFeature_RNA)
  madFeature <- mad(srat$nFeature_RNA)
  lowerBound<- 250
  upperBound <- medianFeature+madFeature*3
  # Filter seurat object
  srat.filter=subset(srat, subset = (nFeature_RNA < medianFeature+madFeature*3) & 
                       (nFeature_RNA>= lowerBound) & (percent.mt < 5))

  # PLOT: Violin Plot Before/After Filtering
  p1 = VlnPlot(srat, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","percent.mt","nFeature_RNA"), ncol = 6, 
               group.by = 'orig.ident', pt.size = 0.1)
  filename <- paste0("./01_QC/3c_filtering_ViolinPlot_before_",srat@project.name, ".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width =6000, height = 2000,p1)
  
  p2 = VlnPlot(srat.filter, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","percent.mt","nFeature_RNA"), ncol = 6, 
               group.by = 'orig.ident', pt.size = 0.1)
  filename <- paste0("./01_QC/3d_filtering_ViolinPlot_after_",srat@project.name, ".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 2000,p2)
  # PLOT: Cells per gene
  p3=ggplot(srat@meta.data, aes(nFeature_RNA)) + geom_histogram(binwidth = 50) +
    geom_vline(xintercept = c(medianFeature), color='black') + coord_cartesian(xlim = c(0,4000)) +
    geom_vline(xintercept = c(lowerBound, upperBound), color='red')+
    ggtitle(paste0(srat@project.name,": median=",medianFeature,
                   ", lowerBound=",lowerBound,", upperBound=",upperBound,
                   ", Cells pre-filter=",nrow(srat@meta.data),
                   ", Cells post-filter=",nrow(srat.filter@meta.data)))
  filename <- paste0("./01_QC/3e_filtering_FeatureHistogram_",srat@project.name, ".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p3)
  objList2[[i]] <-  srat.filter
}

objList <- objList2

# Peak Calling 
# quantify counts in each peak
fragpath<-list(NE_fragpath,Nurse_fragpath,Forager_fragpath)
for(i in seq_len(length(objList))){
	# call peaks using MACS2
    peaks <- CallPeaks(objList[[i]], macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2")
    macs2_counts <- FeatureMatrix(
      fragments = Fragments(objList[[i]]),
      features = peaks,
      cells = colnames(objList[[i]])
    )     
    # create a new assay using the MACS2 peak set and add it to the Seurat object
    objList[[i]][["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = fragpath[[i]],
      annotation = gene.coords
    )
}

#######integrate RNA and ATAC#####################
# Simply merge Seurat objects
merged_obj <- merge(x=objList[[1]],y=c(objList[[2]],objList[[3]]),add.cell.ids = c("NE","Nurse","Forager"),
	project = "honeybee")
Idents(merged_obj) <- gsub("_.*", "", colnames(merged_obj))
merged_obj$orig.ident<-Idents(merged_obj)
honeybee.list <- SplitObject(merged_obj, split.by = "orig.ident")

for (i in 1:length(honeybee.list)) {
	DefaultAssay(honeybee.list[[i]]) <- "RNA"
  honeybee.list[[i]] <- SCTransform(honeybee.list[[i]], verbose = FALSE)
}
for (i in seq_len(length(honeybee.list))) {
  DefaultAssay(honeybee.list[[i]]) <- "SCT"
}
honeybee.features <- SelectIntegrationFeatures(object.list = honeybee.list, nfeatures = 3000)
honeybee.list <- PrepSCTIntegration(object.list = honeybee.list, anchor.features = honeybee.features)
#integrate RNA using rpca

honeybee_list <- lapply(
  X = honeybee.list,
  FUN = RunPCA,
  features = honeybee.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list = honeybee_list,
  normalization.method = "SCT",
  anchor.features = honeybee.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)

honeybee <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA",
  dims = 1:30
)

#run LSI on new seurat object with integrated RNA assay
DefaultAssay(honeybee) <- "ATAC"
honeybee <- RunTFIDF(honeybee)
honeybee <- FindTopFeatures(honeybee, min.cutoff = "q25")
honeybee <- RunSVD(honeybee)

pdf("./02_All_celltype/LSI-depth-correalation.pdf")
# assess the correlation between each LSI component and sequencing depth using the DepthCor() function
DepthCor(honeybee)
dev.off();

#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
honeybee_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI",
  reductions = honeybee@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
honeybee@reductions$integratedLSI <- honeybee_atac@reductions$integratedLSI

#####done integrate ATAC and RNA ################
# RNA analysis
DefaultAssay(honeybee) <- "integratedRNA"
honeybee <- RunPCA(honeybee) %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(honeybee) <- "ATAC"
honeybee <- RunUMAP(honeybee, reduction = 'integratedLSI', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
# We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. 
# We use this graph for UMAP visualization and clustering
#honeybee <- RunPCA(honeybee, npcs = 50, verbose = FALSE)
honeybee <- FindMultiModalNeighbors(
  object = honeybee,
  reduction.list = list("pca", "integratedLSI"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
honeybee <- RunUMAP(honeybee, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
honeybee <- FindClusters(honeybee, graph.name = "wsnn", resolution =1.2, algorithm = 3, verbose = FALSE)

###reorder the level of sample#####
Idents(honeybee)<-honeybee$orig.ident
honeybee$orig.ident<-factor(honeybee$orig.ident,levels=c("NE","Nurse","Forager"))
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./02_All_celltype/honeybee_cluster_WNN.pdf.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(honeybee, cols=my47colors, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(honeybee, cols=my47colors, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(honeybee, cols=my47colors, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
###sample
p1 <- DimPlot(honeybee, cols=my47colors, reduction = "umap.rna", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(honeybee, cols=my47colors, reduction = "umap.atac",group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(honeybee, cols=my47colors, reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))& NoLegend()
dev.off()

DefaultAssay(honeybee) <- "RNA"
Idents(honeybee)<-honeybee$seurat_clusters
pdf('./02_All_celltype/coreceptor_VlnPlot_WNN.pdf',width=15, height=10)
print( VlnPlot(honeybee, features = c("Or2","LOC100577715","LOC100576282","LOC107963999","LOC102656838"), ncol = 1, pt.size = 0) )
dev.off()
pdf('./02_All_celltype/coceptor_FeaturePlot_WNN.pdf', width=12, height=8)
print(FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 10, features = c("Or2","LOC100577715","LOC100576282","LOC107963999","LOC102656838"), ncol = 3))
dev.off()
saveRDS(honeybee,"./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
