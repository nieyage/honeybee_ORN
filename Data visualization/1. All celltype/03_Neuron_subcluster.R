## create Seurat by gex matrix 
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(DropletUtils)
library(Matrix)
set.seed(1234)
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
DefaultAssay(honeybee)<-"RNA"
Neuron <- subset(honeybee,idents=c("Neuron"));
Neuron.list <- SplitObject(Neuron, split.by = "orig.ident")

for (i in 1:length(Neuron.list)) {
  DefaultAssay(Neuron.list[[i]]) <- "RNA"
  Neuron.list[[i]] <- SCTransform(Neuron.list[[i]], verbose = FALSE)
}
for (i in seq_len(length(Neuron.list))) {
  DefaultAssay(Neuron.list[[i]]) <- "SCT"
}
Neuron.features <- SelectIntegrationFeatures(object.list = Neuron.list, nfeatures = 3000)
Neuron.list <- PrepSCTIntegration(object.list = Neuron.list, anchor.features = Neuron.features)
#integrate RNA using rpca
Neuron_list <- lapply(
  X = Neuron.list,
  FUN = RunPCA,
  features = Neuron.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list = Neuron_list,
  normalization.method = "SCT",
  anchor.features = Neuron.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)

Neuron <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA",
  dims = 1:30
)

#run LSI on new seurat object with integrated RNA assay
DefaultAssay(Neuron) <- "ATAC"
Neuron <- RunTFIDF(Neuron)
Neuron <- FindTopFeatures(Neuron, min.cutoff = "q25")
Neuron <- RunSVD(Neuron)

pdf("./03_Neuron/LSI-depth-correalation.pdf")
# assess the correlation between each LSI component and sequencing depth using the DepthCor() function
DepthCor(Neuron)
dev.off();

#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
Neuron_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI",
  reductions = Neuron@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
Neuron@reductions$integratedLSI <- Neuron_atac@reductions$integratedLSI

#####done integrate ATAC and RNA ################
# RNA analysis

DefaultAssay(Neuron) <- "integratedRNA"
Neuron <- RunPCA(Neuron) %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(Neuron) <- "ATAC"
Neuron <- RunUMAP(Neuron, reduction = 'integratedLSI', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
# We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. 
# We use this graph for UMAP visualization and clustering
#Neuron <- RunPCA(Neuron, npcs = 50, verbose = FALSE)
Neuron <- FindMultiModalNeighbors(
  object = Neuron,
  reduction.list = list("pca", "integratedLSI"), 
  dims.list = list(1:50, 2:25),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
Neuron <- RunUMAP(Neuron, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Neuron <- FindClusters(Neuron, graph.name = "wsnn", resolution =3, algorithm = 3, verbose = FALSE)
###reorder the level of sample#####
Idents(Neuron)<-Neuron$orig.ident
Neuron$orig.ident<-factor(Neuron$orig.ident,levels=c("NE","Nurse","Forager"))
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')

pdf("./03_Neuron/Neuron_cluster_WNN.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(Neuron, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(Neuron, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(Neuron, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
###sample
p1 <- DimPlot(Neuron, reduction = "umap.rna", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(Neuron, reduction = "umap.atac",group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(Neuron, reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))& NoLegend()
dev.off()

DefaultAssay(Neuron) <- "RNA"
Idents(Neuron)<-Neuron$seurat_clusters
pdf('./03_Neuron/coreceptor_VlnPlot_WNN.pdf',width=15, height=10)
print( VlnPlot(Neuron, features = c("Or2","LOC100577715","LOC100576282","LOC107963999","LOC102656838"), ncol = 1, pt.size = 0) )
dev.off()
pdf('./03_Neuron/coceptor_FeaturePlot_WNN.pdf', width=12, height=8)
print(FeaturePlot(Neuron, reduction = 'wnn.umap',max.cutoff = 10, features = c("Or2","LOC100577715","LOC100576282","LOC107963999","LOC102656838"), ncol = 3))
dev.off()

Annotation<- c(rep("Orco+Neuron",7332))
Annotation[which(Neuron$seurat_clusters%in%c("31"))]<-"Orco-Neuron"
Neuron$Annotation<- Annotation
Idents(Neuron)<- Neuron$Annotation

DefaultAssay(Neuron) <- "RNA";
saveRDS(Neuron,"./03_Neuron/WNN_Neuron_integrated.rds")

# add Neuron info in Allcelltype 
ORN_barcode <-  names(which(Neuron$Annotation=="Orco+Neuron"))
Neuron_barcode<- names(which(Neuron$Annotation=="Orco-Neuron"))

Annotation_subcluster<- as.character(honeybee$Annotation)
Annotation_subcluster[which(colnames(honeybee)%in% ORN_barcode)]="Orco+Neuron";
Annotation_subcluster[which(colnames(honeybee)%in% Neuron_barcode)]="Orco-Neuron";

honeybee$Annotation_subcluster <- factor(Annotation_subcluster,levels=c("Orco+Neuron","Orco-Neuron",'Epithelial cell','Glial cell',
  'Sheath cell','Obp4+support cell','Obp5+support cell',"Unannotated"))

saveRDS(honeybee,"./02_All_celltype/WNN_honeybee_integrated_all_celltype_Neuron_subcluster.rds")
