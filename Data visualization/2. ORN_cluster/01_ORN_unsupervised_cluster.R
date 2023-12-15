#Step1: ORN cluster by OR gene 
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype_Neuron_subcluster.rds")
DefaultAssay(honeybee) <- "raw_RNA"
#Step1: extract Orco highly exp cluster#####
Idents(honeybee)<- honeybee$Annotation_subcluster
ORN <- subset(honeybee,idents=as.character("Orco+Neuron"))
ORN.list <- SplitObject(ORN, split.by = "orig.ident");
for (i in 1:length(ORN.list)) {
  DefaultAssay(ORN.list[[i]]) <- "raw_RNA"
  ORN.list[[i]] <- SCTransform(ORN.list[[i]],return.only.var.genes = FALSE,variable.features.n=10000,verbose = FALSE)
}
ORN.features <- SelectIntegrationFeatures(object.list = ORN.list, nfeatures = 4000);
ORN.list <- PrepSCTIntegration(object.list = ORN.list, anchor.features = ORN.features)
#integrate RNA using rpca
Neuron_list <- lapply(
  X = ORN.list,
  FUN = ScaleData,
  features = ORN.features,
  verbose = FALSE
)
Neuron_list <- lapply(
  X = Neuron_list,
  FUN = RunPCA,
  features = ORN.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list = Neuron_list,
  normalization.method = "SCT",
  anchor.features = ORN.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)
all_features <- lapply(Neuron_list, row.names) %>% Reduce(intersect,.) 

ORN <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA_ORN",
  features.to.integrate = all_features,
  dims = 1:30
)

#run LSI on new seurat object with integrated RNA assay
DefaultAssay(ORN) <- "ATAC"
ORN <- RunTFIDF(ORN)
ORN <- FindTopFeatures(ORN, min.cutoff = "q25")
ORN <- RunSVD(ORN,n=25)
#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
ORN_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI_ORN",
  reductions = ORN@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
ORN@reductions$integratedLSI_ORN <- ORN_atac@reductions$integratedLSI_ORN
#####done integrate ATAC and RNA ################
# RNA analysis
library(dplyr)
#first cluster for ORN
# vst
DefaultAssay(ORN) <- "RNA"
ORN <- FindVariableFeatures(ORN, selection.method = "vst")
top200 <- head(VariableFeatures(ORN),200)
hvf.info <- HVFInfo(object = ORN,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster2/01_first_cluster/Find_var_RNA.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
DefaultAssay(ORN) <- "integratedRNA_ORN"
ORN <- ScaleData(ORN,features=rownames(ORN))
ORN <- RunPCA(ORN)
ORN <- RunUMAP(ORN,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
ORN <- RunTSNE(ORN,dims = 1:50, min.dist = 0.001,check_duplicates = FALSE, reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')
ORN <- FindNeighbors(object = ORN, reduction = 'pca', dims = 1:50)
ORN <- FindClusters(object = ORN, verbose = FALSE,  resolution =6, algorithm = 3)
table(ORN$seurat_clusters)

Idents(ORN)<-ORN$orig.ident
ORN$orig.ident<-factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )

pdf("./05_ORN_cluster2/01_first_cluster/raw_ORN_cluster.pdf",width=6,height=5)
DimPlot(ORN, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
DimPlot(ORN, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "umap.rna",group.by = "seurat_clusters")
#DimPlot(ORN, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "wnn.umap",group.by = "seurat_clusters")
dev.off()

# plot the dotplot by all receptor gene 
DefaultAssay(ORN) <- "SCT"
Idents(ORN)<-ORN$seurat_clusters
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="IR",]$gene_name)
#IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
pdf("./05_ORN_cluster2/01_first_cluster/first_raw_OR_dotplot_rawRNA-New-filtration-standard.pdf",width=30, height=12)
p1<-DotPlot(ORN, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p1;
dev.off()

DefaultAssay(ORN) <- "SCT"
Idents(ORN)<-ORN$seurat_clusters
pdf('./05_ORN_cluster2/01_first_cluster/first_coreceptor_VlnPlot-New-filtration-standard.pdf',width=15, height=10)
print( VlnPlot(ORN, features = Orco, ncol = 1, pt.size = 0) )
dev.off()

saveRDS(ORN,"./05_ORN_cluster2/01_first_cluster/ORN_integrated_antenna_first_cluster.rds")

#Step2: remove cluster without Orco 
ORN <- readRDS("./05_ORN_cluster2/01_first_cluster/ORN_integrated_antenna_first_cluster.rds")
all_cluster<-levels(ORN$seurat_clusters)
rm_cluster<-c("36")
onecluster <- subset(ORN,idents=setdiff(all_cluster,rm_cluster))
onecluster.list <- SplitObject(onecluster, split.by = "orig.ident")
for (i in 1:length(onecluster.list)) {
  DefaultAssay(onecluster.list[[i]]) <- "raw_RNA"
  onecluster.list[[i]] <- SCTransform(onecluster.list[[i]],return.only.var.genes = FALSE,variable.features.n=3000,verbose = FALSE)
}
onecluster.features <- SelectIntegrationFeatures(object.list = onecluster.list, nfeatures = 3000)
onecluster.list <- PrepSCTIntegration(object.list = onecluster.list, anchor.features = onecluster.features)

#integrate RNA using rpca
onecluster.list <- lapply(
  X = onecluster.list,
  FUN = ScaleData,
  features = onecluster.features,
  verbose = FALSE
)
Neuron_list <- lapply(
  X = onecluster.list,
  FUN = RunPCA,
  features = onecluster.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list = Neuron_list,
  normalization.method = "SCT",
  anchor.features = onecluster.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)

all_features <- lapply(Neuron_list, row.names) %>% Reduce(intersect,.) 
onecluster <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA_onecluster",
  features.to.integrate = all_features,
  dims = 1:30
)

#run LSI on new seurat object with integrated RNA assay
DefaultAssay(onecluster) <- "ATAC"
onecluster <- RunTFIDF(onecluster)
onecluster <- FindTopFeatures(onecluster, min.cutoff = "q15")
onecluster <- RunSVD(onecluster)
#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
onecluster_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI_onecluster",
  reductions = onecluster@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
onecluster@reductions$integratedLSI_onecluster <- onecluster_atac@reductions$integratedLSI_onecluster

DefaultAssay(onecluster) <- "integratedRNA_onecluster"
#onecluster <- FindVariableFeatures(onecluster, selection.method = "vst")
top500 <- head(VariableFeatures(onecluster),500)
#hvf.info <- HVFInfo(object = onecluster,status = TRUE)
#hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
#hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
#hvf.info$variable$vst.variable[501:length(hvf.info$variable$vst.variable)]=0
#hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
#var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
#hvf.info$var.status <- var.status
#pdf("./05_ORN_cluster2/02_second_cluster/Find_var_RNA_top500.pdf",width=20,height=6)
#ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
#geom_point(aes(color=var.status))+xlim(0,10)
#dev.off()

# RNA analysis
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
onecluster <- ScaleData(onecluster, features =rownames(onecluster),verbose = FALSE)
onecluster <- RunPCA(onecluster,features=top500) 
pdf("./05_ORN_cluster2/02_second_cluster/ElbowPlot.pdf")
ElbowPlot(onecluster,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
onecluster <- FindNeighbors(object = onecluster, reduction = 'pca', dims = 1:50)
onecluster <- FindClusters( object = onecluster, verbose = FALSE, resolution =6,algorithm = 3)
table(onecluster$seurat_clusters)
onecluster <- RunTSNE(
  object = onecluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)

pdf("./05_ORN_cluster2/02_second_cluster/second_ORN_cluster_WNN_top500.pdf",width=6,height=5)
DimPlot(onecluster, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(onecluster, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
# plot the dotplot by all receptor gene 
DefaultAssay(onecluster) <- "SCT"
Idents(onecluster)<-onecluster$seurat_clusters
pdf("./05_ORN_cluster2/02_second_cluster/second_raw_OR_dotplot_rawRNA_top500.pdf",width=25, height=12)
p<-DotPlot(onecluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
pdf('./05_ORN_cluster2/02_second_cluster/second_coreceptor_VlnPlot_WNN_top500.pdf',width=15, height=10)
print( VlnPlot(onecluster, features =Orco, ncol = 1, pt.size = 0) )
dev.off()
saveRDS(onecluster,"./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")

onecluster <- readRDS("./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")

# Step3: distinguish the OR cluster 
dotplot_data<- p$data
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled > 1){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- c("27","28","31","32","33","35","37","38","40","41","42","43")

# each multi OR cluster 
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster2/02_second_cluster/distinguish_multi_OR_for_secondcluster.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(onecluster,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
# add max_exp OR label for each cell
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
  for (i in 1:length(colnames(ORN_matrix))){
    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
}
  }
barcode_label<-barcode_label[barcode_label$label!="NA",]
barcode_label<-barcode_label[order(barcode_label$label),]
##cell cosine simility heatmap 
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col <- myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
names(col)<-unique(barcode_label_pheatmap$OR)
ann_colors= list(OR = col)
p1<-pheatmap(trans_dist,
       cluster_cols = TRUE,
       cluster_rows = TRUE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
#calculate the cosine simility within group and between groups;
rownames(barcode_label)<-barcode_label$barcode
within_group<-c()
between_group<-c()
for (i in 1:nrow(trans_dist)){
  for (j in 1:ncol(trans_dist)){
    if(i!=j){
      if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
        within_group<-c(within_group,trans_dist[i,j]);
      }
      else{between_group<-c(between_group,trans_dist[i,j])}
    }
  }
}
# calculate the FC 
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame(cluster,log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)

# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
              geom_density(alpha=.25) + theme_classic() 
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#method = "t.test"
# plot OR heatmap 
DefaultAssay(obj)<-"raw_RNA";
obj<-ScaleData(obj,features=all_receptor_gene);
#p4<-dittoHeatmap(obj,obj_features,slot ="scale.data",cluster_cols=T,scaled.to.max = TRUE)

top_right<-plot_grid(p2,p3,labels = c("B","C"),ncol = 1)
#right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
last<-plot_grid(p1$gtable, top_right, labels = c('A', ''),rel_widths = c(1, 1),label_size = 12, ncol = 2)
title <- ggdraw() + 
  draw_label(
    paste("Cluster",cluster,":","log2FC=",log2FC),
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
print(add_title)
}
dev.off()

# Step4: select the cluster to subcluster 
# first subcluster 
# multiple_stop_cluster #
> log2FCdata
   cluster     log2FC        pvalue
1        1 0.04520271  9.085006e-14
2        2 0.06147687  2.189160e-12
3        5 0.30049926  0.000000e+00
4        7 0.14424370  1.216990e-50
5        8 0.16565959 2.026697e-107
6       10 0.11568475  4.797983e-17
7       11 0.05662896  3.874407e-14
8       12 0.11587880  4.172628e-22
9       15 0.19443849 2.139598e-183
10      18 0.28002031 6.026136e-134
11      19 0.25138096 2.171005e-228
12      20 0.19689640  2.869033e-98
13      21 0.01873434  3.036023e-03
14      22 0.01626887  2.258964e-01
15      23 0.18505486  3.131805e-17
16      24 0.52144780 1.782721e-209
17      25 0.07560773  6.449481e-06
18      30 0.35909064 8.850859e-150
19      34 0.09992559  1.471053e-12
20      36 0.09661868  6.915426e-10
21      39 0.01798649  2.439830e-02
22      44 0.14497994  1.785379e-04
23      45 0.23904754  4.790468e-01


one_classes <- c("27","28","31","32","33","35","37","38","40","41","42","43")
multiple_stop_cluster <- as.character(c(34,39))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))
