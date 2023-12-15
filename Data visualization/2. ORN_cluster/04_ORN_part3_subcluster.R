#Step1: ORN cluster by OR gene 
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))

onecluster <- readRDS("./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
Idents(onecluster) <- "seurat_clusters";
all_cluster<-levels(onecluster$seurat_clusters)
one_classes <- c("27","28","31","32","33","35","37","38","40","41","42","43")
multiple_stop_cluster <- as.character(c(34,39))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))


# make the trans dist tree 
object <- subset(onecluster,idents=need2subcluster)

embeddings <- Embeddings(object = object, reduction = "pca")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("./05_ORN_cluster2/02_second_cluster/second_cluster-ORN-tree-cosine_need2subcluster.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()

m<- ggtree(data.tree) + geom_tiplab()+ geom_treescale()

cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

> cluster_order
 [1] "1"  "15" "19" "20" "3"  "4"  "25" "30" "36" "24" "29" "23" "12" "17" "0" 
[16] "10" "5"  "2"  "14" "26" "13" "21" "11" "18" "16" "22" "7"  "8"  "6"  "9" 


# 4 parts to subcluster 
part1 <- cluster_order[1:8]
part2 <- cluster_order[9:16]
part3 <- cluster_order[17:23]
part4 <- cluster_order[24:32]
need2subcluster 
part1 <- part1[which(part1%in%need2subcluster)]
part2 <- part2[which(part2%in%need2subcluster)]
part3 <- part3[which(part3%in%need2subcluster)]
part4 <- part4[which(part4%in%need2subcluster)]

part1_subcluster <- subset(onecluster,idents=part1)
part2_subcluster <- subset(onecluster,idents=part2)
part3_subcluster <- subset(onecluster,idents=part3)
part4_subcluster <- subset(onecluster,idents=part4)
part1_subcluster
part2_subcluster
part3_subcluster
part4_subcluster

# For part3_subcluster
# vst
DefaultAssay(part3_subcluster) <- "RNA"
part3_subcluster <- FindVariableFeatures(part3_subcluster, selection.method = "vst",features=200)
top100 <- head(VariableFeatures(part3_subcluster),100)
all_receptor_gene[which(all_receptor_gene%in% top200)]
hvf.info <- HVFInfo(object = part3_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/Find_var_RNA_top100.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
#RNA analysis

DefaultAssay(part3_subcluster) <- "integratedRNA_onecluster"
part3_subcluster <- RunPCA(part3_subcluster,features= top100 )
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/ElbowPlot_top100.pdf")
ElbowPlot(part3_subcluster,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
part3_subcluster <- RunTSNE(
  object = part3_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)
part3_subcluster <- FindNeighbors(object = part3_subcluster, reduction = 'pca', dims = 1:50)
part3_subcluster <- FindClusters( object = part3_subcluster, verbose = FALSE, resolution =2.5,algorithm = 3)
table(part3_subcluster$seurat_clusters)

pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/part3_subcluster.pdf",width=6,height=5)
DimPlot(part3_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part3_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
DefaultAssay(part3_subcluster) <- "SCT"
Idents(part3_subcluster)<-part3_subcluster$seurat_clusters
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/part3_subcluster_OR_dotplot_rawRNA.pdf",width=30, height=8)
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()


saveRDS(part3_subcluster,"./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/second_multiple_classes_part3_subcluster.rds");

# Step4: select the cluster to subcluster 
# part3_subcluster distinguish OR pipeline 
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp > 1){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- as.character(cluster_info[cluster_info$Freq==1,1])

# distinguish_multi_OR
DefaultAssay(part3_subcluster)<- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/part3_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(part3_subcluster,idents=cluster);
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
col<-myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
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

> log2FCdata


## perform sub-clustering on cluster to find additional structure
Idents(part3_subcluster)<- part3_subcluster$seurat_clusters
part3_subcluster<-FindSubCluster(part3_subcluster,4,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.3,algorithm = 1)
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
part3_subcluster<-FindSubCluster(part3_subcluster,5,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
table(part3_subcluster$sub.cluster)
DefaultAssay(part3_subcluster) <- "SCT"
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/OR_dotplot_rawRNA_part3_subcluster_recluster.pdf",width=30, height=8)
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(part3_subcluster,"./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/second_multiple_classes_part3_subcluster.rds");

# FigS2 B 
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/part3_subcluster_last.pdf",width=6,height=5)
DimPlot(part3_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part3_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna")
dev.off()

cluster_cellnumber<-as.data.frame(table(Idents(part3_subcluster)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(part3_subcluster))
cluster_cellnumber$color<-myUmapcolors[1:length(levels(part3_subcluster))]
cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/part3_subcluster_cellnumber.pdf",width=4,height=8)
p<-ggplot(data = cluster_cellnumber, aes_string(x = "cluster", y = "number", 
        fill = "cluster")) +  xlab(" ") + ylab("# of cells") + 
        scale_fill_manual(values = cluster_cellnumber$color) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5))+coord_flip();
p
#add gene number in plot 
p+geom_text(aes(label = number), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

DefaultAssay(part3_subcluster) <- "raw_RNA"
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp >= 1){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]

DefaultAssay(part3_subcluster) <- "raw_RNA"
pdf("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/dotplot_part3_subcluster_last.pdf",width=10, height=8)
p<-DotPlot(part3_subcluster, features = unique(c(Orco,rev(as.character(dotplot_data$features.plot))))) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

