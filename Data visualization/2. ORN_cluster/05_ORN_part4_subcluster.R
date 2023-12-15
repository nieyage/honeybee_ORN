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
# make the trans dist tree 
object <- onecluster
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

pdf("./05_ORN_cluster2/02_second_cluster/second_cluster-ORN-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()
pdf("./05_ORN_cluster2/02_second_cluster/second_cluster-ORN-tree-cosine-nocircular.pdf",width=8,height=12)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()
m<- ggtree(data.tree) + geom_tiplab()+ geom_treescale()

cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

 [1] "11" "24" "23" "26" "37" "1"  "4"  "38" "39" "20" "7"  "30" "12" "34" "6" 
[16] "10" "28" "29" "21" "36" "31" "25" "2"  "22" "16" "18" "0"  "3"  "8"  "42"
[31] "13" "19" "27" "41" "15" "14" "5"  "33" "32" "35" "43" "40" "9"  "17"

one_classes <- c("29 ","28","31","33","36","39","41")
multiple_stop_cluster <- as.character(c(25,34,38))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))

# 4 parts to subcluster 
# For part4_subcluster
# vst
DefaultAssay(part4_subcluster) <- "RNA"
part4_subcluster <- FindVariableFeatures(part4_subcluster, selection.method = "vst")
top100 <- head(VariableFeatures(part4_subcluster),100)
all_receptor_gene[which(all_receptor_gene%in% top100)]
hvf.info <- HVFInfo(object = part4_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/Find_var_RNA_top100.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
#RNA analysis
DefaultAssay(part4_subcluster) <- "integratedRNA_onecluster"
part4_subcluster <- RunPCA(part4_subcluster,features= top100 )
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/ElbowPlot_top100.pdf")
ElbowPlot(part4_subcluster,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
part4_subcluster <- RunTSNE(
  object = part4_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:40
)
part4_subcluster <- FindNeighbors(object = part4_subcluster, reduction = 'pca', dims = 1:40)
part4_subcluster <- FindClusters( object = part4_subcluster, verbose = FALSE, resolution =2,algorithm = 3)
table(part4_subcluster$seurat_clusters)
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/part4_subcluster.pdf",width=6,height=5)
DimPlot(part4_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part4_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
DefaultAssay(part4_subcluster) <- "SCT"
Idents(part4_subcluster)<-part4_subcluster$seurat_clusters
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/part4_subcluster_OR_dotplot_rawRNA.pdf",width=30, height=8)
p<-DotPlot(part4_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
saveRDS(part4_subcluster,"./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/second_multiple_classes_part4_subcluster.rds");

# Step4: select the cluster to subcluster 
# part4_subcluster distinguish OR pipeline 
p<-DotPlot(part4_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp >= 1){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- as.character(cluster_info[cluster_info$Freq==1,1])

# distinguish_multi_OR
DefaultAssay(part4_subcluster)<- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/part4_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(part4_subcluster,idents=cluster);
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
  cluster     log2FC        pvalue
1       0  0.2461361 1.711015e-228
2       3  0.5134898  0.000000e+00
3       5  0.1562849  4.898015e-51
4       6  0.3487910  3.753938e-60
5      11  0.3476158  3.715402e-20
6      12 -0.3740238  1.652260e-05


## perform sub-clustering on cluster to find additional structure
Idents(part4_subcluster)<- part4_subcluster$seurat_clusters
part4_subcluster<-FindSubCluster(part4_subcluster,0,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.6,algorithm = 1)
Idents(part4_subcluster)<-part4_subcluster$sub.cluster
part4_subcluster<-FindSubCluster(part4_subcluster,3,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
Idents(part4_subcluster)<-part4_subcluster$sub.cluster
part4_subcluster<-FindSubCluster(part4_subcluster,4,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
table(part4_subcluster$sub.cluster)
DefaultAssay(part4_subcluster) <- "SCT"
Idents(part4_subcluster)<-part4_subcluster$sub.cluster
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/OR_dotplot_rawRNA_part4_subcluster_recluster.pdf",width=30, height=8)
p<-DotPlot(part4_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(part4_subcluster,"./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/second_multiple_classes_part4_subcluster.rds");

# FigS2 B 
Idents(part4_subcluster)<-part4_subcluster$sub.cluster
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/part4_subcluster_last.pdf",width=6,height=5)
DimPlot(part4_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part4_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna")
dev.off()

cluster_cellnumber<-as.data.frame(table(Idents(part4_subcluster)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(part4_subcluster))
cluster_cellnumber$color<-myUmapcolors[1:length(levels(part4_subcluster))]
cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/part4_subcluster_cellnumber.pdf",width=4,height=8)
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

DefaultAssay(part4_subcluster) <- "SCT"
p<-DotPlot(part4_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp >= 1){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];


DefaultAssay(part4_subcluster) <- "raw_RNA"
pdf("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/dotplot_part4_subcluster_last.pdf",width=10, height=8)
p<-DotPlot(part4_subcluster, features = unique(c(Orco,rev(as.character(dotplot_data$features.plot))))) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

