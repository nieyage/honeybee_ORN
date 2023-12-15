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
ORN <- readRDS("./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")

all_cluster<-levels(ORN$seurat_clusters)
one_classes <- c("27","28","31","32","33","35","37","38","40","41","42","43")
multiple_stop_cluster <- as.character(c(34,39))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))

# All ORN dotplot map #
barcode=colnames(ORN)
seurat_clusters=as.character(ORN$seurat_clusters)
cell_info<-data.frame(barcode,seurat_clusters)
cell_info$subcluster<-"None";
no_recluster<- c(one_classes,multiple_stop_cluster)
for (i in which(cell_info$seurat_clusters %in% no_recluster)){
    cell_info[i,]$subcluster <- cell_info[i,]$seurat_clusters
}

part1_subcluster<-readRDS("./05_ORN_cluster2/02_second_cluster/01_part1_subcluster/second_multiple_classes_part1_subcluster.rds")
cell_info_part1<-data.frame(barcode=colnames(part1_subcluster),subcluster=as.character(part1_subcluster$sub.cluster))
cell_info_part1$last_cluster<- paste("p1:",cell_info_part1$subcluster,sep="")
for (i in match(cell_info_part1$barcode,cell_info$barcode)){
  cell_info$subcluster[i] <- cell_info_part1[match(cell_info$barcode[i],cell_info_part1$barcode),]$last_cluster
}

part2_subcluster<-readRDS("./05_ORN_cluster2/02_second_cluster/02_part2_subcluster/second_multiple_classes_part2_subcluster.rds")
cell_info_part2<-data.frame(barcode=colnames(part2_subcluster),subcluster=as.character(part2_subcluster$sub.cluster))
cell_info_part2$last_cluster<- paste("p2:",cell_info_part2$subcluster,sep="")
for (i in match(cell_info_part2$barcode,cell_info$barcode)){
  cell_info$subcluster[i] <- cell_info_part2[match(cell_info$barcode[i],cell_info_part2$barcode),]$last_cluster
}

part3_subcluster<-readRDS("./05_ORN_cluster2/02_second_cluster/03_part3_subcluster/second_multiple_classes_part3_subcluster.rds")
cell_info_part3<-data.frame(barcode=colnames(part3_subcluster),subcluster=as.character(part3_subcluster$sub.cluster))
cell_info_part3$last_cluster<- paste("p3:",part3_subcluster$sub.cluster,sep="")
for (i in match(cell_info_part3$barcode,cell_info$barcode)){
  cell_info$subcluster[i] <- cell_info_part3[match(cell_info$barcode[i],cell_info_part3$barcode),]$last_cluster
}
part4_subcluster<-readRDS("./05_ORN_cluster2/02_second_cluster/04_part4_subcluster/second_multiple_classes_part4_subcluster.rds")
cell_info_part4<-data.frame(barcode=colnames(part4_subcluster),subcluster=as.character(part4_subcluster$sub.cluster))
cell_info_part4$last_cluster<- paste("p4:",part4_subcluster$sub.cluster,sep="")
for (i in match(cell_info_part4$barcode,cell_info$barcode)){
  cell_info$subcluster[i] <- cell_info_part4[match(cell_info$barcode[i],cell_info_part4$barcode),]$last_cluster
}

rownames(cell_info)<-cell_info$barcode;
cell_info<-cell_info[colnames(ORN),]
ORN$subcluster<-cell_info$subcluster;
table(ORN$subcluster)

saveRDS(ORN,"./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/Unsupervised_ORN_cluster_WNN_add_subcluster.rds")

ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/Unsupervised_ORN_cluster_WNN_add_subcluster.rds")
# remake tree for dotplot order
DefaultAssay(ORN) <- "integratedRNA_onecluster"
Idents(ORN)<-ORN$subcluster
object<-ORN;
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
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/cluster-ORN-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/cluster-ORN-tree-cosine-nocircular.pdf",width=8,height=12)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()

#cluster order by tree
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

# plot dotplot by SCT value 
Idents(ORN)<-factor(ORN$subcluster,levels=cluster_order)
DefaultAssay(ORN)<-"raw_RNA"
ORN<- SCTransform(ORN,return.only.var.genes = FALSE,variable.features.n=11794,verbose = FALSE)
DefaultAssay(ORN) <- "SCT"

# Plot receptor 
receptor.dot <- DotPlot(ORN, features = all_receptor_gene) + #scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
#p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-receptor.dot$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled >= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))

DefaultAssay(ORN)<-"SCT"
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/UnSupervised_WNN_dotplot-allfeature_orderbytree.pdf",width=28, height=20)
p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
p2<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') +  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) +
scale_color_gradientn(colours = c('#008080', '#FF00FF'),  name = 'Average\nexpression', oob = scales::squish)
p2
dev.off()

pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/UnSupervised_WNN_dotplot-signif-feature_orderbytree.pdf",width=25, height=20)
p<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
p2<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') +  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) +
scale_color_gradientn(colours = c('#008080', '#FF00FF'),  name = 'Average\nexpression', oob = scales::squish)
p2
dev.off()

# merge some subcluster (same OR feature and can not distinguish by trans dist )

# show a typical combination;
# select a beautiful track to show :
library(ggplot2)
library(gggenes)
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
library(RColorBrewer)
log2FCdata<-data.frame();
DefaultAssay(ORN)<- "integratedRNA_onecluster"

 cluster_info<- c("p1:2","p1:18")
  obj<-subset(ORN,idents=cluster_info);

# p1 cell cosine simility heatmap 
    embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
    trans_dist <- 1-cosine(t(embeddings))
    barcode_label_pheatmap<-data.frame(label=obj$subcluster)
    rownames(barcode_label_pheatmap)<- colnames(obj)
    col<-brewer.pal(12,"Set3")[1:length(unique(barcode_label_pheatmap$label))]
    names(col)<-unique(barcode_label_pheatmap$label)
    ann_colors= list(label = col)
 

pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/combination_data_multi_OR_with_powerful.pdf",width=6,height=6)
 pheatmap(trans_dist,
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=F,
             show_colnames=F)
dev.off()
# 
ORN$subcluster[which(ORN$subcluster=="p4:14")]="p3:22"
ORN$subcluster[which(ORN$subcluster=="p4:4_1")]="p1:12"
ORN$subcluster[which(ORN$subcluster=="p2:9")]="p2:8"
ORN$subcluster[which(ORN$subcluster=="p1:18")]="p1:2"

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )
ORN$subcluster<-factor(ORN$subcluster,levels=cluster_order)
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/ORN_cluster_WNN_all_add_subcluster_info.pdf",width=12,height=6)
###cluster
DimPlot(ORN, cols=c(myUmapcolors,myUmapcolors,myUmapcolors), reduction = "tsne.rna",pt.size=0.02,  label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
dev.off();
saveRDS(ORN,"./05_ORN_cluster2/02_second_cluster/Unsupervised_ORN_cluster_WNN_add_subcluster_merge_same.rds")

ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/Unsupervised_ORN_cluster_WNN_add_subcluster_merge_same.rds")
DefaultAssay(ORN)<- "SCT"
receptor.dot <- DotPlot(ORN, features = all_receptor_gene) + #scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
#p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-receptor.dot$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>40&&dotplot_data[i,]$avg.exp.scaled >= 2.5){dotplot_data[i,]$state="Yes"};
}

dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[which(dotplot_data$avg.exp>1),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))

# remove the cluster have no signif dotplot feature and cell number is lower than 
cluster_number<- as.data.frame(table(ORN$subcluster))
cluster_number<- cluster_number[order(cluster_number$Freq,decreasing=T),]
cluster_number<- cluster_number[cluster_number$Var1%in%dotplot_data$id,]
cluster_number<- cluster_number[cluster_number$Freq>30,]

keep_cluster<- as.character(cluster_number$Var1)
## cluster trans_dist tree color by useful 
# remake tree for dotplot order

#cluster tree 
Idents(ORN)<- ORN$subcluster
ORN_withpower <-subset(ORN,idents=keep_cluster)
# merge same OR cluster
ORN_withpower$subcluster[which(ORN_withpower$subcluster=="p4:0_2")]="p4:0_1"

DefaultAssay(ORN_withpower) <- "integratedRNA_onecluster"

object<- ORN_withpower
DefaultAssay(object)<-"RNA"
obj<-FindVariableFeatures(object, selection.method = "vst")
top <- head(VariableFeatures(obj),500)
obj<-ScaleData(obj,rownames(obj))
ORN_withpower <- RunPCA(obj,features=c(all_receptor_gene,top),reduction.name="obj_features_pca") 
object<- ORN_withpower
embeddings <- Embeddings(object = object, reduction = "obj_features_pca")[,1:50]
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

#cluster order by tree
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)
cluster_order
ORN_withpower$subcluster<-factor(ORN_withpower$subcluster,levels=cluster_order)
Idents(ORN_withpower)<-factor(ORN_withpower$subcluster,levels=cluster_order)

ORN_withpower$cell_group <- as.numeric(factor(ORN_withpower$subcluster))
Idents(ORN_withpower)<- factor(ORN_withpower$cell_group,levels=as.character(1:60))


# plot the dotplot 
DefaultAssay(ORN_withpower)<-"SCT"
p<-DotPlot(ORN_withpower,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>40&&dotplot_data[i,]$avg.exp.scaled>= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
#dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% dotplot_data$features.plot])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]

dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))
write.csv(dotplot_data,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
DefaultAssay(ORN_withpower)<-"SCT"
pdf("./00_Figure/Fig2/Fig2D-b-remove_nopower-dotplot-orderbytree_lastest.pdf",width=25, height=14)
p<-DotPlot(ORN_withpower,features = dotplot_feature,cols=c("lightgrey","#0000CC")) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
DefaultAssay(ORN_withpower)<-"raw_RNA"

saveRDS(ORN_withpower,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_latest.rds")


# part VF
# For part1_subcluster
DefaultAssay(part1_subcluster) <- "RNA"
part1_subcluster <- FindVariableFeatures(part1_subcluster, selection.method = "vst",features=100)
part1_subcluster_top100 <- head(VariableFeatures(part1_subcluster),100)

# For part2_subcluster
DefaultAssay(part2_subcluster) <- "RNA"
part2_subcluster <- FindVariableFeatures(part2_subcluster, selection.method = "vst",features=100)
part2_subcluster_top100 <- head(VariableFeatures(part2_subcluster),100)

# For part3_subcluster
DefaultAssay(part3_subcluster) <- "RNA"
part3_subcluster <- FindVariableFeatures(part3_subcluster, selection.method = "vst",features=200)
part3_subcluster_top100 <- head(VariableFeatures(part3_subcluster),100)
# 4 parts to subcluster 
# For part4_subcluster
# vst
DefaultAssay(part4_subcluster) <- "RNA"
part4_subcluster <- FindVariableFeatures(part4_subcluster, selection.method = "vst")
part4_subcluster_top100 <- head(VariableFeatures(part4_subcluster),100)

part_VF<- c(part1_subcluster_top100,part2_subcluster_top100,part3_subcluster_top100,part4_subcluster_top100)
#RNA analysis
DefaultAssay(ORN_withpower) <- "integratedRNA_onecluster"
ORN_withpower <- RunPCA(ORN_withpower,features= part_VF )
pdf("./05_ORN_cluster2/02_second_cluster/ORN_withpower_ElbowPlot_top100.pdf")
ElbowPlot(ORN_withpower,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
ORN_withpower <- RunTSNE(
  object = ORN_withpower,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:30
)
ORN_withpower <- FindNeighbors(object = ORN_withpower, reduction = 'pca', dims = 1:30)
ORN_withpower <- FindClusters( object = ORN_withpower, verbose = FALSE, resolution =2,algorithm = 3)
table(ORN_withpower$seurat_clusters,ORN_withpower$subcluster)
pdf("./05_ORN_cluster2/02_second_cluster/last_ORN_withpower_subcluster.pdf",width=9,height=5)
#DimPlot(ORN_withpower, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(ORN_withpower, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "subcluster")
dev.off()

Idents(ORN_withpower)<- ORN_withpower$subcluster
#cluster tree 
DefaultAssay(ORN_withpower) <- "integratedRNA_onecluster"
object<- ORN_withpower
embeddings <- Embeddings(object = object, reduction = "obj_features_pca")[,1:50]
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
#cluster order by tree
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

ORN_withpower$subcluster<-factor(ORN_withpower$subcluster,levels=cluster_order)
Idents(ORN_withpower)<-factor(ORN_withpower$subcluster,levels=cluster_order)

ORN_withpower$cell_group <- as.numeric(factor(ORN_withpower$subcluster))
Idents(ORN_withpower)<- factor(ORN_withpower$cell_group,levels=as.character(1:60))

DefaultAssay(ORN_withpower)<-"raw_RNA"
ORN_withpower$cell_group[which(ORN_withpower$cell_group=="46")]="44"
ORN_withpower$cell_group <- as.numeric(factor(ORN_withpower$cell_group))
ORN_withpower$cell_group<- factor(ORN_withpower$cell_group,level=as.character(1:59))

Idents(ORN_withpower)<- ORN_withpower$cell_group
# plot the dotplot 
DefaultAssay(ORN_withpower)<-"SCT"
p<-DotPlot(ORN_withpower,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>40&&dotplot_data[i,]$avg.exp.scaled>= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
#dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% dotplot_data$features.plot])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]

dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))
write.csv(dotplot_data,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
DefaultAssay(ORN_withpower)<-"SCT"
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

pdf("./00_Figure/Fig2/Fig2D-b-remove_nopower-dotplot-orderbytree_lastest.pdf",width=25, height=14)
p<-DotPlot(ORN_withpower,features = dotplot_feature,cols=c(solarExtra[4],solarExtra[8])) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
saveRDS(ORN_withpower,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_latest.rds")




