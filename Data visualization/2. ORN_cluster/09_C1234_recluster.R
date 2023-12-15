# Fig3: Cis-Elements Orchestrating Olfactory Receptor Co-Expression 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)

set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C"  , "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
# get the OR gene and the promoter reads 
obj<- readRDS("./05_ORN_cluster2/05_combination_group_recluster/obj_recluster.rds")
obj$group<- Idents(obj)

# recall peak 
DefaultAssay(obj)<-"ATAC"
peak<-CallPeaks(
       obj,
       group.by = "group",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="./00_Figure/Fig4/",
       combine.peaks=TRUE
)

macs2_counts <- FeatureMatrix(
     fragments = Fragments(obj),
     features = peak,
     cells = colnames(obj)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# create a new assay using the MACS2 peak set and add it to the Seurat object
obj[["peaks_obj"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(obj),
  annotation = Annotation(obj)
)

ORN_count<-obj@assays$RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)

DefaultAssay(obj)<-"peaks_ORN_subcluster"
obj_peaks<- c( rownames(obj)[grep("Group15-698....-.......",rownames(obj))],rownames(obj)[grep("Group15-697....-.......",rownames(obj))],rownames(obj)[grep("Group15-695....-.......",rownames(obj))],
  rownames(obj)[grep("Group15-696....-.......",rownames(obj))])

ORN_count<-obj@assays$peaks_ORN_subcluster
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_peaks),]
ORN_matrix2<-as.matrix(ORN_count)

library(GenomicRanges)

# 创建一个包含你提供的字符串的向量
strings <- c("Group15-6982829-6983491", "Group15-6970090-6970517",
             "Group15-6971874-6972439", "Group15-6956351-6957167",
             "Group15-6962048-6962247", "Group15-6964370-6964569","Group15-6966000-6967000")

# 将字符串转换成 GRanges 对象
granges_obj <- GRanges(
  seqnames = rep("Group15", length(strings)),  # 假设所有序列都来自 "Group15"
  ranges = IRanges(
    start = c(6982829,6970090,6971874,6956351,6962048,6964370,6966000),  # 提取起始位置
    end = c(6983491,6970517,6972439,6957167,6962247,6964569,6967000)     # 提取结束位置
  )
)


macs2_counts <- FeatureMatrix(
     fragments = Fragments(obj),
     features = granges_obj,
     cells = colnames(obj)
     )   
last_data<- as.data.frame(t(rbind(ORN_matrix,ORN_matrix2,macs2_counts)))

last_data$cluster<- obj$group

write.csv(last_data,"./00_Figure/Fig4/Fig4-last-data.csv")

last_data<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig4/Fig4-last-data-V2.csv")
rownames(last_data)<- last_data$X
last_data<- last_data[colnames(obj),]
obj$group_manully<- last_data$cluster

Idents(obj)<- obj$group_manully
# show a typical combination;
# select a beautiful track to show :
log2FCdata<-data.frame();
DefaultAssay(obj) <- "peaks_ORN_subcluster"
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./00_Figure/Fig4/combination_group2_cluster_WNN.pdf",width=15,height=5)
###cluster
DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "tsne.rna", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
dev.off()

# color set 
color_for_cluster<- c("#4BA9D1",my47colors[6:7])
color_for_group <- c("#476D87","#E95C59")

ORN_count<-obj@assays$SCT
barcode_label<-data.frame(barcode=colnames(obj),label=obj$group_manully)
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-barcode_label[colnames(ORN_matrix),]
barcode_label<-barcode_label[order(barcode_label$label),]
# p1 cell cosine simility heatmap 
#DefaultAssay(obj)<-"integratedRNA_onecluster"
#obj<- ScaleData(obj)
#obj <- RunPCA(obj,reduction.name="pca") 
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:20]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(label=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode

#trans_dist<-trans_dist[rownames(barcode_label_pheatmap),rownames(barcode_label_pheatmap)]
p1<-pheatmap(trans_dist,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         #annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F)
# calculate the cosine simility within group and between groups;
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
if(!is.null(median(between_group))){
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame("1",log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)
}



# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+guides(color='none')+
              geom_density(alpha=.25) + theme_classic()+ scale_fill_manual(values=color_for_group)+ scale_color_manual(values=color_for_group)
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type",width=0.6,) +stat_compare_means()+guides(color = "none")+ scale_color_manual(values=color_for_group)
#raw counts heatmap 
# Heat map of expression  value 
    DefaultAssay(obj)<-"raw_RNA"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[c( "Or63-b","LOC410603" ,   "LOC107963999"   ,    "LOC100578045"),])))
    obj_data<-obj_data[rownames(barcode_label),]
    p4<-pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "#CC0000"))(100),
             #annotation_col = barcode_label,
             #annotation_colors = ann_colors,
             #annotation_row = barcode_label,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
    top_right<-plot_grid(p2,p3,labels = c(" "," "),rel_widths = c(2, 1))
    right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" "," "))
    last<-plot_grid(p1$gtable, right, labels = c(' ', ''), label_size = 12, ncol = 2)
    title <- ggdraw() + 
      draw_label(
        paste("combination","log2FC=",log2FC),
        fontface = 'bold',
        x = 0,
        hjust = 0
      )
    add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )

pdf("./00_Figure/Fig4/combination_group2_trans_exp.pdf",width=12,height=5)
print(add_title)
dev.off()

# Fig4H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(obj)<-"peaks_ORN_subcluster"
obj$group_manully<- factor(obj$group_manully,levels=c("C1","C2","C3","C4"))
Idents(obj)<- obj$group_manully
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp######
  idents.plot <- Idents(obj)
  # plot region 
  start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
  ranges.show <- paste(seq,start,end,sep="-")
  col<- color_for_cluster
  p1<-CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 50,
    links=F
  )

pdf("./00_Figure/Fig4/combination_group2_trackplot.pdf",width=10,height=10)
#p1<-p1& scale_fill_manual(values=col)
CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 50,
    links=F
  )
CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 10,
    links=F
  )
CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 200,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 50,
    links=F
  )

dev.off()
saveRDS(obj,"./00_Figure/Fig4/Fig4-last-data-obj.rds")

# show the UMAP and the transcript dist heatmap 
# Fig4B 
color_for_cluster<- c("#4BA9D1",my47colors[7:9])


pdf("./00_Figure/Fig4/Fig4B-UMAP-cluster.pdf", width=12, height=4)
p1 <- DimPlot(obj, cols=color_for_cluster,pt.size = 0.8, reduction = "wnn.umap", label = TRUE, label.size = 3.5, repel = TRUE)
p2 <- DimPlot(obj, cols=color_for_cluster,pt.size = 0.8, reduction = "umap.atac", label = TRUE, label.size = 3.5, repel = TRUE)
p3 <- DimPlot(obj, cols=color_for_cluster,pt.size = 0.8, reduction = "umap.rna", label = TRUE, label.size = 3.5, repel = TRUE)
p3 + p2 + p1 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# show a typical combination;
# select a beautiful track to show :
library(ggplot2)
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
library(RColorBrewer)
log2FCdata<-data.frame();

obj$Annotation<- obj$group_manully
ORN_count<-obj@assays$SCT
barcode_label<-data.frame(barcode=colnames(obj),label=obj$group_manully)
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-barcode_label[colnames(ORN_matrix),]
barcode_label<-barcode_label[order(barcode_label$label),]
# p1 cell cosine simility heatmap 
#DefaultAssay(obj)<-"integratedRNA_onecluster"
#obj<- ScaleData(obj)
#obj <- RunPCA(obj,reduction.name="pca") 
embeddings <- Embeddings(object = obj, reduction = "obj_features_pca")
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(label=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)
#trans_dist<-trans_dist[rownames(barcode_label_pheatmap),rownames(barcode_label_pheatmap)]
p1<-pheatmap(trans_dist,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F)
# calculate the cosine simility within group and between groups;
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
if(!is.null(median(between_group))){
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame("1",log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)
}
# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+guides(color='none')+
              geom_density(alpha=.25) + theme_classic()+ scale_fill_manual(values=color_for_group)+ scale_color_manual(values=color_for_group)
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type",width=0.6,) +stat_compare_means()+guides(color = "none")+ scale_color_manual(values=color_for_group)
#raw counts heatmap 
# Heat map of expression  value 
    Idents(obj)<-obj$Annotation
    DefaultAssay(obj)<-"raw_RNA"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
    obj_data<-obj_data[rownames(barcode_label),]
    p4<-pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "#CC0000"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
    top_right<-plot_grid(p2,p3,labels = c(" "," "),rel_widths = c(2, 1))
    right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" "," "))
    last<-plot_grid(p1$gtable, right, labels = c(' ', ''), label_size = 12, ncol = 2)
    title <- ggdraw() + 
      draw_label(
        paste("combination","log2FC=",log2FC),
        fontface = 'bold',
        x = 0,
        hjust = 0
      )
    add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )

pdf("./00_Figure/Fig4/Fig4CD-C1234_trans_dist-within-between.pdf",width=10,height=5)
p2|p3
dev.off()

pdf("./00_Figure/FigS4/FigS4A-C1234_recluster_trans_dist_heatmap.pdf",width=6,height=5)
p1
dev.off()

obj$Annotation<- factor(obj$Annotation,levels=c("C1","C2","C3","C4"))

library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
# Fig4H coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  Idents(obj)<-obj$Annotation
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp######
  idents.plot <- Idents(obj)
  # plot region 
  start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
  ranges.show <- paste(seq,start,end,sep="-")
  col<- color_for_cluster
  p1<-CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 200,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 20,
    links=F
  )

pdf("./00_Figure/Fig4/Fig4H-C1234_recluster_trackplot.pdf",width=10,height=10)
p1
dev.off()

C1_barcode<-rownames(obj@meta.data[obj$group_manully=="C1",])
C2_barcode<-rownames(obj@meta.data[obj$group_manully=="C2",])
C3_barcode<-rownames(obj@meta.data[obj$group_manully=="C3",])
C4_barcode<-rownames(obj@meta.data[obj$group_manully=="C4",])

barcode<- data.frame(barcode=C4_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig4/RNA/C4barcode.txt",row.name=F,col.names=F)

sed -i 's/"//g' *

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes C4barcode.txt --cores 40 --out-bam barcode_RNA.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_RNA.bam  > ./barcode_RNA.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_RNA.bedGraph barcode_RNA.norm.bedGraph &> barcode_RNA.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_RNA.norm.bedGraph > barcode_RNA.norm.sorted.bedGraph
  bedGraphToBigWig barcode_RNA.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt C4_barcode_RNA.norm.bw









