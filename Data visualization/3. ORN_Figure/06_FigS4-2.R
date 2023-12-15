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
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));


obj<- readRDS("./00_Figure/Fig4/Fig4-last-data-obj.rds")
Idents(obj)<-obj$group_manully
obj_features<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")


# Find the DEG among the 4 cluster:
DefaultAssay(obj)<- "SCT"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
markers<- markers[!duplicated(rownames(markers)),]
write.csv(markers,"./00_Figure/Fig4/Fig4-DEG_markers-DEG_heatmap_need2select_showmarkers.csv")

# tau index find DEG 

# 2. tau cluster specific:
library(VGAM)
ORN_matrix<-as.matrix(GetAssayData(obj));
#filter the gene only appear in a few cells 
cell_pct = function(data){
    pct<-length(data[data!=0])/length(data)
    return(pct)
}
gene_pct<-apply(ORN_matrix,1,cell_pct)
gene_pass_pct<-names(gene_pct[gene_pct>0.005])
obj<-NormalizeData(obj)
obj$group_manully<- factor(obj$group_manully,levels=levels(obj))
ORN_avg<-AverageExpression(
       obj,
       assays = "raw_RNA",
       features = gene_pass_pct,
       return.seurat = FALSE,
       group.by = "group_manully",
       #add.ident = NULL,
       slot = "data")
ORN_avg<-ORN_avg$raw_RNA
colnames(ORN_avg)<- levels(obj)
#https://fmicompbio.github.io/swissknife/reference/specificityScore-methods.html#value-1
source("/md01/nieyg/project/honeybee/add_antenna/swissknife-master/R/tissue_specificity_score.R")
library(matrixStats)
gene_tau<-specificityScore(
  ORN_avg,
  method = c("tau", "TSI", "counts"),
  #group = ORN$subcluster,
  thresh = 0,
  expr_values = "logcounts",
  na.rm = FALSE
)
names(gene_tau)<-rownames(ORN_avg)
# plot the density plot for gene_tau 
pdf("./00_Figure/Fig4/DEG_tau_density.pdf",width=10,height=5)
data<- as.data.frame(gene_tau)
ggplot(data, aes(x=data[,1])) + xlab("")+
              geom_density(alpha=.25) + theme_classic() 

#Kmeans
dev.off()

#plot the cluster specific pheatmap
gene_specific<-names(which(gene_tau>=1))
gene_specific_data<-as.data.frame(ORN_avg[gene_specific,])
data<-data.frame()
for (gene in gene_specific){
    gene_cluster_avg<-gene_specific_data[gene,]
    gene_specific_cluster<-names(gene_cluster_avg[which(gene_cluster_avg==max(gene_cluster_avg))])
    data_subset<-data.frame(gene,cluster=gene_specific_cluster,tau=gene_tau[gene]);
    data<-rbind(data,data_subset)
}
data$cluster<-factor(data$cluster,levels=colnames(gene_specific_data))
data<-data[order(data$cluster),]
tau1_gene<-data$gene
write.csv(data,"./00_Figure/Fig4/DEG_tau.csv")
C2_C3<- data[which(data$cluster%in%c("C2","C3")),]$gene
tau_Avg <-ORN_avg[tau1_gene,]
library(pheatmap)
count=t(scale(t(tau_Avg),scale = T,center = T))
pdf("./00_Figure/Fig4/DEG_tau_density.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
dev.off();

pdf("./00_Figure/Fig4/Fig4-DEG_tau-DEG_heatmap.pdf",width=5,height=5)
DoHeatmap(object = obj,features=C2_C3,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = obj,features=C2_C3,label=T, group.colors =color_for_cluster,
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

# downsample C1 to 50
C1_barcode<-rownames(obj@meta.data[obj$group_manully=="C1",])
C2_barcode<-rownames(obj@meta.data[obj$group_manully=="C2",])
C3_barcode<-rownames(obj@meta.data[obj$group_manully=="C3",])
C4_barcode<-rownames(obj@meta.data[obj$group_manully=="C4",])
down_C1<- sample(C1_barcode,50)


down_obj<-subset(obj,cell=c(down_C1,C2_barcode,C3_barcode,C4_barcode))
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
markers$cluster<- factor(markers$cluster,levels=c("C1","C2","C3","C4"))
#top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
markers<- markers[order(markers$cluster),]
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[7:9])

pdf("./00_Figure/Fig4/Fig4-DEG_markers-DEG_heatmap.pdf",width=5,height=4)
DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = obj,features=top10$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])

DoHeatmap(object = down_obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = down_obj,features=top10$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

# show select gene exp and ATAC signal
gene<- unique(markers$gene)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp#####
  idents.plot <- Idents(obj)
  # plot region 
pdf("./00_Figure/Fig4/Fig4-DEG-track.pdf",width=10,height=5)
for(i in gene){
    p1<-CoveragePlot(
    object = obj,
    region = i,
    features=i,
    window = 150,
    expression.assay = "raw_RNA",
    expression.slot = "data",
    extend.upstream = 500,
    annotation = TRUE,
    extend.downstream = 500
  )
   p2<-  p1 + scale_fill_manual(color_for_cluster)
print(p2)
}
dev.off()



# Find the DEP among the 4 cluster:
DefaultAssay(obj)<- "peaks_ORN_subcluster"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
write.csv(markers,"./00_Figure/Fig4/Fig4-DEG_markers-DEP_heatmap.csv")

sambaNight<- c("#1873CC", "#1798E5" ,"#00BFFF", "#4AC596" ,"#00CC00" ,"#A2E700", "#FFFF00" ,"#FFD200","#FFA500")
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

obj<- ScaleData(obj)
DefaultAssay(down_obj)<- "peaks_ORN_subcluster"
down_obj<- ScaleData(down_obj)
pdf("./00_Figure/Fig4/Fig4-DEP_markers-DEP_heatmap.pdf",width=7,height=5)
DoHeatmap(object = down_obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
DoHeatmap(object = down_obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])

DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
dev.off();


obj_features<- markers$gene
barcode_label<-data.frame(barcode=colnames(obj),label=Idents(obj))
DefaultAssay(obj)<-"raw_RNA"
obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
barcode_label<- barcode_label[order(barcode_label$label),]
obj_data<- obj_data[rownames(barcode_label),]
C39_data<-obj_data[rownames(barcode_label[barcode_label$label=="C1",]),]
C40_data<-obj_data[rownames(barcode_label[barcode_label$label=="C2",]),]
C41_data<-obj_data[rownames(barcode_label[barcode_label$label=="C3",]),]
C42_data<-obj_data[rownames(barcode_label[barcode_label$label=="C4",]),]

library(pheatmap)

clusterMatrix <- function(input_matrix) {
  # Define the clustering method and other parameters
  clustering_method <- "complete"  # You can change this to other methods like "ward.D", "single", etc.
  # Perform clustering
  p <- pheatmap(
    input_matrix,
    clustering_method = clustering_method,
    cluster_cols = F,
    cluster_rows = T,
  )
  clustered_matrix <- input_matrix[p$tree_row$order,]
  # Return the clustered matrix
  return(clustered_matrix)
}
smooth_column <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 2:(length(col) - 1)) {
    smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
  }
  smoothed[1] <- (col[1] + col[2]) / 2
  smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
  return(smoothed)
}

C39_data_clustered<- clusterMatrix(C39_data)
C39_data_clustered_smoothed_data <- as.data.frame(lapply(C39_data_clustered, smooth_column))
rownames(C39_data_clustered_smoothed_data)<- rownames(C39_data_clustered)

C40_data_clustered<- clusterMatrix(C40_data)
C40_data_clustered_smoothed_data <- as.data.frame(lapply(C40_data_clustered, smooth_column))
rownames(C40_data_clustered_smoothed_data)<- rownames(C40_data_clustered)
C41_data_clustered<- clusterMatrix(C41_data)
C41_data_clustered_smoothed_data <- as.data.frame(lapply(C41_data_clustered, smooth_column))
rownames(C41_data_clustered_smoothed_data)<- rownames(C41_data_clustered)
C42_data_clustered<- clusterMatrix(C42_data)
C42_data_clustered_smoothed_data <- as.data.frame(lapply(C42_data_clustered, smooth_column))
rownames(C42_data_clustered_smoothed_data)<- rownames(C42_data_clustered)


clustered_smoothed_data<- rbind(C39_data_clustered_smoothed_data,
	C40_data_clustered_smoothed_data,
	C41_data_clustered_smoothed_data,
	C42_data_clustered_smoothed_data)

# 对每一列进行平滑
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[6:8])
barcode_label_pheatmap<-data.frame(label=c(rep("C39",nrow(C39_data)),rep("C40",nrow(C40_data)),rep("C41",nrow(C41_data)),rep("C42",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-rownames(clustered_smoothed_data)
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)

smoothed_data <- as.data.frame(lapply(obj_data, smooth_column))
rownames(smoothed_data)<- rownames(obj_data)
barcode_label_pheatmap<-data.frame(label=c(rep("C39",nrow(C39_data)),rep("C40",nrow(C40_data)),rep("C41",nrow(C41_data)),rep("C42",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-rownames(smoothed_data)

clustered_smoothed_data<-  t(scale(t(clustered_smoothed_data),scale=T))

pdf("./00_Figure/FigS4/FigS4_DEG_smoothed_heatmap.pdf",height=8,width=8)
pheatmap(t(clustered_smoothed_data),
             cluster_cols = F,
             cluster_rows = TRUE,
             #color = colorRampPalette(c("white", "#CC0000"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()




# for C1: promoter1 
obj_C1<- subset(obj,idents="C1")
ORN_count<-obj_C1@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
library(UpSetR)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
library(ComplexUpset)
pdf("./00_Figure/Fig4/Fig4F-promoter1-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
#upset(data, sets=c("Or63_b","LOC410603","LOC107963999","LOC100578045"), order.by=c("degree","freq"),empty.intersections=TRUE,  mb.ratio = c(0.5, 0.5), keep.order = F)
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[3:4],obj_upset[3],obj_upset[4]#'Outside of known sets'
         ),
    queries=list(
        upset_query(intersect=obj_upset[1],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:2],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:3],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:4],color="#4BA9D1",fill="#4BA9D1")
    )
)
dev.off()

# four gene exp violin plot in C1 
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C1"])
data <- FetchData(obj,vars = obj_features,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter1-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()


# for C2: promoter2 
obj_C2<- subset(obj,idents="C2")
ORN_count<-obj_C2@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./00_Figure/Fig4/Fig4F-promoter2-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[2],  color=color_for_cluster[2],fill=color_for_cluster[2]),
        upset_query(intersect=obj_upset[2:3],color=color_for_cluster[2],fill=color_for_cluster[2]),
        upset_query(intersect=obj_upset[2:4],color=color_for_cluster[2],fill=color_for_cluster[2])
    ))
dev.off()


# four gene exp violin plot in C2
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C2"])
data <- FetchData(obj,vars = obj_features,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter2-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()



# for C3: promoter3 
obj_C3<- subset(obj,idents="C3")
ORN_count<-obj_C3@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./00_Figure/Fig4/Fig4F-promoter3-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[3],  color=color_for_cluster[3],fill=color_for_cluster[3]),
        #upset_query(intersect=obj_upset[3:3],color=color_for_cluster[3],fill=color_for_cluster[3]),
        upset_query(intersect=obj_upset[3:4],color=color_for_cluster[3],fill=color_for_cluster[3])
    ))
dev.off()



# four gene exp violin plot in C3
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C3"])
data <- FetchData(obj,vars = obj_features,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter3-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()


# for C4: promoter4
obj_C4<- subset(obj,idents="C4")
ORN_count<-obj_C4@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./00_Figure/Fig4/Fig4F-promoter4-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[4],  color=color_for_cluster[4],fill=color_for_cluster[4])#,
        #upset_query(intersect=obj_upset[3:3],color=color_for_cluster[3],fill=color_for_cluster[3]),
        #upset_query(intersect=obj_upset[3:4],color=color_for_cluster[3],fill=color_for_cluster[3])
    ))
dev.off()


# four gene exp violin plot in C3
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C4"])
data <- FetchData(obj,vars = obj_features,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter4-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()

