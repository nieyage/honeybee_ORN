# combination OR recluster by ATAC+RNA signal 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
DefaultAssay(ORN)<-"raw_RNA"
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
   Idents(ORN)<- ORN$subcluster
   obj<-subset(ORN,idents="p1:14");
   obj_features<- unique(dotplot_data[dotplot_data$id%in%"p1:14",]$features.plot)

    DefaultAssay(obj)<-"SCT"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
    obj_data<-obj_data[order(obj_data$Or25,obj_data$Or26,obj_data$Or27,decreasing=T),]
# 定义平滑函数
smooth_column <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 2:(length(col) - 1)) {
    smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
  }
  smoothed[1] <- (col[1] + col[2]) / 2
  smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
  return(smoothed)
}

# 对每一列进行平滑
smoothed_data <- as.data.frame(lapply(obj_data, smooth_column))

pdf("./00_Figure/Fig4/Fig4A-OR25_27_heatmap-smooth.pdf",width=8,height=2)
    pheatmap(t(smoothed_data),
             cluster_cols = F,
             cluster_rows = F,
             border_color=NA,
             color = colorRampPalette(c("white", "#CA0002"))(100),
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()


   obj<-subset(ORN,idents="p4:5");
   obj_features<- unique(dotplot_data[dotplot_data$id%in%"p4:5",]$features.plot)

    DefaultAssay(obj)<-"SCT"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
colnames(obj_data)<-c("Or40","Or39")
obj_data<- obj_data[,c(2,1)]
obj_data<-obj_data[order(obj_data$Or39,obj_data$Or40,decreasing=T),]
smoothed_data <- as.data.frame(lapply(obj_data, smooth_column))

pdf("./00_Figure/Fig4/Fig4A-Or39_40_heatmap_smooth.pdf",width=8,height=2)
    pheatmap(t(smoothed_data),
             cluster_cols = F,
             cluster_rows = F,
             border_color=NA,
             color = colorRampPalette(c("white", "#CA0002"))(100),
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()

