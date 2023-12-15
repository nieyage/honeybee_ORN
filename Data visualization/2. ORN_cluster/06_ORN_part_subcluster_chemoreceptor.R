
#part1_subcluster <- subset(onecluster,idents=part1)
#part2_subcluster <- subset(onecluster,idents=part2)
#part3_subcluster <- subset(onecluster,idents=part3)
#part4_subcluster <- subset(onecluster,idents=part4)
#part1_subcluster
#part2_subcluster
#part3_subcluster
#part4_subcluster
##RNA analysis
#DefaultAssay(part1_subcluster) <- "integratedRNA_onecluster"
#part1_subcluster <- RunPCA(part1_subcluster,features= all_receptor_gene )
#pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/ElbowPlot_all_receptor_gene.pdf")
#ElbowPlot(part1_subcluster,ndims = 50, reduction = "pca")
#dev.off()
##build a tSNE visualization
#part1_subcluster <- RunTSNE(
#  object = part1_subcluster,
#  assay = "integratedRNA_onecluster",
#  min.dist = 0.001,
#  verbose = TRUE,
#  check_duplicates = FALSE,
#  reduction.name = "tsne.rna",
#  reduction.key = "rnatSNE_",
#  dims = 1:50
#)
#part1_subcluster <- FindNeighbors(object = part1_subcluster, reduction = 'pca', dims = 1:50)
#part1_subcluster <- FindClusters( object = part1_subcluster, verbose = FALSE, resolution =4,algorithm = 3)
#table(part1_subcluster$seurat_clusters)
#pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_ORN_cluster_WNN_part1_subcluster_all_receptor_gene.pdf",width=6,height=5)
#DimPlot(part1_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
#DimPlot(part1_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
#dev.off()
#DefaultAssay(part1_subcluster) <- "raw_RNA"
#Idents(part1_subcluster)<-part1_subcluster$seurat_clusters
#pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_raw_OR_dotplot_rawRNA_part1_subcluster_all_receptor_gene.pdf",width=30, height=8)
#p<-DotPlot(part1_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
#p
#dev.off()
#
#DefaultAssay(part2_subcluster) <- "integratedRNA_onecluster"
#part2_subcluster <- RunPCA(part2_subcluster,features= all_receptor_gene )
#pdf("./05_ORN_cluster/02_second_cluster/02_part2_subcluster/ElbowPlot_all_receptor_gene.pdf")
#ElbowPlot(part2_subcluster,ndims = 50, reduction = "pca")
#dev.off()
##build a tSNE visualization
#part2_subcluster <- RunTSNE(
#  object = part2_subcluster,
#  assay = "integratedRNA_onecluster",
#  min.dist = 0.001,
#  verbose = TRUE,
#  check_duplicates = FALSE,
#  reduction.name = "tsne.rna",
#  reduction.key = "rnatSNE_",
#  dims = 1:25
#)
#part2_subcluster <- FindNeighbors(object = part2_subcluster, reduction = 'pca', dims = 1:25)
#part2_subcluster <- FindClusters( object = part2_subcluster, verbose = FALSE, resolution =4,algorithm = 3)
#table(part2_subcluster$seurat_clusters)
#
#pdf("./05_ORN_cluster/02_second_cluster/02_part2_subcluster/part2_subcluster_all_receptor_gene.pdf",width=6,height=5)
#DimPlot(part2_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
#DimPlot(part2_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
#dev.off()
#DefaultAssay(part2_subcluster) <- "raw_RNA"
#Idents(part2_subcluster)<-part2_subcluster$seurat_clusters
#pdf("./05_ORN_cluster/02_second_cluster/02_part2_subcluster/part2_subcluster_OR_rawRNA_dotplot_all_receptor_gene.pdf",width=30, height=8)
#p<-DotPlot(part2_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
#p
#dev.off()
#
#DefaultAssay(part3_subcluster) <- "integratedRNA_onecluster"
#part3_subcluster <- RunPCA(part3_subcluster,features= all_receptor_gene )
#pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/ElbowPlot_all_receptor_gene.pdf")
#ElbowPlot(part3_subcluster,ndims = 50, reduction = "pca")
#dev.off()
##build a tSNE visualization
#part3_subcluster <- RunTSNE(
#  object = part3_subcluster,
#  assay = "integratedRNA_onecluster",
#  min.dist = 0.001,
#  verbose = TRUE,
#  check_duplicates = FALSE,
#  reduction.name = "tsne.rna",
#  reduction.key = "rnatSNE_",
#  dims = 1:40
#)
#part3_subcluster <- FindNeighbors(object = part3_subcluster, reduction = 'pca', dims = 1:40)
#part3_subcluster <- FindClusters( object = part3_subcluster, verbose = FALSE, resolution =4,algorithm = 3)
#table(part3_subcluster$seurat_clusters)
#
#pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_all_receptor_gene.pdf",width=6,height=5)
#DimPlot(part3_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
#DimPlot(part3_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
#dev.off()
#DefaultAssay(part3_subcluster) <- "raw_RNA"
#Idents(part3_subcluster)<-part3_subcluster$seurat_clusters
#pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_OR_rawRNA_dotplot_all_receptor_gene.pdf",width=30, height=8)
#p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
#p
#dev.off()
#p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
#dotplot_data<-p$data;
#dotplot_data$state<-"No";
#for(i in 1:nrow(dotplot_data)){
#   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
#   if(dotplot_data[i,]$pct.exp>30&&dotplot_data[i,]$avg.exp.scaled > 1){dotplot_data[i,]$state="Yes"};
#   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
# }
#dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
#dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
#cluster_info<-as.data.frame(table(dotplot_data$id))
#all_cluster<-cluster_info$Var1;
#multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
#one_classes <- as.character(cluster_info[cluster_info$Freq==1,1])
#
## distinguish_multi_OR
#library(pheatmap)
#library(dittoSeq)
#library(cowplot)
#library(lsa)
#library(ggpubr)
#DefaultAssay(part3_subcluster)<- "integratedRNA_onecluster"
#log2FCdata<-data.frame()
#pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_distinguish_multi_OR_all_receptor_gene.pdf",width=14,height=6)
#for (cluster in multiple_classes){
#print(cluster)
#obj<-subset(part3_subcluster,idents=cluster);
#obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
## add max_exp OR label for each cell
#ORN_count<-obj@assays$raw_RNA
#ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
#ORN_matrix<-as.matrix(ORN_count)
#ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
#ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
#barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
#  for (i in 1:length(colnames(ORN_matrix))){
#    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
#    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
#}
#  }
#barcode_label<-barcode_label[barcode_label$label!="NA",]
#barcode_label<-barcode_label[order(barcode_label$label),]
###cell cosine simility heatmap 
#embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
#embeddings <- embeddings[barcode_label$barcode,]
#trans_dist <- 1-cosine(t(embeddings))
#barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
#rownames(barcode_label_pheatmap)<-barcode_label$barcode
#col<-myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
#names(col)<-unique(barcode_label_pheatmap$OR)
#ann_colors= list(OR = col)
#p1<-pheatmap(trans_dist,
#         cluster_cols = TRUE,
#         cluster_rows = TRUE,
#         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
#         annotation_col = barcode_label_pheatmap,
#         annotation_colors = ann_colors,
#         annotation_row = barcode_label_pheatmap,
#         annotation_legend = TRUE,
#         show_rownames=F,
#         show_colnames=F
#  )
##calculate the cosine simility within group and between groups;
#rownames(barcode_label)<-barcode_label$barcode
#within_group<-c()
#between_group<-c()
#for (i in 1:nrow(trans_dist)){
#  for (j in 1:ncol(trans_dist)){
#    if(i!=j){
#      if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
#        within_group<-c(within_group,trans_dist[i,j]);
#      }
#      else{between_group<-c(between_group,trans_dist[i,j])}
#    }
#  }
#}
## calculate the FC 
#log2FC<-log2(median(between_group))-log2(median(within_group));
#test<-wilcox.test(within_group,between_group);
#pvalue<-test$p.value;
#data_subset<-data.frame(cluster,log2FC,pvalue)
#log2FCdata<-rbind(log2FCdata,data_subset)
## plot density line 
## manage data
#type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
#var<-c(within_group,between_group)
#data<-data.frame(type,var)
#data$type<-factor(data$type,levels=c("within-OR","between-OR"))
#p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
#              geom_density(alpha=.25) + theme_classic() 
## t-test
#p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
##method = "t.test"
## plot OR heatmap 
#DefaultAssay(obj)<-"raw_RNA";
#obj<-ScaleData(obj,features=all_receptor_gene);
##p4<-dittoHeatmap(obj,obj_features,slot ="scale.data",cluster_cols=T,scaled.to.max = TRUE)
#top_right<-plot_grid(p2,p3,labels = c("B","C"),ncol = 1)
##right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
#last<-plot_grid(p1$gtable, top_right, labels = c('A', ''),rel_widths = c(1, 1),label_size = 12, ncol = 2)
#title <- ggdraw() + 
#  draw_label(
#    paste("Cluster",cluster,":","log2FC=",log2FC),
#    fontface = 'bold',
#    x = 0,
#    hjust = 0
#  )
#add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
#print(add_title)
#}
#dev.off()
#
#> log2FCdata
#   cluster     log2FC        pvalue
#1        0 0.40654417  0.000000e+00
#2        1 0.32655189 2.389511e-107
#3        2 0.29165070 2.871736e-205
#4        3 0.34572612 6.044935e-203
#5        4 0.16297276  5.222954e-18
#6        5 0.55178488 8.962644e-279
#7        6 0.49250457  5.299571e-51
#8        7 0.42744528 9.762644e-220
#9        8 0.56813982 1.194533e-256
#10       9 0.48361880  0.000000e+00
#11      10 0.31425618  8.004350e-25
#12      15 0.17125858  1.523529e-34
#13      17 0.17421430  2.988831e-21
#14      19 0.81674666  0.000000e+00
#15      21 0.30859834  2.458680e-13
#16      24 0.29765000  3.398723e-33
#17      26 0.93242431  6.118863e-58
#18      27 0.03305013  4.159962e-01
#DefaultAssay(part4_subcluster) <- "integratedRNA_onecluster"
#part4_subcluster <- RunPCA(part4_subcluster,features= all_receptor_gene )
#pdf("./05_ORN_cluster/02_second_cluster/04_part4_subcluster/ElbowPlot_all_receptor_gene.pdf")
#ElbowPlot(part4_subcluster,ndims = 50, reduction = "pca")
#dev.off()
##build a tSNE visualization
#part4_subcluster <- RunTSNE(
#  object = part4_subcluster,
#  assay = "integratedRNA_onecluster",
#  min.dist = 0.001,
#  verbose = TRUE,
#  check_duplicates = FALSE,
#  reduction.name = "tsne.rna",
#  reduction.key = "rnatSNE_",
#  dims = 1:30
#)
#part4_subcluster <- FindNeighbors(object = part4_subcluster, reduction = 'pca', dims = 1:30)
#part4_subcluster <- FindClusters( object = part4_subcluster, verbose = FALSE, resolution =4,algorithm = 3)
#table(part4_subcluster$seurat_clusters)
#
#pdf("./05_ORN_cluster/02_second_cluster/04_part4_subcluster/part4_subcluster_all_receptor_gene.pdf",width=6,height=5)
#DimPlot(part4_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
#DimPlot(part4_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
#dev.off()
#DefaultAssay(part4_subcluster) <- "raw_RNA"
#Idents(part4_subcluster)<-part4_subcluster$seurat_clusters
#pdf("./05_ORN_cluster/02_second_cluster/04_part4_subcluster/part4_subcluster_OR_rawRNA_dotplot_all_receptor_gene.pdf",width=30, height=8)
#p<-DotPlot(part4_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
#p
#dev.off()#