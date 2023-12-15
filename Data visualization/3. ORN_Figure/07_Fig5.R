# Fig5A plot the tree and label 39-42:

colors_for_exp_pattern<- c("#476D87","#E95C59")
obj<- ORN
DefaultAssay(obj)<-"integratedRNA_onecluster"
embeddings <- Embeddings(object = obj, reduction = "obj_features_pca")[,1:50]
data.dims <- lapply(X = levels(x = obj), FUN = function(x) {
    cells <- WhichCells(object = obj, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = obj)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))

library(ggtree)
pdf("./00_Figure/Fig5/Fig5A-ORN-tree-cosine.pdf",width=5,height=14)
ggtree(data.tree,layout="circular") + 
geom_tiplab()+ 
geom_hilight(node=c(68,98,99),fill = "blue",alpha = 0.6)
dev.off()


obj<- readRDS("./00_Figure/Fig4/Fig4-last-data-obj.rds")
Idents(obj)<-obj$group_manully
obj_features<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
barcode_label<-data.frame(barcode=colnames(obj),label=Idents(obj))
DefaultAssay(obj)<-"raw_RNA"
obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
barcode_label<- barcode_label[order(barcode_label$label),]
obj_data<- obj_data[rownames(barcode_label),]


C39_data<-obj_data[rownames(barcode_label[barcode_label$label=="C1",]),]
C40_data<-obj_data[rownames(barcode_label[barcode_label$label=="C2",]),]
C41_data<-obj_data[rownames(barcode_label[barcode_label$label=="C3",]),]
C42_data<-obj_data[rownames(barcode_label[barcode_label$label=="C4",]),]

C39_data<- C39_data[which(rowSums(C39_data)>0),]
C40_data<- C40_data[which(rowSums(C40_data)>0),]
C41_data<- C41_data[which(rowSums(C41_data)>0),]
C42_data<- C42_data[which(rowSums(C42_data)>0),]

C42_data<- C42_data[order(C42_data$LOC100578045),]
C41_data<- C41_data[order(C41_data$LOC100578045,C41_data$LOC107963999),]
C39_data<- C39_data[order(C39_data$LOC100578045,C39_data$LOC107963999,C39_data$LOC410603,C39_data$`Or63-b`),]
C40_data<- C40_data[order(C40_data$LOC100578045,C40_data$LOC107963999,C40_data$LOC410603,C40_data$`Or63-b`),]
last_data_heatmap<- rbind(C39_data,C40_data,C41_data,C42_data)

write.csv(last_data_heatmap,"00_Figure/Fig4/Fig4E_4gene_heatmap_data.csv")

last_data_heatmap<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig4/Fig4E_4gene_heatmap_data.csv",row.names=1)
smooth_column <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 2:(length(col) - 1)) {
    smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
  }
  smoothed[1] <- (col[1] + col[2]) / 2
  smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
  return(smoothed)
}

smooth_column_5cell <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 3:(length(col) - 2)) {
    smoothed[i] <- (col[i - 2] + col[i - 1] + col[i] + col[i + 1] + col[i + 2]) / 5
  }
  smoothed[1] <- (col[1] + col[2] + col[3]) / 3
  smoothed[2] <- (col[1] + col[2] + col[3] + col[4]) / 4
  smoothed[length(col) - 1] <- (col[length(col) - 3] + col[length(col) - 2] + col[length(col) - 1] + col[length(col)]) / 4
  smoothed[length(col)] <- (col[length(col) - 2] + col[length(col) - 1] + col[length(col)]) / 3
  return(smoothed)
}

library(pheatmap)


# 对每一列进行平滑
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[6:8])
barcode_label_pheatmap<-data.frame(label=c(rep("C1",nrow(C39_data)),rep("C2",nrow(C40_data)),rep("C3",nrow(C41_data)),rep("C4",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-rownames(last_data_heatmap)
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)
last_data_heatmap_smoothed_data <- as.data.frame(lapply(last_data_heatmap, smooth_column))
rownames(last_data_heatmap_smoothed_data)<- rownames(last_data_heatmap)


pdf("./00_Figure/Fig4/Fig4E_4gene_heatmap.pdf",height=3,width=8)
pheatmap(t(last_data_heatmap),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "#CC0000"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "#CC0000"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()

# 
# clusterMatrix <- function(input_matrix) {
#   # Define the clustering method and other parameters
#   clustering_method <- "complete"  # You can change this to other methods like "ward.D", "single", etc.
#   # Perform clustering
#   p <- pheatmap(
#     input_matrix,
#     clustering_method = clustering_method,
#     cluster_cols = F,
#     cluster_rows = T,
#   )
#   clustered_matrix <- input_matrix[p$tree_row$order,]
#   # Return the clustered matrix
#   return(clustered_matrix)
# }
# smooth_column <- function(col) {
#   smoothed <- numeric(length(col))
#   for (i in 2:(length(col) - 1)) {
#     smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
#   }
#   smoothed[1] <- (col[1] + col[2]) / 2
#   smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
#   return(smoothed)
# }
# 
# C39_data_clustered<- clusterMatrix(C39_data)
# C39_data_clustered_smoothed_data <- as.data.frame(lapply(C39_data_clustered, smooth_column))
# rownames(C39_data_clustered_smoothed_data)<- rownames(C39_data_clustered)
# C40_data_clustered<- clusterMatrix(C40_data)
# C40_data_clustered_smoothed_data <- as.data.frame(lapply(C40_data_clustered, smooth_column))
# rownames(C40_data_clustered_smoothed_data)<- rownames(C40_data_clustered)
# C41_data_clustered<- clusterMatrix(C41_data)
# C41_data_clustered_smoothed_data <- as.data.frame(lapply(C41_data_clustered, smooth_column))
# rownames(C41_data_clustered_smoothed_data)<- rownames(C41_data_clustered)
# C42_data_clustered<- clusterMatrix(C42_data)
# C42_data_clustered_smoothed_data <- as.data.frame(lapply(C42_data_clustered, smooth_column))
# rownames(C42_data_clustered_smoothed_data)<- rownames(C42_data_clustered)

# clustered_smoothed_data<- rbind(C39_data_clustered_smoothed_data,
# 	C40_data_clustered_smoothed_data,
# 	C41_data_clustered_smoothed_data,
# 	C42_data_clustered_smoothed_data)
# 


colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/Fig4/Fig4E_4gene_heatmap_pink.pdf",height=3,width=8)
pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colors_list[10:100],
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colors_list[5:100],
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colors_list[1:100],
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()




# Fig5 UCSC track:
# LG2 gene exp dotplot 
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/Unsupervised_ORN_cluster_WNN_add_subcluster.rds")
Idents(ORN)<- ORN$subcluster
DefaultAssay(ORN) <- "integratedRNA_onecluster"
object<- ORN
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
#cluster order by tree
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

ORN$subcluster<-factor(ORN$subcluster,levels=cluster_order)
Idents(ORN)<-factor(ORN$subcluster,levels=cluster_order)


G2<- chemoreceptor_info_data[chemoreceptor_info_data$seqnames=="Group2",]$gene_name
receptor.dot <- DotPlot(ORN, features = G2) + #scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
dotplot_data<-receptor.dot$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>5&&dotplot_data[i,]$avg.exp.scaled >= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
#dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% dotplot_data$features.plot])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]

dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))


DefaultAssay(ORN)<- "SCT"
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/G2_OR_UnSupervised_WNN_dotplot-signif-feature_orderbytree.pdf",width=15, height=20)
p<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p

dev.off()

data<- ORN@assays$RNA[c("Or9","Or10","Or11"),]
Or9 <- colnames(data[,which(data[1,]>0)])
Or10 <- colnames(data[,which(data[2,]>0)])
Or11 <- colnames(data[,which(data[3,]>0)])

Or9_10<- intersect(Or9,Or10)


library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN)<-"peaks"
clear_multiple_classes_order_id<- unique(dotplot_data$id)
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/OR_without_nopower_trackplot.pdf",width=10,height=10)
for (cluster in clear_multiple_classes_order_id){
print(cluster)
obj<-subset(ORN,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 100,
  annotation = TRUE,
  extend.downstream = 100,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
}
dev.off()


chrLG2:1056100-2513400


chrLG2:9976700-10160000

region1:chrLG2:1056100-1186600
library(ggplot2)
# 读取GTF文件，假设GTF文件名为"yourfile.gtf"
gtf_data <- read.csv("/Users/fraya/Documents/project/honeybee/chemoreceptor/OR_gene_df_gene_exp_pattern_Group2.csv", header = FALSE)

# 指定要绘制的区域范围，这里假设绘制chr1的1-10000范围内的基因
chr <- "Group2"
start_pos <- 1056100
end_pos <- 1186600
head(gtf_data)
## 筛选符合区域范围的基因
genes <- subset(gtf_data, V3 == chr & V4 >= start_pos & V5 <= end_pos)

ggplot(genes, aes(xmin = V4, xmax = V5, 
                  y = V3, fill = V7, label = V1, forward = 1)) +
  geom_gene_arrow() +
  facet_wrap(~ V3, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  geom_gene_label(align = "centre",min.size = 2,reflow=T)

genes


chr <- "Group2"
start_pos <- 9976700
end_pos <- 10160000
head(gtf_data)
## 筛选符合区域范围的基因
genes <- subset(gtf_data, V3 == chr & V4 >= start_pos & V5 <= end_pos)

ggplot(genes, aes(xmin = V4, xmax = V5, 
                  y = V3, fill = V7, label = V1, forward = 1)) +
  geom_gene_arrow() +
  facet_wrap(~ V3, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  geom_gene_label(align = "centre",min.size = 2,reflow=T)

genes


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
markers<- read.csv("./00_Figure/Fig4/Fig4-DEG_markers-DEG_heatmap_need2select_showmarkers.csv")



chemosensory_behavior<- c("Orco","Snmp1","a10","kuz")
phasphorous_metabolic<- c("Gapdh1","Ald","AcCoAS","COX7C","retm","kdn","Adk1","ND-ASHI","Cyt-c-p")
neurontransmission<- c("Ca-alpha1D","CG43066","slo","Gs2","Slob","Gclc","to","Eaat2","Cirl")
respose_stimulus<- c("Peblll","antdn","Drsl5")
biological_regulation<- c("SelR","CG4928","CG1358","Iscu")

nervous_system_development<- c("fra","RhoGEF64C","dlp","fz2","Con","Fas3","otk","Drl-2","Nrt","sls","miple1","Ace","stan","nerfin-1","robo1","robo2","Sema2a","IA-2","N")
axon_development<- c("Btk","N","cib","brat","Antp","stg","kuz","ft","Myo81F","Act42A","Act87E","mam","Rbfox1",#"A2bp1",
  "ps")
developmental_process<-c("sano","nw","shep","CG2852","Hr4","ftz-f1","rn","ru","krn")

fly_gene<- list(chemosensory_behavior,phasphorous_metabolic,neurontransmission,respose_stimulus,biological_regulation,
  nervous_system_development,axon_development,developmental_process)

fly2honeybee<- read.csv("./10_fly_honeybee_ID_trans/05_fly2honeybee.csv")

for(i in 1:8){
  gene<- fly_gene[[i]]
  tmp<-gene[which(gene%in% fly2honeybee$fly_gene)]
  print(gene)
  print(tmp)
}

grep("Csp",fly2honeybee$fly_gene)


chemosensory_behavior<- c("Or2","LOC413995","LOC725209",#Snmp1
  "CSP3","CSP4",#a10
  "LOC409611"# kuz
  )
phasphorous_metabolic<- c("LOC410122","LOC413924",#"Gapdh1",
  "LOC725455","LOC550785",#"Ald",
  "LOC409624",#"AcCoAS",
  "LOC726297",#"COX7C",
  "LOC551515",#"retm",
  "Pcl","LOC410059",#"kdn",
  "Adk1",
  "LOC551660",#"ND-ASHI",
  "LOC724543","CytC"#"Cyt-c-p"
  )

neurontransmission<- c(
  "Tpcn1","LOC100578899",#"Ca-alpha1D",
  "LOC409142",#"CG43066",
  "LOC413994",#"slo",
  "GlnS",#"Gs2",
  "LOC409492",#"Slob",
  "LOC100578132",#"Gclc",
  "LOC552773",#"to",
  "Eaat-2","LOC408769",#Eaat2
  "LOC552142"#"Cirl"
  )

respose_stimulus<- c(
  "CSP6","CSP3","CSP4","CSP1",#"PebIII",
  "LOC406147"#"antdh",
)
biological_regulation<- c(
  "LOC724494",#"SelR",
  "LOC413134",#"CG4928",
  "LOC408828","LOC107965219",#"CG1358",
  "LOC409130"#"Iscu"
  )

nervous_system_development_honeybee<- fly2honeybee[which(nervous_system_development%in%fly2honeybee$fly_gene),]$honeybee_gene_name
axon_development_honeybee<- fly2honeybee[which(axon_development%in%fly2honeybee$fly_gene),]$honeybee_gene_name
developmental_process_honeybee<- fly2honeybee[which(developmental_process%in%fly2honeybee$fly_gene),]$honeybee_gene_name


LOC408769 Eaat2同源基因


gene<- "LOC408769"
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp#####
  idents.plot <- Idents(obj)
  # plot region 
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "SCT",
  genes.use = gene
)
obj$group_manully<- factor(obj$group_manully,levels=c("C1","C2","C3","C4"))
Idents(obj)<- obj$group_manully
pdf("./00_Figure/Fig4/Fig4-LOC408769-track.pdf",width=10,height=5)
for(i in gene){
    p1<-CoveragePlot(
    object = obj,
    region = i,
    features=i,
    window = 200,
    expression.assay = "SCT",
    expression.slot = "data",
    extend.upstream = 1500,
    annotation = TRUE,
    links = TRUE,
    extend.downstream = 1500
  )
print(p1)
}

# recall peak
DefaultAssay(obj)<-"ATAC"
peak2<-CallPeaks(
       obj,
       group.by = "group_manully",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       #broad = FALSE,
       format = "BED",
       broad=TRUE,
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="./00_Figure/Fig4/",
       combine.peaks=TRUE
)

macs2_counts <- FeatureMatrix(
     fragments = Fragments(obj),
     features = peak2,
     cells = colnames(obj)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# create a new assay using the MACS2 peak set and add it to the Seurat object
obj[["peaks_obj2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(obj),
  annotation = Annotation(obj)
)

library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
gene="LOC409565"
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp#####
  idents.plot <- Idents(obj)
  # plot region 
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "RNA",
  method = "spearman", # pearson
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.5,
  score_cutoff = 0.04,
  genes.use = gene
)
obj$group_manully<- factor(obj$group_manully,levels=c("C1","C2","C3","C4"))
Idents(obj)<- obj$group_manully
pdf("./00_Figure/Fig4/Fig4-LOC409565-track-recallpeak.pdf",width=10,height=5)
for(i in gene){
    p1<-CoveragePlot(
    object = obj,
    region = i,
    features=i,
    window = 200,
    expression.assay = "RNA",
    expression.slot = "data",
    extend.upstream = 1500,
    annotation = TRUE,
    links = TRUE,
    extend.downstream = 1500
  )
print(p1)
}
dev.off()

dev.off()
DefaultAssay(obj)<-"SCT"
pdf("./00_Figure/Fig4/Fig4-LOC408769-boxplot.pdf",width=5,height=5)
VlnPlot(obj,"LOC408769",pt.size=0)+ 
 geom_boxplot(width=.2,col="black",fill="white")+  
 NoLegend()
dev.off()



Neuronpeptide<- c("sNPF-R") # 
Ion_channels<- c("Ca-alpha1T","Ca-beta","Cngl","GluClalpha","Hk","Ih","KCNQ","para","Piezo","ppk9","ppk13","ppk16","ppk22","ppk23","ppk28","Shab","Shaw","SK","slo","Teh2","tipE")
Hormone_receptor<- c("Dh31-R","Dh44-R2","Dh44","ETHR","hec")
Neurontransmitter_receptors<- c("5-HT2B","Dop1R1","DopEcR","GABA-B-R1","Oamb","Octbeta1R","Octbeta2R","Octbeta3R")
Neuronpeptide_receptors<- c("AstC-R2","CCKLR-17D1","CG13995","Lgr3","MsR1","NPFR","Pdfr","SPR","TkR99D","TrissinR")



fly_gene<- list(Neuronpeptide,Ion_channels,Hormone_receptor,Neurontransmitter_receptors,Neuronpeptide_receptors)

fly2honeybee<- read.csv("./10_fly_honeybee_ID_trans/05_fly2honeybee.csv")

for(i in 1:5){
  gene<- fly_gene[[i]]
  tmp<-gene[which(gene%in% fly2honeybee$fly_gene)]
  print(gene)
  print(tmp)
}

gene_honeybee<- fly2honeybee[which(Neuronpeptide_receptors%in%fly2honeybee$fly_gene),]$honeybee_gene_name
gene_honeybee %in% markers$gene


axon_development_honeybee<- fly2honeybee[which(axon_development%in%fly2honeybee$fly_gene),]$honeybee_gene_name
developmental_process_honeybee<- fly2honeybee[which(developmental_process%in%fly2honeybee$fly_gene),]$honeybee_gene_name



library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
gene=c("LOC413399","Foxp","LOC410246","LOC409565")
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp#####
  idents.plot <- Idents(obj)
  # plot region 

obj$group_manully<- factor(obj$group_manully,levels=c("C1","C2","C3","C4"))
Idents(obj)<- obj$group_manully
pdf("./00_Figure/Fig4/FigS4-C1234-DEG-track-recallpeak.pdf",width=10,height=5)
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "RNA",
  method = "spearman", # pearson
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.1,
  score_cutoff = 0.15,
  genes.use = gene[1]
)
    p1<-CoveragePlot(
    object = obj,
    region = gene[1],
    features=gene[1],
    window = 200,
    expression.assay = "RNA",
    expression.slot = "data",
    extend.upstream = 1500,
    annotation = TRUE,
    links = TRUE,
    extend.downstream = 1500
  )
print(p1)
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "RNA",
  method = "spearman", # pearson
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.1,
  score_cutoff = 0.2,
  genes.use = gene[2]
)
    p1<-CoveragePlot(
    object = obj,
    region = gene[2],
    features=gene[2],
    window = 200,
    expression.assay = "RNA",
    expression.slot = "data",
    extend.upstream = 1500,
    annotation = TRUE,
    links = TRUE,
    extend.downstream = 1500
  )
print(p1)
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "RNA",
  method = "spearman", # pearson
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.1,
  score_cutoff = 0.1,
  genes.use = gene[3]
)
    p1<-CoveragePlot(
    object = obj,
    region = gene[3],
    features=gene[3],
    window = 200,
    expression.assay = "RNA",
    expression.slot = "data",
    extend.upstream = 1500,
    annotation = TRUE,
    links = TRUE,
    extend.downstream = 1500
  )
print(p1)
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "RNA",
  method = "spearman", # pearson
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.1,
  score_cutoff = 0.1,
  genes.use = gene[4]
)
    p1<-CoveragePlot(
    object = obj,
    region = gene[4],
    features=gene[4],
    window = 200,
    expression.assay = "RNA",
    expression.slot = "data",
    extend.upstream = 1500,
    annotation = TRUE,
    links = TRUE,
    extend.downstream = 1500
  )
print(p1)
dev.off()


# MP OR pair vs single exp OR pair 
MP_gene1<- c("LOC102656904","LOC100577101","LOC107965761","LOC725052")
MP_gene2<- c("LOC102656221","Or41","LOC102655285","LOC100576839")
MP_OR_pair<- data.frame(gene1=MP_gene1,gene2=MP_gene2)
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
singleOR_cluster<-as.character(cluster_info[cluster_info$Freq==1,1])
multiOR<- unique(dotplot_data[which(dotplot_data$id %in% multiOR_cluster),]$features.plot)
singleOR<- unique(dotplot_data[which(dotplot_data$id %in% singleOR_cluster),]$features.plot)

singleOR<- colnames(OR_kmer_cor)[which(singleOR%in% colnames(OR_kmer_cor))]
remaining_gene<-singleOR
Single_OR_pair<- data.frame()
for(gene1 in remaining_gene){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    tmp<- data.frame(gene1=gene1,gene2=gene2)
    Single_OR_pair<- rbind(Single_OR_pair,tmp)
  }
}

# promoter similarity 
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_promoter_fasta<-readAAStringSet("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_location/OR_promoter_merged.fa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(OR_promoter_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
tree<-as.phylo(clust)
# the similarity of OR pair in the same cluster and random OR pair 
# All OR pair similarity matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100

MP_pair_seq_data<- data.frame()
for (i in 1:nrow(MP_OR_pair)){
  gene1<- MP_OR_pair[i,1]
  gene2<- MP_OR_pair[i,2]
  tmp_data<- data.frame(gene1=gene1,gene2=gene2,promoter_sequence_similarity=dist_percent[gene1,gene2],type="MP_coexp")
  MP_pair_seq_data<- rbind(MP_pair_seq_data,tmp_data)
}
Single_pair_seq_data<- data.frame()
for (i in 1:nrow(Single_OR_pair)){
  gene1<- Single_OR_pair[i,1]
  gene2<- Single_OR_pair[i,2]
  tmp_data<- data.frame(gene1=gene1,gene2=gene2,promoter_sequence_similarity=dist_percent[gene1,gene2],type="Single_exp")
  Single_pair_seq_data<- rbind(Single_pair_seq_data,tmp_data)
}
last_pair_seq_data<- rbind(Single_pair_seq_data,MP_pair_seq_data)
library(ggpubr)
colors_for_exp_pattern<- c("#476D87","#E95C59")
last_pair_seq_data$type<-factor(last_pair_seq_data$type,levels=c("Single_exp","MP_coexp"))
pdf("./00_Figure/Fig5/Fig5I-MPvsSingle_sequence_similarity.pdf",width=3,height=3)
ggboxplot(last_pair_seq_data, x="type", y="promoter_sequence_similarity", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("sequence similarity of OR pair promoter")+
scale_color_manual(values = colors_for_exp_pattern)
dev.off()


# kmer correlation
library(stringr)
a=read.csv('/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/matrix.count',header=F,sep="\t")
colnames(a)=c('OR','kmer','count')
library(reshape2)
counts=dcast(a,formula=kmer~OR)
counts[is.na(counts)]=0
rownames(counts)<- counts[,1]
counts=counts[,-1]
# k=5 
k_5_data <- counts[which(str_length(rownames(counts))==5),]
# correlation coefficient 
OR_kmer_cor <- cor(k_5_data,method="pearson")
MP_pair_kmer_cor_data<- data.frame()
for (i in 1:nrow(MP_OR_pair)){
  gene1<- MP_OR_pair[i,1]
  gene2<- MP_OR_pair[i,2]
  tmp_data<- data.frame(gene1=gene1,gene2=gene2,kmer_cor=OR_kmer_cor[gene1,gene2],type="MP_coexp")
  MP_pair_kmer_cor_data<- rbind(MP_pair_kmer_cor_data,tmp_data)
}
Single_pair_kmer_cor_data<- data.frame()
for (i in 1:nrow(Single_OR_pair)){
  gene1<- Single_OR_pair[i,1]
  gene2<- Single_OR_pair[i,2]
  tmp_data<- data.frame(gene1=gene1,gene2=gene2,kmer_cor=OR_kmer_cor[gene1,gene2],type="Single_exp")
  Single_pair_kmer_cor_data<- rbind(Single_pair_kmer_cor_data,tmp_data)
}
last_pair_kmer_cor_data<- rbind(Single_pair_kmer_cor_data,MP_pair_kmer_cor_data)
library(ggpubr)
colors_for_exp_pattern<- c("#476D87","#E95C59")
last_pair_kmer_cor_data$type<-factor(last_pair_kmer_cor_data$type,levels=c("Single_exp","MP_coexp"))
pdf("./00_Figure/Fig5/Fig5J-MPvsSingle_kmer_cor.pdf",width=3,height=3)
ggboxplot(last_pair_kmer_cor_data, x="type", y="kmer_cor", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("kmer correlation of OR pair promoter")+
scale_color_manual(values = colors_for_exp_pattern)
dev.off()

kaks_data<- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/result_dir_muscle2/OR_pair_kaks.out")
colnames(kaks_data)<- c("OR_pair","method","KA","Ks","Ka/Ks","p_value","length","S-sites","N-sites","Folo-sites","Substitutions","S-sub","N-sub","Folo-S-sub","Folo-N-sub","Divergence-Time","Substitution-Rate-Ratio","GC","ML_score","AICc","Akaike-Weight","Model")
dN_data <- kaks_data[,c(1,3)]
dS_data <- kaks_data[,c(1,4)]
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
# corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)
sequence_dist<- as.matrix(sdist)
sequence_data<- melt(sequence_dist);
sequence_data$OR_pair<- paste(sequence_data$Var1,sequence_data$Var2,sep="-");
sequence_data$rate <- 100-(sequence_data$value/411)*100
sequence_data<- sequence_data[sequence_data$rate<100,]
OR_pair_data <- sequence_data[sequence_data$rate>40,1:2]
OR_pair_data_last<- rbind(OR_pair_data,OR_pair_data[,c(2,1)])
OR_pair_data_last$last<- paste(OR_pair_data_last$Var1,OR_pair_data_last$Var2,sep="-")


# label MP point in dN of CDS and kmer_cor of promoter 
 promoter_matrix<- melt(OR_kmer_cor);
 promoter_matrix$OR_pair<- paste(promoter_matrix$Var1,promoter_matrix$Var2,sep="-");
 data<- merge(dN_data,promoter_matrix,by="OR_pair")
 data<- data[,c("OR_pair","KA","value")]
 colnames(data)<- c("OR_pair","CDS_dN","promoter_kmer")
 data<- data[!duplicated(data),]
 
data<- data[data$OR_pair%in% OR_pair_data_last$last,]
data$type<-"other"
data$type[which(data$OR_pair%in%c("LOC102656904-LOC102656221","LOC102656221-LOC102656904",
  "LOC100577101-Or41","Or41-LOC100577101",
  "LOC107965761-LOC102655285", "LOC102655285-LOC107965761",
  "LOC100576839-LOC725052","LOC725052-LOC100576839"))]="MP_OR_pair"
library(ggpmisc)
pdf("./00_Figure/Fig5/Fig5G-OR_CDS_dN_vs_promoter_kmer_label_MP.pdf",width=6,height=8)
p1<- ggplot(data, aes(CDS_dN,promoter_kmer,col=type))+ 
       geom_point()+ 
       scale_color_manual(values=c("red","black"))+
       geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
       theme_classic()+ggtitle("The correlation between kmer of promoter and dS of CDS")+
       stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)

p2<- ggplot(data, aes(CDS_dN,promoter_kmer))+ 
       geom_point()+ 
       geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
       theme_classic()+ggtitle("The correlation between kmer of promoter and dN of CDS")+
       stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
p1/p2
dev.off()


# label MP point in dS of CDS and kmer_cor of promoter 
  data<- merge(dS_data,promoter_matrix,by="OR_pair")
 data<- data[,c("OR_pair","Ks","value")]
 colnames(data)<- c("OR_pair","CDS_dS","promoter_kmer")
 data<- data[!duplicated(data),]
data<- data[data$OR_pair%in% OR_pair_data_last$last,]

Single_MP_ORpair<- c(paste(MP_OR_pair$gene1,MP_OR_pair$gene2,sep="-"),paste(Single_OR_pair$gene1,Single_OR_pair$gene2,sep="-"))
data<- data[data$OR_pair%in% Single_MP_ORpair,]

data$type<-"other"
data$type[which(data$OR_pair%in%c("LOC102656904-LOC102656221","LOC102656221-LOC102656904",
  "LOC100577101-Or41","Or41-LOC100577101",
  "LOC107965761-LOC102655285", "LOC102655285-LOC107965761",
  "LOC100576839-LOC725052","LOC725052-LOC100576839"))]="MP_OR_pair"
library(ggpmisc)
pdf("./00_Figure/Fig5/Fig5H-OR_CDS_dS_vs_promoter_kmer_label_MP.pdf",width=6,height=4)
p<- ggplot(data, aes(CDS_dS,promoter_kmer,col=type))+ 
       geom_point()+ 
       scale_color_manual(values=c("red","black"))+
       geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
       theme_classic()+ggtitle("The correlation between kmer of promoter and dS of CDS")+
       stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
print(p)
p<- ggplot(data, aes(CDS_dS,promoter_kmer))+ 
       geom_point()+ 
       geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
       theme_classic()+ggtitle("The correlation between kmer of promoter and dS of CDS")+
       stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
print(p)
dev.off()











