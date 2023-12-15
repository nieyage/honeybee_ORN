library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(RColorBrewer)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
set.seed(1234)

honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
# Fig1B: 
pdf("./00_Figure/Fig1/Fig1B-honeybee_cluster_WNN_top.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(honeybee,cols=myUmapcolors,reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(honeybee,cols=myUmapcolors,reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(honeybee,cols=myUmapcolors,reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off();

pdf("./00_Figure/Fig1/Fig1B-honeybee_cluster_WNN_bottom.pdf",width=15,height=5)
###sample
p1 <- DimPlot(honeybee, cols=c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]), reduction = "umap.rna", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(honeybee, cols=c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]), reduction = "umap.atac",group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(honeybee, cols=c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]), reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))pdf("./00_Figure/Fig2/Fig2A-Unsupervised_ORN_cluster_WNN.pdf",width=9,height=6)
DimPlot(onecluster, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = F, label.size = 5, repel = TRUE)
dev.off();

dev.off()



# Fig1C: 
Idents(honeybee)<-honeybee$Annotation;
pdf("./00_Figure/Fig1/Fig1C-honeybee_annotation_allcelltype_WNN-Annotation_first.pdf",width=7,height=5)
DimPlot(honeybee, label = T, repel = TRUE, cols=c("#8ECC92","#49A5CC",myUmapcolors[6:10]), reduction = "wnn.umap",group.by = "Annotation")+ ggtitle("")
DimPlot(honeybee, label = F, repel = TRUE, cols=c("#8ECC92","#49A5CC",myUmapcolors[6:10]), reduction = "wnn.umap",group.by = "Annotation")+ ggtitle("")
dev.off()

# Fig1D:
# Major celltype track and violin plot 
##Track for Marker genes promoters
# remove unannotated 
honeybee<- subset(honeybee,idents=setdiff(levels(honeybee),"Unannotated"))

DefaultAssay(honeybee) <- "peaks"
# first compute the GC content for each peak
honeybee <- RegionStats(honeybee, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(honeybee)$tx_id <-Annotation(honeybee)$gene_name
#features<-c("Or2","LOC411079","LOC410151","LOC406073","LOC409780","Obp5","Obp11","Obp4")
# link peaks to genes
honeybee <- LinkPeaks(
  object = honeybee,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Syt1","LOC411079","LOC410151","5-ht7","Obp4","Obp5","LOC406073")
)
######Visulize track and RNA exp######
idents.plot <- Idents(honeybee)
# Neuron:
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group1-5132500-5137500",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group1-5125000-5137500"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Syt1",
  assay = "RNA"
)
p1<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Epithelial cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group5-6387000-6396000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group5-6387000-6396000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC411079",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

# glial cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group15-2484000-2489000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group15-2484000-2489000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC410151",
  assay = "RNA"
)
p3<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Sheath cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group6-730000-740000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group6-730000-740000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "5-ht7",
  assay = "RNA"
)
p4<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Obp4
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group9-11980000-11983000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group9-11980000-11983000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Obp4",
  assay = "SCT"
)
p5<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# obp5
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group9-11944000-11946000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group9-11944000-11946000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Obp5",
  assay = "RNA"
)
p6<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

set<- c("#8ECC92","#49A5CC",myUmapcolors[6:10])
p1<-p1& scale_fill_manual(values=set)&labs(title="Syt1")
p2<-p2& scale_fill_manual(values=set)&labs(title="LOC411079(GRH)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)&labs(title="LOC410151(repo)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p4<-p4& scale_fill_manual(values=set)&labs(title="5-ht7") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p5<-p5& scale_fill_manual(values=set)&labs(title="Obp4") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p6<-p6& scale_fill_manual(values=set)&labs(title="Obp5") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())

pdf("./00_Figure/Fig1/Fig1D-Marker_gene-select-peaktrack-WNN.pdf",height=8,width=24) 
p1|p2|p3|p4|p5|p6
dev.off()


######Visulize track and RNA exp######
idents.plot <- Idents(honeybee)
# Neuron:
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group1-5132500-5137500",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group1-5125000-5137500"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Syt1",
  assay = "RNA"
)
p1<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Epithelial cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group5-6387000-6396000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group5-6387000-6396000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC411079",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

# glial cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group15-2484000-2489000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group15-2484000-2489000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC410151",
  assay = "RNA"
)
p3<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Sheath cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group5-3430000-3440000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group5-3430000-3440000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC406073",
  assay = "RNA"
)
p4<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Obp4
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group9-11980000-11983000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group9-11980000-11983000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Obp4",
  assay = "SCT"
)
p5<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# obp5
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group9-11944000-11946000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group9-11944000-11946000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Obp5",
  assay = "RNA"
)
p6<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

set<- c("#8ECC92","#49A5CC",myUmapcolors[6:10])
p1<-p1& scale_fill_manual(values=set)&labs(title="Syt1")
p2<-p2& scale_fill_manual(values=set)&labs(title="LOC411079(GRH)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)&labs(title="LOC410151(repo)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p4<-p4& scale_fill_manual(values=set)&labs(title="PROS") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p5<-p5& scale_fill_manual(values=set)&labs(title="Obp4") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p6<-p6& scale_fill_manual(values=set)&labs(title="Obp5") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())

pdf("./00_Figure/Fig1/Fig1D-PROS-Marker_gene-select-peaktrack-WNN.pdf",height=8,width=24) 
p1|p2|p3|p4|p5|p6
dev.off()








# Fig1E:
Neuron<-readRDS("./03_Neuron/WNN_Neuron_integrated.rds")
markers <- FindAllMarkers(Neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Neuron_top<-markers[markers$cluster=="Orco-Neuron",7]
# Or2 and LOC408554
DefaultAssay(Neuron)<-"RNA"
color_for_neuron<-c("#6580A8","#63A8A6","#A8A7A7")
# ORN: 
# Non-ORN: 
pdf('./00_Figure/Fig1/Fig1E-Neuron_marker_FeaturePlot_WNN.pdf', width=13, height=4)
p1<- FeaturePlot(Neuron,reduction = 'wnn.umap',max.cutoff = 10,features = c("LOC413063") ,order=TRUE, ncol = 1)+ggtitle("LOC413063 (pepple)")
p2<- FeaturePlot(Neuron,cols =c("lightgrey", color_for_neuron[1]), reduction = 'wnn.umap',max.cutoff = 7,features = c("LOC551837") ,order=FALSE, ncol = 1)+ggtitle("LOC551837 (bgm)")
p3<- FeaturePlot(Neuron,cols =c("lightgrey", color_for_neuron[2]), reduction = 'wnn.umap',max.cutoff = 10,features = c("Or2") ,order=TRUE, ncol = 1)+ggtitle("Or2 (Orco)")
p1|p2|p3
dev.off()

# Fig1F:
DefaultAssay(Neuron) <- "peaks"
# first compute the GC content for each peak
Neuron <- RegionStats(Neuron, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(Neuron)$tx_id <- Annotation(Neuron)$gene_name

######Visulize track and RNA exp######
idents.plot <- Idents(Neuron)
# Neuron
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group5-12940000-12944000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group5-12940000-12944000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "LOC413063",
  assay = "RNA"
)
p1<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# LOC551837
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group1-22596500-22598000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group1-22596500-22598000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "LOC551837",
  assay = "RNA"
)
p3<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Orco
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group1-5722000-5725000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group1-5722000-5725000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "Or2",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

# Orco
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group1-5722000-5725000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group1-5722000-5725000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "Or2",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

set<- color_for_neuron[1:2]
p1<-p1& scale_fill_manual(values=set)&labs(title="LOC413063 (pepple)")
p2<-p2& scale_fill_manual(values=set)&labs(title="Or2 (Orco)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)&labs(title="LOC551837 (bgm)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())

pdf("./00_Figure/Fig1/Fig1F-Neuron-Orco-track.pdf",width=13,height=4)
p1|p3|p2
dev.off()


# dotplot in Orco-/+ neuron and non-neuron cluster 
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
Neuron<-readRDS("./03_Neuron/WNN_Neuron_integrated.rds")
Orco_p_barcode<- rownames(Neuron@meta.data[Neuron$Annotation=="Orco+Neuron",])
Orco_n_barcode<- rownames(Neuron@meta.data[Neuron$Annotation=="Orco-Neuron",])

type<- rep("Non-neuron",length(honeybee$Annotation))
for (i in 1:nrow(honeybee@meta.data)){
  if(colnames(honeybee)[i] %in% Orco_p_barcode){type[i]="Orco+Neuron"}
  if(colnames(honeybee)[i] %in% Orco_n_barcode){type[i]="Orco-Neuron"}
}
table(type)
honeybee$Annotation_type<- type
Idents(honeybee)<- honeybee$Annotation_type
markers <- FindAllMarkers(honeybee, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Orco_n_Neuron_top<-markers[markers$cluster=="Orco-Neuron",7]
Orco_p_Neuron_top<-markers[markers$cluster=="Orco+Neuron",7]
Non_Neuron_top<-markers[markers$cluster=="Non-Neuron",7]


Neuron_marker <- c("LOC726238","LOC724243","LOC100820634","LOC726770")
#Neuron_marker <- intersect(Orco_n_Neuron_top,Orco_p_Neuron_top)
Orco_marker <- c("Or2","LOC552552","LOC726019","LOC551704")
Orco_n_marker <- Orco_n_Neuron_top[1:4]

Idents(honeybee)<- factor(Idents(honeybee),levels=c("Orco-Neuron","Orco+Neuron","Non-neuron"))
DefaultAssay(honeybee)<-"SCT"
pdf("./00_Figure/Fig1/Fig1F-b-dotplot.pdf",width=8,height=4)
DotPlot(honeybee,features = c(Neuron_marker,Orco_marker,Orco_n_marker)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dev.off()



#Fig1G: proportion of cell: ORN vs Non ORN;
ORN<-c(2405,2032,2840)
all<-c(2773,2555,3540)
Non_ORN<- all-ORN;
data<-data.frame(group=c(rep("ORN",3),rep("Non-ORN",3)),
  Sample=rep(c("NE","Nurse","Forager"),2),
  cellnumber=c(ORN,Non_ORN))
data$group<-factor(data$group,levels=c("ORN","Non-ORN"))
data$Sample<-factor(data$Sample,levels=c("NE","Nurse","Forager"))
pdf("./00_Figure/Fig1/Fig1G-ORNvsNonORN_proportion.pdf",width=4,height=4)
ggplot(data = data, aes_string(x = "group", y = "cellnumber", 
        fill = "Sample")) +  xlab("orig.ident") + ylab("Percent of cells") + 
        scale_fill_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2])) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_bw();
ggplot(data = data, aes_string(x = "group", y = "cellnumber", 
        fill = "Sample")) +  xlab("orig.ident") + ylab("cell number") + 
        scale_fill_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2])) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()
dev.off()


#cross species
#Fig1H:
#OR vs Non-ORN porportion 
# change to 3 types
# Orco+Neuron+;Orco-Neuron+;Orco-Neuron-;
Apis_mellifera<-c(7277,55,1536);
#Fly_Joint.integrated <- readRDS("/md01/liyh526/project/Fly-cooperate/shiny_app/WNN_fly_integrated_annotation_antenna.rds")
#Fly_ORN<-readRDS("/md01/liyh526/project/Fly-cooperate/7.18run/WNN_ORN_integrated_antenna.rds")
#mosquito_all<-readRDS("/data/R02/nieyg/project/honeybee/data/publish_data/mosquito/SeuratObject1_Antenna_mergedBatches_AllCells.rds")
#mosquito_neuron<-readRDS("/data/R02/nieyg/project/honeybee/data/publish_data/mosquito/SeuratObject2_Antenna_mergedBatches_Neurons.rds")
#Drosophila_melanogaster<-c(2808,154,3756) #our data 


Drosophila_melanogaster<-c(7940,4266,25048)
Aedes_aegypti<-c(4742,433,8704);
cross_species_cellnumber<-data.frame(species=c(rep("Apis mellifera",3),rep("D.melanogaster",3),rep("Ae.aegypti",3)),
  celltype=rep(c("Orco+Neuron","Orco-Neuron","Non-neuron"),3),
  cellnumber=c(Apis_mellifera,Drosophila_melanogaster,Aedes_aegypti))
cross_species_cellnumber$celltype<-factor(cross_species_cellnumber$celltype,levels=c("Orco+Neuron","Orco-Neuron","Non-neuron"));
cross_species_cellnumber$species<-factor(cross_species_cellnumber$species,levels=c("Apis mellifera","D.melanogaster","Ae.aegypti"))
pdf("./00_Figure/Fig1/Fig1H-cross_species_ORNvsNonORN_proportion_3types_publishdata.pdf",width=4,height=4)
ggplot(data = cross_species_cellnumber, aes_string(x = "species", y = "cellnumber", 
        fill = "celltype")) +  xlab(" ") + ylab("% Percent of cells") + 
        scale_fill_manual(values = c(color_for_neuron[1:2],"#CACCCC")) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_bw()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
dev.off();

# Fig1I:
# OR,IR,GR gene number
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="IR",]$gene_name)
Apis_mellifera<-c(length(OR_gene),length(IR_gene),length(GR_gene));
Drosophila_melanogaster<-c(60,66,60)
Aedes_aegypti<-c(114,135,107)
cross_species_genenumber<-data.frame(species=c(rep("Apis_mellifera",3),rep("Drosophila_melanogaster",3),rep("Aedes_aegypti",3)),
  genetype=rep(c("ORs","IRs","GRs"),3),
  genenumber=c(Apis_mellifera,Drosophila_melanogaster,Aedes_aegypti))
cross_species_genenumber$genetype<-factor(cross_species_genenumber$genetype,levels=c("ORs","IRs","GRs"));
cross_species_genenumber$species<- factor(cross_species_genenumber$species,levels=c("Apis_mellifera","Drosophila_melanogaster","Aedes_aegypti"))
pdf("./00_Figure/Fig1/Fig1I-cross_species_ORGRIRgene_proportion.pdf",width=4,height=4)
p<-ggplot(data = cross_species_genenumber, aes_string(x = "species", y = "genenumber", 
        fill = "genetype")) +  xlab(" ") + ylab("chemosensory receptor genes") + 
        scale_fill_manual(values = c("#B781CC","#7A9BCC","#CC6E7B")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = genenumber), size = 3, hjust = 0.5, vjust = 3, position = "stack") 
dev.off();

#Fig1J:
# the OB barplot
cross_species_glomeruli<-data.frame(species=c("Apis mellifera","D. melanogaster","Ae. aegypti"),
  glomeruli=c(160,55,65));
cross_species_glomeruli$species<- factor(cross_species_glomeruli$species,levels=c("Apis mellifera","D. melanogaster","Ae. aegypti"))
pdf("./00_Figure/Fig1/Fig1J-cross_species_glomeruli.pdf",width=3,height=4)
p<-ggplot(data = cross_species_glomeruli, aes_string(x = "species", y = "glomeruli", 
        fill = "species")) +  xlab(" ") + ylab("# of glomeruli") + 
        geom_bar( stat = "identity",width=0.6,color = 'black', fill='grey') +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = glomeruli), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

# version 2023.10.7

# Fig1H: "Orco+Neuron","Orco-Neuron","Non-neuron"

ORN<-c(2405,2032,2840)
other_neuron<- c(30,18,7)
all<-c(2773,2555,3540)
Non_ORN<- all-ORN-other_neuron
color_for_neuron<-c("#6580A8","#63A8A6","#A8A7A7")

data<-data.frame(celltype=c(rep("Orco+Neuron",3),rep("Orco-Neuron",3),rep("Non_Neuron",3)),
  Stage=rep(c("NE","Nurse","Forager"),3),
  cellnumber=c(ORN,other_neuron,Non_ORN))
data$celltype<-factor(data$celltype,levels=c("Orco+Neuron","Orco-Neuron","Non_Neuron"))
data$Stage<-factor(data$Stage,levels=c("NE","Nurse","Forager"))

pdf("./00_Figure/Fig1/Fig1H-ORNvsNonORN_proportion_among_stages.pdf",width=4,height=4)
ggplot(data = data, aes_string(x = "Stage", y = "cellnumber", 
        fill = "celltype")) +  xlab(" ") + ylab("% Percent of cells") + 
        scale_fill_manual(values = c(color_for_neuron[1:2],"#CACCCC")) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
dev.off();

