library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(patchwork)
set.seed(1234)
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
# version 2023.10.7
# FigS1B: QC for diff stage 
library(scCustomize)
Idents(honeybee)<- honeybee$orig.ident
pdf('./00_Figure/FigS1/FigS1A-QC for diff stage .pdf', width=6, height=16)
p1<- VlnPlot(honeybee, features = "nCount_RNA",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p2<- VlnPlot(honeybee, features = "nFeature_RNA",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p3<- VlnPlot(honeybee, features = "nCount_ATAC",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p4<- VlnPlot(honeybee, features = "nFeature_ATAC",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p5<- VlnPlot(honeybee, features = "percent.mt",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p6<- VlnPlot(honeybee, features = "nucleosome_signal",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p7<- VlnPlot(honeybee, features = "TSS.enrichment",pt.size = 0,cols = c("#5CC8F2","#009E73","#E69F00"))+geom_boxplot(width=.2,col="black")+ theme(axis.title.x = element_blank())
p1/p2/p3/p4/p5/p6/p7
dev.off()

# FigS1A:
 DefaultAssay(honeybee)<-"RNA"
  pdf('./00_Figure/FigS1A-allcelltype_marker_FeaturePlot_WNN.pdf', width=14, height=8)
  p1<-FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 10,features = c("Syt1") ,order=TRUE, ncol = 1)&labs(title="Neuron: Syt1")
  p3<-FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 10,features = c("LOC411079") ,order=TRUE, ncol = 1)&labs(title="Epithelial cell:LOC411079 (GRH)")
  p4<-FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 10,features = c("LOC410151") ,order=TRUE, ncol = 1)&labs(title="Glial cell:LOC410151 (repo)")
  p5<-FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 10,features = c("LOC406073") ,order=TRUE, ncol = 1)&labs(title="Sheath cell:LOC406073 (PROS)")
  p7<-FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 30,features = c("Obp4") ,order=TRUE, ncol = 1)
  p8<-FeaturePlot(honeybee, reduction = 'wnn.umap',max.cutoff = 50,features = c("Obp5") ,order=TRUE, ncol = 1)
  f1<-p1|p3|p4
  f2<-p5|p7|p8;
  f1/f2
  dev.off()

# FigS1B:
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(honeybee) <- "peaks"
Idents(honeybee)<- honeybee$Annotation
# first compute the GC content for each peak
honeybee <- RegionStats(honeybee, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(honeybee)$tx_id <-Annotation(honeybee)$gene_name
######Visulize track and RNA exp######
# Obp4
idents.plot <- Idents(honeybee)
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group9-11981000-11982000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group9-11981000-11982000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Obp4",
  assay = "RNA"
)
p1<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# obp5
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group9-11944500-11945500",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group9-11944500-11945500"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Obp5",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#B3B3B3")
p1<-p1& scale_fill_manual(values=set)
p2<-p2& scale_fill_manual(values=set)& theme(strip.text.y.left = element_blank(),strip.background = element_blank())
pdf("./00_Figure/FigS1B-Obp-peaktrack-WNN.pdf",height=8,width=6) 
p1|p2
dev.off()

# FigS1C
DefaultAssay(honeybee)<-"RNA"
markers <- FindAllMarkers(honeybee, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers<-markers[which(markers$p_val_adj<0.05),]
table(markers$cluster);
write.csv(markers,"./00_Figure/FigS1C-Allcelltype_top_markers_gene_description.csv");

# downsample Neuron to 400
Neuron_barcode<-colnames(Neuron);
down_neuron<- sample(Neuron_barcode,400)
Otherbarcode<-setdiff(colnames(honeybee),Neuron_barcode);

down_obj<-subset(honeybee,cell=c(down_neuron,Otherbarcode))
markers <- FindAllMarkers(down_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);

down_obj<-ScaleData(down_obj,features=rownames(down_obj));
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

pdf("./00_Figure/FigS1/FigS1G-Allcelltype_top_markers-DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = down_obj,features=top10$gene[1:70],label=T, group.colors =c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494","#B3B3B3" ),
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = down_obj,features=top10$gene[1:70],label=T, group.colors =c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494","#B3B3B3" ),
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])

dev.off();

# FigS1D:
markers <- FindAllMarkers(honeybee, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Epithelial<-markers[markers$cluster=="Epithelial cell",7]
glial<-markers[markers$cluster=="Glial cell",7]
Sheath<-markers[markers$cluster=="Sheath cell",7]
Neuron<-markers[markers$cluster=="Neuron",7]
OBP5<-markers[markers$cluster=="Obp5+support cell",7]
OBP4<-markers[markers$cluster=="Obp4+support cell",7]

library(AnnotationHub)
library(biomaRt)
library(dplyr)
library(goseq)
library(DOSE)
library(stringr)
library(clusterProfiler)
library(GenomicRanges)
library(AnnotationDbi)


Apis_mellifera.OrgDb <-loadDb("/md01/nieyg/project/honeybee/antenna_gtf/Apis_mellifera.OrgDb")
columns(Apis_mellifera.OrgDb)
#GO for marker genes in OBP5 cells
  gene_unique<-OBP5
  id <- mapIds(x = Apis_mellifera.OrgDb,
               keys = gene_unique,
               keytype = "SYMBOL",
               column = "ENTREZID");
  id<-na.omit(id) 
  OBP5_ego<-enrichGO(gene     = id,
                OrgDb    = Apis_mellifera.OrgDb,
                ont      = "BP",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                readable = TRUE,
                pool = FALSE)
  OBP5_ego<-OBP5_ego[OBP5_ego$pvalue<0.05,]

  #GO for marker genes in OBP4 cells
  gene_unique<-OBP4
  id <- mapIds(x = Apis_mellifera.OrgDb,
               keys = gene_unique,
               keytype = "SYMBOL",
               column = "ENTREZID");
  id<-na.omit(id) 
  OBP4_ego<-enrichGO(gene     = id,
                OrgDb    = Apis_mellifera.OrgDb,
                ont      = "BP",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                readable = TRUE,
                pool = FALSE)
  OBP4_ego<-OBP4_ego[OBP4_ego$pvalue<0.05,]
#GO for marker genes in Epithelial cells
  gene_unique<-Epithelial
  id <- mapIds(x = Apis_mellifera.OrgDb,
               keys = gene_unique,
               keytype = "SYMBOL",
               column = "ENTREZID");
  id<-na.omit(id) 
  #ggo <- groupGO(gene     = id,
  #               OrgDb    = Apis_mellifera.OrgDb,
  #               ont      = "BP",
  #               level    = 3,
  #               readable = TRUE)
  #ggo<-ggo[order(ggo$Count,decreasing = T),]
  Epithelial_ego<-enrichGO(gene     = id,
                OrgDb    = Apis_mellifera.OrgDb,
                ont      = "BP",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                readable = TRUE,
                pool = FALSE)
  Epithelial_ego<-Epithelial_ego[Epithelial_ego$pvalue<0.05,]
#GO for marker genes in glial cells
  gene_unique<-glial
  id <- mapIds(x = Apis_mellifera.OrgDb,
               keys = gene_unique,
               keytype = "SYMBOL",
               column = "ENTREZID");
  id<-na.omit(id) 
  glial_ego<-enrichGO(gene     = id,
                OrgDb    = Apis_mellifera.OrgDb,
                ont      = "BP",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                readable = TRUE,
                pool = FALSE)
  glial_ego<-glial_ego[glial_ego$pvalue<0.05,]
#GO for marker genes in Neuron cells
  gene_unique<- Neuron
  id <- mapIds(x = Apis_mellifera.OrgDb,
               keys = gene_unique,
               keytype = "SYMBOL",
               column = "ENTREZID");
  id<-na.omit(id) 
  Neuron_ego<-enrichGO(gene     = id,
                OrgDb    = Apis_mellifera.OrgDb,
                ont      = "BP",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                readable = TRUE,
                pool = FALSE)
  Neuron_ego<-Neuron_ego[Neuron_ego$pvalue<0.05,]
#GO for marker genes in Sheath cells
  gene_unique<-unique(gsub("-g.*","",Sheath));
  id <- mapIds(x = Apis_mellifera.OrgDb,
               keys = gene_unique,
               keytype = "SYMBOL",
               column = "ENTREZID");
  id<-na.omit(id) 
  Sheath_ego<-enrichGO(gene     = id,
                OrgDb    = Apis_mellifera.OrgDb,
                ont      = "BP",
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                readable = TRUE,
                pool = FALSE)
  Sheath_ego<-Sheath_ego[Sheath_ego$pvalue<0.05,]
Epithelial_ego$celltype="Epithelial cell";
glial_ego$celltype="Glial cell"
Sheath_ego$celltype="Sheath cell"
Neuron_ego$celltype="Neuron"
OBP5_ego$celltype="Obp5+support cell"
OBP4_ego$celltype="Obp4+support cell"
all_ego<-rbind(Neuron_ego,Epithelial_ego,glial_ego,Sheath_ego,OBP4_ego,OBP5_ego)
write.csv(all_ego,"./00_Figure/FS1D-Allcelltype_GO.csv")

all_ego<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS1/FS1D-Allcelltype_GO-20231113.csv")
top5_term <- all_ego 
library(ggplot2)
top5_term$celltype<-factor(top5_term$celltype,levels=c("Neuron","Epithelial cell","Glial cell","Sheath cell",
  "Obp4+support cell","Obp5+support cell"))
pdf("./00_Figure/FigS1/FigS1D-Allcelltype_GO.pdf",width=10,height=12)
p <- ggplot(top5_term,aes(y=Count,x=Description,fill=pvalue)) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(celltype~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 14),
            legend.position="right",
            legend.title = element_text(size=18),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
p
dev.off()

# FigS1F:Orco+Neuron and Orco-Neuron:
Neuron<-readRDS("./03_Neuron/WNN_Neuron_integrated.rds")
Idents(Neuron) <- Neuron$Annotation
pdf("./00_Figure/FigS1E-Neuron_cluster_WNN.pdf",width=15,height=5)
p1 <- DimPlot(Neuron, reduction = "wnn.umap", group.by = "seurat_clusters",label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("cluster")
p2 <- DimPlot(Neuron, cols=c("#F0A04B", "#183A1D"),reduction = "wnn.umap", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("annotation")
p3 <- DimPlot(Neuron, cols=c("#5CC8F2","#009E73","#E69F00"),reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("stage")
p1 +p2+p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

