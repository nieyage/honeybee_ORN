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
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype_Neuron_subcluster.rds")

# find the DEG among all celltypes and get the specifical TF in honeybees
Idents(honeybee)<- honeybee$Annotation_subcluster
DefaultAssay(honeybee)<-"peaks"
da_peaks <- FindAllMarkers(
  object = honeybee,
  test.use = 'LR',
  logfc.threshold = 0.1,
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
#markers <- FindAllMarkers(honeybee, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
table(da_peaks$cluster)
da_peaks<-da_peaks[da_peaks$p_val_adj<0.05,]
write.csv(da_peaks,"./02_All_celltype/DEP_FindAllMarkers_all_celltype.csv")
#verification
da_peaks<- da_peaks[da_peaks$avg_log2FC>1,]
peak2show<- rownames(da_peaks)

honeybee<-ScaleData(honeybee,features=rownames(honeybee))
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
library(ArchR)
pdf("./02_All_celltype/DEP_FindAllMarkers_all_celltype_heatmap.pdf",width=20,height=10)
DoHeatmap(object = honeybee,features=peak2show,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = honeybee,features=peak2show,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = honeybee,features=peak2show,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = ArchRPalettes$solarExtra)+NoLegend()
dev.off();
# avg 
peak_Avg <- AverageExpression(honeybee,features=peak2show,assays = "peaks")
library(pheatmap)
count=t(scale(t(peak_Avg$peaks),scale = T,center = F))
pdf("./02_All_celltype/DEP_FindAllMarkers_all_celltype_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();

# motif enrichment 
# Get a list of motif position frequency matrices from the JASPAR database

library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'insects', all_versions = FALSE)
)
DefaultAssay(honeybee) <- 'peaks'
# add motif information
honeybee <- AddMotifs(
  object = honeybee,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = pfm
)

All_motif_info <- data.frame()
for (cluster in levels(honeybee)){
  cluster_peak <- da_peaks[da_peaks$cluster==cluster,]$gene;
  enriched.motifs <- FindMotifs(honeybee,features = cluster_peak);
  enriched.motifs$cluster <- cluster;
  All_motif_info <- rbind(All_motif_info,enriched.motifs)
}
library(dplyr)
All_motif_info<- All_motif_info[All_motif_info$p.adjust<0.05,]
FC1<- All_motif_info[All_motif_info$fold.enrichment>1,]
#top3 <- All_motif_info %>% group_by(cluster) %>% top_n(n = 3, wt = fold.enrichment)
motif2show<-unique(FC1$motif.name)

# motif enrich matrix
motif_matrix<-matrix(ncol=length(motif2show),nrow=length(levels(honeybee)))
colnames(motif_matrix)<-motif2show
rownames(motif_matrix)<-levels(honeybee)

last_motif_info<-All_motif_info[which(All_motif_info$motif.name%in%motif2show),]
for (i in 1:nrow(last_motif_info)){
  cluster=last_motif_info[i,]$cluster;
  motif=last_motif_info[i,]$motif.name;
  motif_matrix[cluster,motif]<-last_motif_info[i,]$fold.enrichment
}
motif_matrix[is.na(motif_matrix)]<- 0

library(pheatmap)
#count=t(scale(t(motif_matrix),scale = T,center = T))
pdf("./02_All_celltype/FindAllMarker_DEP_motif_heatmap.pdf",width=20,height=20)
pheatmap(motif_matrix,cluster_cols = T,cluster_rows = T,
              color = colorRampPalette(c("white", "firebrick3"))(100),
              cellwidth = 10, cellheight = 10,
              show_rownames=T,show_colnames=T)
dev.off()

# show ORN specific motif
#cluster_closest_open <- ClosestFeature(ORN, cluster_peak)
ORN_cluster_peak <- rownames(da_peaks[da_peaks$cluster=="Orco+Neuron",])
enriched.motifs <- FindMotifs(
  object = honeybee,
  features = ORN_cluster_peak
)
honeybee <- RunChromVAR(
  object = honeybee,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor
)
pdf("./02_All_celltype/ORN_cluster_Motif.pdf",width=12,height=8)
DefaultAssay(honeybee) <- 'peaks'
MotifPlot(
  object = honeybee,
  motifs = head(rownames(enriched.motifs))
)
dev.off()
pdf("./02_All_celltype/ORN_cluster_Motif_chromVar.pdf",width=7,height=8)
DefaultAssay(honeybee) <- 'chromvar'
FeaturePlot(
  object = honeybee,reduction = 'wnn.umap',
  features = head(rownames(enriched.motifs)),
  min.cutoff = 'q10',order=TRUE
  )
dev.off()

differential.activity <- FindAllMarkers(
  object = honeybee,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
pdf("./02_All_celltype/ORN_cluster_Motif_differential.activity.pdf",width=12,height=8)
MotifPlot(
  object = honeybee,
  motifs = differential.activity[differential.activity$cluster=="Orco+Neuron",]$gene,
  assay = 'peaks'
)
dev.off()
pdf("./02_All_celltype/ORN_cluster_Motif_differential.activity_featureplot.pdf",width=16,height=12)
DefaultAssay(honeybee) <- 'chromvar'
FeaturePlot(
  object = honeybee,reduction = 'wnn.umap',
  features = differential.activity[differential.activity$cluster=="Orco+Neuron",]$gene,
  min.cutoff = 'q10',order=TRUE
  )
dev.off()

motif_differential.activity<-  differential.activity[differential.activity$cluster=="Orco+Neuron",]$gene
motif_peak_enrich<- rownames(enriched.motifs)

overlap_motif<- intersect(motif_differential.activity,motif_peak_enrich)

# Motif footprinting
# Now we can footprint any motif that we have positional information for. By default, this includes every instance of the motif in the genome. We can instead use the in.peaks = TRUE parameter to include only those motifs that fall inside a peak in the assay. The Footprint() function gathers all the required data and stores it in the assay. 
# We can then plot the footprinted motifs using the PlotFootprint() function.
# gather the footprinting information for sets of motifs
# some peak length beyond the boundary 
gr <- granges(honeybee[['peaks']])
gr <- gr[seqnames(gr) == "Group8"]
# peak end
min(start(gr))
# genome end 
honeybee_genome<- BSgenome.Amel.HAv3.1.update.chemoreceptor
end(honeybee_genome$`Group8`)

DefaultAssay(honeybee) <- 'peaks'
honeybee <- Footprint(
  object = honeybee,
  motif.name =overlap_motif,     
  upstream = 1,
  downstream = 100,
  in.peaks=T,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor
)
 

# plot the footprint data for each group of cells
pdf("./02_All_celltype/overlap_motif_Footprint.pdf",width=12,height=12)
PlotFootprint(honeybee, features = overlap_motif)
dev.off()

# the expression of the overlap motif list 
TF_info<- last_motif_info[last_motif_info$ %in% overlap_motif,c()]
TF_name<- 
label<- 
DefaultAssay(honeybee)<-"raw_RNA"
pdf("./02_All_celltype/overlap_motif_TF-ORN-exp.pdf",width=6,height=12)
p<-DotPlot(honeybee, features = TF,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()

# the peak in ORN 

open_peak <- ORN_cluster_peak[ORN_cluster_peak$avg_log2FC > 0.25, ]$gene
close_peak <- ORN_cluster_peak[ORN_cluster_peak$avg_log2FC < -0.25, ]$gene
closest_open <- ClosestFeature(honeybee, unique(open_peak))
closest_close <- ClosestFeature(honeybee, close_peak)

# save the ORN specific peak by bed format 



# TF exp 
TF<-c("LOC552100",#ovo
  "LOC410499",#prd
  "LOC411079",#grh
  "Usp","Dl",
  "LOC411207",#slbo
  #"LOC10057220",#ro(rough)
  "LOC100576147", "LOC724740"  ,  "LOC725966"  ,  "LOC726165"#fkh
  )
label<-c("ovo",#ovo
  "prd",#prd
  "grh",#grh
  "Usp","Dl",
  "slbo",#slbo
  #"ro",#ro(rough)
  "fkh FD4", "slp2"  ,  "fkh crocodile"  ,  "slp1"#fkh
  )
DefaultAssay(ORN)<-"raw_RNA"
pdf("./ORN/remove_nopower/TF-ORN-exp.pdf",width=6,height=12)
p<-DotPlot(ORN, features = TF,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()
 




