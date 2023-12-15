
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
set.seed(1234)

######honey marker gene annotation ######
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")
#####Annotate cells by RNA assay################
DefaultAssay(honeybee) <- "RNA"
Idents(honeybee)<-honeybee$seurat_clusters
#UMAP for all cell types 
features <- c( "LOC411079",####Epithelial cell
  "LOC413466","LOC410151",#glial cells
  "LOC411597","LOC413021",####GB46795=Ppn,Hemocyte
  "LOC412455",#GB54378=MHC,Muscle cells
  "SsRbeta","LOC406073","LOC726353","LOC724520",#sheath cells,#Auxiliary cells--thecogen (Th)
  "LOC413870",#socket cells/Tormogen
  "LOC727035",#shaft cells /Trichogen
  "Obp4" ,"Obp5" ,"Obp11","Obp12",###OBPs,Auxiliary cells--trichogen (Tr)
  "LOC410657",#听觉神经
  "Syt1","LOC408554","brp","LOC552555","LOC725310",#GB42866=Brp,YKT6/GB52698=Syt,Neuron
  "LOC413063","LOC410689","LOC410353","LOC409780",
  "LOC408777",###=nompC,Johnston organ neuron
  "LOC412949","LOC727431",#IR,ionotropic receptor 
  "Or2"
  )

label<-c("GRH",#Epithelial cell
  "OAZ","repo",#glial cells
  "hemocytin","Ppn",####Hemocyte
  "MHC",#Muscle cells
  "SSRBETA","PROS","nompA","nompA-1",#sheath cells,#Auxiliary cells--thecogen (Th)
  "SU(H)",#socket cells
  "POXN",#shaft cells
  "OBP4" ,"OBP5" ,"OBP11","OBP12",###OBPs,Auxiliary cells--trichogen (Tr)
  "ACJ6",#听觉神经
  "Syt1","CADN","BRP","YKT6","SYT",#Neuron
  "pebble","Elav-like2","Elav-likem4","Elav-likem2",#neuron
  "nompC",#Johnston organ neuron
  "IR25a","GR64f",
  "Or2")
pdf("./02_All_celltype/honeybee_cluster-annotation-all_celltype.pdf",width=12,height=8)
p<-DotPlot(honeybee, features = features,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()

Idents(honeybee)<-honeybee$seurat_clusters
#####further annotation########
honeybee <- RenameIdents(
  object = honeybee,
  '0' = 'Neuron',
  '1' = 'Neuron',
  '2' = 'Neuron',
  '3' = 'Neuron',
  '4' = 'Neuron',
  '5' = 'Neuron',
  '6' = 'Neuron',
  '7' = 'Neuron',
  '8' = 'Neuron',
  '9' = 'Neuron',
  '10' = 'Sheath cell',
  '11' = 'Obp5+support cell',
  '12' = 'Neuron',
  '13' = 'Epithelial cell',
  '14' = 'Neuron',
  '15' = 'Unannotated',
  '16' = 'Glial cell',
  '17' = 'Glial cell',
  '18' = 'Neuron',
  '19' = 'Obp4+support cell',
  '20' = 'Glial cell',
  '21' = 'Neuron',
  '22' = 'Neuron',
  '23' = 'Neuron',
  '24' = 'Neuron',
  '25' = 'Neuron',
  '26' = 'Neuron',
  '27' = 'Obp5+support cell'
  )
honeybee@meta.data$Annotation<-Idents(honeybee)
table(honeybee$Annotation,honeybee$orig.ident)

                      NE Nurse Forager
  Neuron            2435  2050    2847
  Sheath cell         63   101     151
  Obp5+support cell   86    86      96
  Epithelial cell     25    63     120
  Unannotated         46    63      66
  Glial cell          90   141     206
  Obp4+support cell   28    51      54

honeybee$Annotation<-factor(honeybee$Annotation,levels=c("Neuron",'Epithelial cell','Glial cell',
  'Sheath cell','Obp4+support cell','Obp5+support cell',"Unannotated"))
Idents(honeybee)<-honeybee$Annotation;

pdf("./02_All_celltype/honeybee_annotation_allcelltype_UMAP.pdf",width=6,height=5)
DimPlot(honeybee, label = T, repel = TRUE, cols=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#B3B3B3"), reduction = "wnn.umap",group.by = "Annotation")
DimPlot(honeybee, label = F, repel = TRUE, cols=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#B3B3B3"), reduction = "wnn.umap",group.by = "Annotation")+ ggtitle("")
dev.off()

# Save rds have annotation information 
DefaultAssay(honeybee) <- "RNA"
saveRDS(honeybee,"./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")

# Marker gene for different cell types (color on UMAP ) 

pdf('./02_All_celltype/All_Marker_gene_FeaturePlot_WNN.pdf', width=5.5, height=5)
for (i in 1:28){
  p<-FeaturePlot(honeybee,order=T, reduction = 'wnn.umap',max.cutoff = 10, features = features[i], ncol = 1)
  print(p)
}
dev.off()

##Track for Marker genes promoters
Idents(honeybee)<-honeybee$Annotation
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(honeybee) <- "peaks"
# first compute the GC content for each peak
honeybee <- RegionStats(honeybee, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(honeybee)$tx_id <- Annotation(honeybee)$gene_name 
#features<-c("Or2","LOC411079","LOC410151","LOC406073","LOC409780","Obp5","Obp11","Obp4")

# link peaks to genes
honeybee <- LinkPeaks(
  object = honeybee,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = features
)
######Visulize track and RNA exp######
idents.plot <- Idents(honeybee)

pdf("./02_All_celltype/Marker_gene-peaktrack-RNAexp-WNN.pdf",height=8,width=8)
for(i in features){
  print(i)
  p1 <- CoveragePlot(
  object = honeybee,
  region = i,
  features = i,
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 500,
  annotation=TRUE,
  extend.downstream = 500
)
print(p1)}

dev.off()

#Epithelial cell 
p1 <- CoveragePlot(
  object = honeybee,
  region = "Group5-6519000-6520000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)
#glial cell 
p2 <- CoveragePlot(
  object = honeybee,
  region = "Group15-2481500-2482500",
  #extend.upstream = 20,
  annotation=TRUE,peaks = F,
  #extend.downstream = -1000,
  links=F
)
#Sheath 
p3 <- CoveragePlot(
  object = honeybee,
  region = "Group5-3433000-3433900",
  #extend.upstream = 500,
  annotation=TRUE,peaks = F,
  #extend.downstream = -1000,
  links=F
)

#Neuron 
p5 <- CoveragePlot(
  object = honeybee,
  region = "Group1-3299000-3300000",
  #extend.upstream = 0,
  annotation=TRUE,peaks = F,
  #extend.downstream = -10000,
  links=F
)
#ORN
p6 <- CoveragePlot(
  object = honeybee,
  region = "Group1-5722000-5725000",
  #extend.upstream = 500,
  annotation=TRUE,peaks = F,
  #extend.downstream = -23000,
  links=F
)
#Obp5
p7 <- CoveragePlot(
  object = honeybee,
  region = "Group9-11944500-11945500",
  #extend.upstream = 0,
  annotation=TRUE,peaks = F,
  #extend.downstream = -10000,
  links=F
)


#Obp12
p9 <- CoveragePlot(
  object = honeybee,
  region = "Group9-11981000-11982000",
  #extend.upstream = 0,
  annotation=TRUE,peaks = F,
  #extend.downstream = -10000,
  links=F
)

set<-c('#F1BB72', '#764AF1', "#B03B2A","#E95C58",'#749F82','#548CA8','#E59CC4', '#AB3282', '#23452F','#E5D2DD' )

p1<-p1& scale_fill_manual(values=set)&labs(title="GRH")
p2<-p2& scale_fill_manual(values=set)&labs(title="repo") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)&labs(title="PROS") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
#p4<-p4& scale_fill_manual(values=set)&labs(title="Ppn") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p5<-p5& scale_fill_manual(values=set)&labs(title="Elavl2") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p6<-p6& scale_fill_manual(values=set)&labs(title="Orco") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p7<-p7& scale_fill_manual(values=set)&labs(title="OBP5") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p8<-p8& scale_fill_manual(values=set)&labs(title="OBP11") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p9<-p9& scale_fill_manual(values=set)&labs(title="OBP4") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
pdf("./02_All_celltype/Marker_gene-select-peaktrack-WNN.pdf",height=8,width=16) 
p1|p2|p3|p5|p6|p7|p8|p9
dev.off()

###
pdf("./02_All_celltype/Marker_select_violin.pdf")
p<-ExpressionPlot(
  object = honeybee,
  features = c("LOC411079","LOC410151","LOC406073","LOC409780","Or2","Obp5","Obp11","Obp4"),
  assay = "RNA"
)
p& scale_fill_manual(values=set)
dev.off()
