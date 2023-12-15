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
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
DefaultAssay(ORN)<-"raw_RNA"


library(scCustomize)
colors_list <- c(myUmapcolors,myUmapcolors)
Orco<- c("Or2","LOC552552","LOC551704","LOC726019")
pdf('./00_Figure/FigS2/FigS2B-Orcocoreceptor_VlnPlot_RNA.pdf',width=25, height=8)
VlnPlot(ORN,col=colors_list,add.noise=F,log=T,stack = F, features = Orco, ncol = 1, pt.size = 0)+ scale_y_log10()
#Stacked_VlnPlot(seurat_object = ORN,pt.size=0.1, features = Orco, x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = colors_list)+ scale_y_log10()
dev.off()

DefaultAssay(ORN)<-"RNA"
pdf('./00_Figure/FigS2/FigS2C-Orco-receptor-FeaturePlot.pdf', width=17, height=4)
p4<-FeaturePlot(ORN, reduction = 'tsne.rna',features = c("Or2") ,max.cutoff = 10,order=TRUE, ncol = 1)
p1<-FeaturePlot(ORN, reduction = 'tsne.rna',features = c("LOC552552"),max.cutoff =2.5, order=TRUE,ncol = 1)
p2<-FeaturePlot(ORN, reduction = 'tsne.rna',features = c("LOC726019"),max.cutoff =2.5, order=TRUE,ncol = 1)
p3<-FeaturePlot(ORN, reduction = 'tsne.rna',features = c("LOC551704"),max.cutoff =2.5, order=TRUE,ncol = 1)
p4|p1|p2|p3
dev.off()


# FigS2D  SCATTER PLOTS (Or2 Vs potential Orco)
DefaultAssay(ORN) <- "RNA"
scatterColors <- c('#A06CB4', '#DF6C78', '#911A2E', '#CD9139', '#B4B4B6', '#21918c', '#3b528b', '#440154')
### Plotting scatter plot: 1. single-cell level, 2. cluster level
ORN$Or2_UMIs <- ORN@assays$SCT@counts['Or2',]
ORN$LOC552552_UMIs <- ORN@assays$SCT@counts['LOC552552',]
ORN$LOC726019_UMIs <- ORN@assays$SCT@counts['LOC726019',]
ORN$LOC551704_UMIs <- ORN@assays$SCT@counts['LOC551704',]
ORN$Or2_Exp <- ORN@assays$RNA@data['Or2',]
ORN$LOC552552_Exp <- ORN@assays$RNA@data['LOC552552',]
ORN$LOC726019_Exp <- ORN@assays$RNA@data['LOC726019',]
ORN$LOC551704_Exp <- ORN@assays$RNA@data['LOC551704',]

library(tidyr)
library(ggrepel)
library(cowplot)
# 2. in cluster level
DefaultAssay(ORN)<-"SCT"
p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC552552')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC552552, color=id)) + 
  geom_point() +
  scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC552552')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC552552), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS2/FigS2D-1-neuronFigures_Or2-v-LOC552552_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC726019')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 
pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC726019, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC726019')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC726019), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS2/FigS2D-2-neuronFigures_Or2-v-LOC726019_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)


plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC551704')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC551704, color=id)) + 
  geom_point() +scale_color_manual(values=colors_list)+
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC551704')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC551704), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS2/FigS2D-3-neuronFigures_Or2-v-LOC551704_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# fly and honeybee OR and GR RMSD 

# FigS2E Orco and Other Ir gene RMSD heatmap 

# fly Ir 
cd /md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/Orco_pdb/
***fly and honybee coreceptor 
ls *.pdb > pdb_list.txt
for file1 in $(<pdb_list.txt)
do
 for file2 in $(<pdb_list.txt)
 do 
  python ../RMSD.py $file1 $file2 >> RMSD_result.txt
 done
done
sed -i '/^Matrix/d' RMSD_result.txt 

RMSD <- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/Orco_pdb/RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/Orco_pdb/pdb_list.txt")
ID <- gsub(".pdb","",pdb$V1)
ID<-gsub(".*_","",ID)

data<- matrix(ncol=length(ID),nrow=length(ID))
rownames(data)<- ID
colnames(data)<- ID
i=1
for (row in 1:length(ID)){
    for(col in 1:length(ID)){
        data[row,col]=RMSD[i];
        i=i+1
    }
}

# plot heatmap 
library(pheatmap)
label_pheatmap<- data.frame(species=c(rep("D.melanogaster",4),rep("Apis mellifera",4)))
rownames(label_pheatmap) <- colnames(data)
data[data==0]<- 0.5
data<- log2(data)
ann_colors<-list(
species = c("D.melanogaster"="#7985B6", "Apis mellifera"="#C7B6E1"))

pdf("./00_Figure/FigS2/FigS2E_honeybee_fly_OR_GR_RMSD_heatmap.pdf",width=7,height=5)
pheatmap(data,
      cluster_cols = T,
      cluster_rows = T,
      color = colorRampPalette(c("#1C214F","#3082BD","#F1ECEC"))(100),
      annotation_col = label_pheatmap,
      annotation_colors = ann_colors,
      annotation_row = label_pheatmap,
      annotation_legend = TRUE,
      show_rownames=T,
      show_colnames=T
 )
dev.off()





# FigS2F: Orco track plot 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN) <- "peaks_ORN_subcluster"
# first compute the GC content for each peak
ORN <- RegionStats(ORN, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(ORN)$tx_id <-Annotation(ORN)$gene_name
#features<-c("Or2","LOC411079","LOC410151","LOC406073","LOC409780","Obp5","Obp11","Obp4")
# link peaks to genes

ORN <- LinkPeaks(
  object = ORN,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "RNA",
  genes.use = Orco
)
######Visulize track and RNA exp######
idents.plot <- Idents(ORN)

pdf("./00_Figure/FigS2/FigS2F-Orco_gene-peaktrack-RNAexp-WNN.pdf",height=12,width=6)
# Or2
p1 <- CoveragePlot(
  object = ORN,
  region = "Group1-5723000-5724000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)

# LOC552552
p2 <- CoveragePlot(
  object = ORN,
  region = "Group8-7228000-7230000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)


# LOC726019
p3 <- CoveragePlot(
  object = ORN,
  region = "Group11-16143000-16145000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)
#LOC726019

p4 <- CoveragePlot(
  object = ORN,
  region = "Group14-4996500-4997500",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)
set<-c(myUmapcolors,myUmapcolors)
p1<-p1& scale_fill_manual(values=set)
p2<-p2& scale_fill_manual(values=set) & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set) & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p4<-p4& scale_fill_manual(values=set) & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p1|p2|p3|p4
dev.off()

