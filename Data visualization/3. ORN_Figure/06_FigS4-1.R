# # version 2023.10.7
# Fig two single OR and two coexp OR pairs

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
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

# FigS4K: LOC725052 LOC100576839
pdf('./00_Figure/FigS4/FigS4K-multiple-OR-FeaturePlot.pdf', width=16, height=4)
# Visualize co-expression of two features simultaneously
FeaturePlot(ORN, features = c("LOC100576839", "LOC725052"),cols=c("lightgrey", "#E31A1C", "#4DAE49"), max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("p3_3:Or85b_LOC102655285")
dev.off()

DefaultAssay(ORN) <- "SCT"
### Plotting scatter plot: 1. single-cell level, 2. cluster level
ORN$LOC725052_UMIs <- ORN@assays$SCT@counts['LOC725052',]
ORN$LOC100576839_UMIs <- ORN@assays$SCT@counts['LOC100576839',]
#.....................................................................................
#  LOC725052 vs. LOC100576839 UMI
# ....................................................................................
#  single-cell level
library(cowplot)
pmain <- ORN@meta.data %>%
  ggplot( aes(LOC725052_UMIs, LOC100576839_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = ORN@meta.data, aes(x = LOC725052_UMIs), fill="#B31416") +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ORN@meta.data, aes(x = LOC100576839_UMIs), fill="#4DAE49") + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.3, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.3, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS4/FigS4L-LOC725052vsLOC100576839_UMI.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 1000, height =1000,p3)




