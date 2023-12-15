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

ORN <- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_latest.rds")

DefaultAssay(ORN)<-"raw_RNA"
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C"  , "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
# Fig2A:Unsupervised clustering ORN 
onecluster <- readRDS("./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")
pdf("./00_Figure/Fig2/Fig2A-Unsupervised_ORN_cluster_WNN.pdf",width=6,height=6)
DimPlot(onecluster, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = F, label.size = 5, repel = TRUE) & NoLegend() 
dev.off()

# Fig2B:supervised clustering ORN 
pdf("./00_Figure/Fig2/Fig2B-Supervised_ORN_cluster_WNN_remove_nopower.pdf",width=6,height=6)
DimPlot(ORN, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = F, label.size = 5, repel = TRUE)& NoLegend() 
dev.off();

# Fig2C:the proportion of detected OR
ORN_all<- readRDS("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/Unsupervised_ORN_cluster_WNN_add_subcluster.rds")
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
all_OR_gene<- setdiff(all_OR_gene,Orco)
detected_OR_gene<- rownames(ORN_all)[which(rownames(ORN_all)%in%all_OR_gene)]
power_OR_gene<- dotplot_feature
power_OR_gene<- setdiff(power_OR_gene,Orco)
data<- data.frame(gene=c("total_OR","detected_OR","powerful_OR","powerful_cluster"),
    number=c(length(all_OR_gene),length(detected_OR_gene),length(power_OR_gene),length(levels(ORN))))
data$gene<- factor(data$gene,levels=c("total_OR","detected_OR","powerful_OR","powerful_cluster"))
pdf("./00_Figure/Fig2/Fig2C-number_of_detected_OR_or_cluster.pdf",width=3,height=4)
p<-ggplot(data = data, aes_string(x = "gene", y = "number")) +  
        xlab(" ") +
        ylab("Number of OR or cluster") + 
       # scale_fill_manual(values =  "#89C75F" ) + 
        geom_bar( stat = "identity",width=0.6,color = 'white', fill="#6580A8") +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
round(data$number, 2)
p+geom_text(aes(label = number), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

# Fig2D: tree + downsample heatmap + barplot
colors_for_exp_pattern<- c("#476D87","#E95C59")
DefaultAssay(ORN) <- "integratedRNA_onecluster"
object<- ORN
embeddings <- Embeddings(object = object, reduction = "obj_features_pca")[,1:50]
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
ORN<- object
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])

library(ggtree)
pdf("./00_Figure/Fig2/Fig2D-1-remove_nopower-cluster-ORN-tree-cosine.pdf",width=5,height=8)
tree <- groupOTU(data.tree, .node=as.character(multiOR_cluster))
ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
dev.off()

library(ggtree)
pdf("./00_Figure/Fig4/Fig4A-remove_nopower-cluster-ORN-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
dev.off()

# Fig2D-2: chemoreceptor gene heatmap
# downsample to same number
target_cell_count <- 32 
seurat_object<- ORN
downsampled_clusters <- list()
for (cluster_id in unique(seurat_object$cell_group)) {
  cluster_cells <- which(seurat_object$cell_group == cluster_id)
  if (length(cluster_cells) > target_cell_count) {
    downsampled_cells <- sample(cluster_cells, target_cell_count)
  } else {
    downsampled_cells <- cluster_cells
  }
  downsampled_clusters[[cluster_id]] <- downsampled_cells
}

downsampled_ORN <- subset(ORN, cells = unlist(downsampled_clusters))
saveRDS(downsampled_ORN,"./00_Figure/Fig2/Fig2D_downsampled_ORN.rds")

# change the OR gene name 
downsampled_ORN<- readRDS("./00_Figure/Fig2/Fig2D_downsampled_ORN.rds")
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")
DefaultAssay(downsampled_ORN)<- "SCT"
all_gene<- rownames(downsampled_ORN)
RenameGenesSeurat_SCT <- function(obj = ls.Seurat[[i]], newnames = tmp) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts and @data ")
  SCT <- obj@assays$SCT
  if (nrow(SCT) == length(newnames)) {
    if (length(SCT@counts)) SCT@counts@Dimnames[[1]]            <- newnames
    if (length(SCT@data)) SCT@data@Dimnames[[1]]                <- newnames
    #if (length(SCT@scale.data)) SCT@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(SCT) != nrow(newnames)"}
  obj@assays$SCT <- SCT
  return(obj)
}

for (i in 1:length(all_gene)){
    if(all_gene[i] %in% ORgene_name_trans$OR_gene){
        tmp=ORgene_name_trans[which(ORgene_name_trans$OR_gene==all_gene[i]),]$last_name;
        all_gene[i]=tmp
    }
}

ORN_trans <- RenameGenesSeurat_SCT(downsampled_ORN,newnames = all_gene)

dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)))

for (i in 1:length(dotplot_feature)){
    if(dotplot_feature[i] %in% ORgene_name_trans$OR_gene){
        tmp=ORgene_name_trans[which(ORgene_name_trans$OR_gene==dotplot_feature[i]),]$last_name;
        dotplot_feature[i]=tmp
    }
}

DefaultAssay(ORN_trans)<-"SCT"
ORN_trans<- ScaleData(ORN_trans,features=dotplot_feature)
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

pdf("./00_Figure/Fig2/Fig2D-2-remove_nopower-chemoreceptor_gene_heatmap-orderbytree_yasuo.pdf",width=8, height=13)
#DoHeatmap(ORN_trans,disp.max = 1,slot = "data",size = 2.5,angle = 315, features = dotplot_feature,draw.lines = FALSE,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = solarExtra[4:9])
#DoHeatmap(ORN_trans,disp.max = 1,slot = "data",size = 2,angle = 315, features = dotplot_feature,draw.lines = FALSE,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = solarExtra[4:8])
#DoHeatmap(ORN_trans,disp.max = 2,disp.min = -2,slot = "scale.data",size = 2,angle = 315, features = dotplot_feature,draw.lines = FALSE,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = solarExtra[4:8])
#DoHeatmap(ORN_trans,slot = "scale.data",size = 2,angle = 315, features = dotplot_feature,draw.lines = FALSE,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = solarExtra[2:8])
#DoHeatmap(ORN_trans,slot = "scale.data",size = 2,angle = 315, features = dotplot_feature,draw.lines = FALSE,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(ORN_trans,slot = "scale.data",size = 2,angle = 315, features = dotplot_feature,draw.lines = FALSE,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = solarExtra[4:8])
dev.off()

# Fig2D-3 # of cells 
ORN$cell_group<- factor(ORN$cell_group,levels=as.character(1:60))
cluster_cellnumber<-as.data.frame(table(Idents(ORN)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(ORN))
cluster_cellnumber$color<-c(myUmapcolors,myUmapcolors)[1:length(levels(ORN))]
#cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)
pdf("./00_Figure/Fig2/Fig2D-3-remove_nopower_ORN_cluster_cellnumber.pdf",width=6,height=8)
p<-ggplot(data = cluster_cellnumber, aes_string(x = "cluster", y = "number", 
        fill = "cluster")) +  xlab(" ") + ylab("# of cells") + 
        scale_fill_manual(values = cluster_cellnumber$color) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5))+coord_flip();
p
#add gene number in plot 
p+geom_text(aes(label = number), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

# Fig 2E 
cluster_number<-as.data.frame(table(table(dotplot_data$id)))
cluster_number$Var1<-as.character(cluster_number$Var1)
cluster_number[4,1]<-">3"
cluster_number[4,2]<-5
cluster_number<-cluster_number[1:4,]
cluster_number$percent<- cluster_number$Freq/sum(cluster_number$Freq)*100;
cluster_number$Var1<- factor(cluster_number$Var1,levels=c("1","2","3",">3"))

pdf("./00_Figure/Fig2/Fig2E-multi_OR_percent.pdf",width=3,height=4)
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "percent")) +  
        xlab("# chemosensory receptor expressed") +
        ylab("% percent of cluster") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity",width=0.6,color = 'black', fill='grey') +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
round(cluster_number$percent, 2)
p+geom_text(aes(label = Freq), size = 3, hjust = 0.5, vjust = 3) 
dev.off();
# Fig2F: cross species expression pattern statics 
# single OR  vs multiple OR cluster porportion 
single_OR_cluster<- length(levels(ORN))-length(multiOR_cluster)
multiple_OR_cluster<- length(multiOR_cluster)
Apis_mellifera<-c(single_OR_cluster,multiple_OR_cluster);

#fly need to use the public datasets
Drosophila_melanogaster<-c(40,5)

# total:42
Aedes_aegypti<-c(23,19);
cross_species_cluster_number<-data.frame(species=c(rep("Apis mellifera",2),rep("D.melanogaster",2),rep("Ae.aegypti",2)),
  exp_pattern=rep(c("single_OR","multiple_OR"),3),
  cluster_number=c(Apis_mellifera,Drosophila_melanogaster,Aedes_aegypti))
cross_species_cluster_number$exp_pattern<-factor(cross_species_cluster_number$exp_pattern,levels=c("single_OR","multiple_OR"));
cross_species_cluster_number$species<-factor(cross_species_cluster_number$species,levels=c("Apis mellifera","D.melanogaster","Ae.aegypti"))

pdf("./00_Figure/Fig2/Fig2F-cross_species_expression_pattern_statics.pdf",width=4,height=4)
ggplot(data = cross_species_cluster_number, aes_string(x = "species", y = "cluster_number", 
        fill = "exp_pattern")) +  xlab(" ") + ylab("% Percent of cells") + 
        scale_fill_manual(values = colors_for_exp_pattern) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));

p<-ggplot(data = cross_species_cluster_number, aes_string(x = "species", y = "cluster_number", 
        fill = "exp_pattern")) +  xlab(" ") + ylab("# of cluster") + 
        scale_fill_manual(values = colors_for_exp_pattern) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = cluster_number), size = 3, hjust = 0.5, vjust = 3, position = "stack") 
dev.off();

# Fig2G:
OR_number<-as.data.frame(table(table(dotplot_data$features.plot)));
#OR_number<-OR_number[-1,]
OR_number$percent<- OR_number$Freq/sum(OR_number$Freq)*100;
pdf("./00_Figure/Fig2/Fig2G-OR_exp_in_multicluster_percent.pdf",width=3,height=4)
p<-ggplot(data = OR_number, aes_string(x = "Var1", y = "percent")) +  
        xlab("# of cluster OR expressed") +
        ylab("% percent of chemosensory receptor") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity",width=0.6,color = 'black', fill='grey') +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
p+geom_text(aes(label = Freq), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

# Fig2H: OR pep tree + transcript dist heatmap:
# OR correlation  heatmap and sequence tree 
# sequence tree 
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
#OR2 is placed in the last column;
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

# corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature,]
# color by the group info 
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)
data.tree <- tree
pdf("./00_Figure/Fig2/Fig2H-a-OR_sequence_protein_similarity-tree_add_groupinfo_heatmap.pdf",width=25,height=16)
ggtree(data.tree,ladderize = FALSE, branch.length = "none") + geom_tiplab(size=3) + 
theme(legend.position = "right")+  geom_treescale()
dev.off()
m<-ggtree(data.tree,ladderize = FALSE, branch.length = "none")+ geom_tiplab(size=3) + theme(legend.position = "right")+ scale_color_manual (values =myUmapcolors[5:30]) 
gene_order<-na.omit(m$data[order(m$data$y),]$label)
gene_order<-as.character(gene_order);

# OR color bar in heatmap 

library(pheatmap)
DefaultAssay(ORN)<-"SCT"
matrix<- ORN@assays$SCT[dotplot_feature,]
cor_data<-cor(as.data.frame(t(matrix)))
cosine_dist <- (1-cosine(cor_data))
dist_data <- cor_data[rev(gene_order),rev(gene_order)]
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
label_pheatmap<- data.frame(Group=Group_info$seqnames)
rownames(label_pheatmap) <- Group_info$gene_name
ann_colors<-list(Group =    c("Group1"=myUmapcolors[1],    "Group2"=myUmapcolors[2],    "Group3"=myUmapcolors[3],    "Group4"=myUmapcolors[4],    "Group5"=myUmapcolors[5],    "Group6"=myUmapcolors[6],    "Group7"=myUmapcolors[7],    "Group8"=myUmapcolors[8],    "Group9"=myUmapcolors[9],    "Group10"=myUmapcolors[10],
    "Group11"=myUmapcolors[11],    "Group12"=myUmapcolors[12],    "Group13"=myUmapcolors[13],    "Group14"=myUmapcolors[14],    "Group15"=myUmapcolors[15],    "Group16"=myUmapcolors[16],    "GroupUN3"=myUmapcolors[17],
    "GroupUN226"=myUmapcolors[18],    "GroupUN243"=myUmapcolors[19],    "GroupUN248"=myUmapcolors[20]))
pdf("./00_Figure/Fig2/Fig2H-b-OR_in_dotplot_sequence_correlation.pdf",width=17,height=16)
pheatmap(dist_data,
        border = F,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("white", "#3082BD","#1C214F"))(100),
         annotation_legend = TRUE,
         annotation_colors = ann_colors,
         annotation_row = label_pheatmap,
         show_rownames=T,
         show_colnames=T
    )
pheatmap(dist_data,
        border = F,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("white", "#3082BD","#1C214F"))(100),
         annotation_legend = TRUE,
         annotation_colors = ann_colors,
         annotation_row = label_pheatmap,
         show_rownames=T,
         show_colnames=T
    )
dev.off()

# Fig2I and Fig2G: genomic distance 

#the distance of OR pair in the same cluster and random OR pair 
DefaultAssay(ORN)<-"ATAC"
gene_transcript<-Annotation(ORN)[which(Annotation(ORN)$type=="transcript"),];
OR_gene_transcript<-gene_transcript[gene_transcript$gene_name%in%dotplot_feature,]
OR_gene_transcript<-OR_gene_transcript[!duplicated(OR_gene_transcript$gene_name),]
AmelchrNameLength <- read.table("/md01/nieyg/ref/10X/honeybee/de_novo_antenna/Amel_antenan/star/chrNameLength.txt", sep="\t", stringsAsFactors=F) 
AmelchrNameLength<-AmelchrNameLength[which(AmelchrNameLength$V1%in%as.character(unique(OR_gene_transcript@seqnames@values))),]
data<-as.data.frame(sort(OR_gene_transcript));
data<-data[,c(1:3,12,4,5)]
write.table(data,"OR_gene_transcript.bed",sep="\t",row.names=F,col.names=F)

# #Shell
# sed -i 's/"//g' OR_gene_transcript.bed;
# bedtools sort -i OR_gene_transcript.bed > OR_gene_transcript_sorted.bed
# bedtools closest -s -d -io -N -a OR_gene_transcript_sorted.bed -b OR_gene_transcript_sorted.bed > output.bed
# awk '{print $NF,"\t",$1,"\t",$4,"\t",$10}' output.bed > closestOlfrGenes.txt
# sort -n closestOlfrGenes.txt | awk '$1 > 0 {print $0}' > sortedClosestOlfrGenes.txt

# back R 
# make OR pair 
OR_pair<-data.frame()
remaining_gene<-dotplot_feature
for(gene1 in dotplot_feature){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    data_subset<-data.frame(gene1=gene1,gene2=gene2);
    OR_pair<-rbind(OR_pair,data_subset)
  }
}
OR_pair$gene1_seqname <- data[match(OR_pair$gene1,data$gene_name),1]
OR_pair$gene2_seqname <- data[match(OR_pair$gene2,data$gene_name),1]

#Same seqname?
for (i in 1:nrow(OR_pair)){
    if(OR_pair[i,]$gene1_seqname==OR_pair[i,]$gene2_seqname){
        OR_pair$same_seqname[i]="Yes";}
        else{OR_pair$same_seqname[i]="No"}
};

#Same cluster?
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$gene1;
    gene2<-OR_pair[i,]$gene2;
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id;
    if(length(which(duplicated(cluster_id)))){
        OR_pair$same_cluster[i]="co-exp";
    }else{OR_pair$same_cluster[i]="non co-exp";}
    
};

# The probability of OR pair appearing on the same chromosome
#> table(OR_pair$same_seqname,OR_pair$same_cluster)
      co-exp non co-exp
  No       0       2341
  Yes     60        759


probability<-data.frame(OR_pair_type=c("non co-exp","co-exp"),
    probability=c(759/2341,60/60),count=c(60,759));
probability$OR_pair_type<-factor(probability$OR_pair_type,levels=c("non co-exp","co-exp"))
pdf("./00_Figure/Fig2/Fig2I-2G-The_probability_OR_pair_appearing_on_same_chromosome.pdf",width=6,height=3)
p<-ggplot(data = probability, aes_string(x = "OR_pair_type", y = "probability",color="OR_pair_type")) +  
        xlab("OR pair type") +
        ylab("The probability on the same chromosome") + 
        geom_bar(stat = "identity", width = 0.6,fill="white") +
        theme_classic()+guides(color='none')+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)) + 
        scale_color_manual(values = colors_for_exp_pattern);
#add number in plot 
a<- p+geom_text(aes(label = count), size = 3, hjust = 0.5, vjust = 3)
# the distance 
same_chrom<-OR_pair[OR_pair$same_seqname=="Yes",]
for (i in 1:nrow(same_chrom)){
    gene1<-same_chrom[i,]$gene1;
    gene2<-same_chrom[i,]$gene2;
    OR_pair_transcript<-OR_gene_transcript[OR_gene_transcript$gene_name%in%c(gene1,gene2),]
    OR_pair_transcript<-sort(OR_pair_transcript);
    same_chrom$dist[i]<-width(gaps(OR_pair_transcript))[2];
}

same_chrom$same_cluster<-factor(same_chrom$same_cluster,levels=c("non co-exp","co-exp"))
b<- ggboxplot(same_chrom, x="same_cluster", y="dist",color = "same_cluster",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("Genomic distance")+
scale_color_manual(values = colors_for_exp_pattern)
a|b 
dev.off()


# Fig2K: OR transcriptomic distance 
# the transcript distance between coexp and uncoexp
library(pheatmap)
DefaultAssay(ORN)<-"SCT"
matrix<- ORN@assays$SCT[dotplot_feature,]
cor_data<-cor(as.data.frame(t(matrix)))
cosine_dist <- (1-cosine(cor_data));

same_cluster <- c()
not_same_cluster <- c()
remaining_gene<-rownames(cosine_dist)
data<- cosine_dist
for(gene1 in remaining_gene){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id
    dist<-data[gene1,gene2]
    if(length(which(duplicated(cluster_id)))){
        same_cluster<-c(same_cluster,dist)
    }else{not_same_cluster<-c(not_same_cluster,dist)}
  }
}
type<-c(rep("co-exp",length(same_cluster)),rep("non co-exp",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data2<-data.frame(type,var)
data2$type<-factor(data2$type,levels=c("non co-exp","co-exp"))
library(ggpubr)
pdf("./00_Figure/Fig2/Fig2K-OR_transcript_dist_distribution.pdf",width=3,height=3)
ggboxplot(data2, x="type", y="var", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("transcript distance")+
scale_color_manual(values = colors_for_exp_pattern)
dev.off()

# Fig2L: sequence similarity
colors_for_exp_pattern<- c("#476D87","#E95C59")
#Fig 3D: Sequence similarity
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature,]
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
tree<-as.phylo(clust)
# the similarity of OR pair in the same cluster and random OR pair 
# All OR pair similarity matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
same_cluster <- c()
not_same_cluster <- c()
remaining_gene<-dotplot_feature
for(gene1 in dotplot_feature){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id
    dist<-dist_percent[gene1,gene2]
    if(length(which(duplicated(cluster_id)))){
        same_cluster<-c(same_cluster,dist)
    }else{not_same_cluster<-c(not_same_cluster,dist)}
  }
}
type<-c(rep("co-exp",length(same_cluster)),rep("non co-exp",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("non co-exp","co-exp"))
library(ggpubr)
library(cowplot)

pdf("./00_Figure/Fig2/Fig2L-OR_sequence_similarity_distribution.pdf",width=3,height=3)
p1<-ggplot(data, aes(x=var, fill=type)) + xlab("% OR sequence similarity")+
          geom_density(alpha=.25) + theme_classic()+theme(legend.position="top")+
scale_color_manual(values =colors_for_exp_pattern)
p2<-ggboxplot(data, x="type", y="var", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("% OR sequence similarity")+
scale_color_manual(values =colors_for_exp_pattern)
p2
dev.off()

# Fig2M: OR RMSD
# the RMSD (structure similarity)
RMSD <- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/rename_pdb/RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/rename_pdb/pdb_list.txt")
ID <- gsub(".pdb","",pdb$V1)
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
#data[data==0]<- 0.5
#data<- log2(data)
same_cluster <- c()
not_same_cluster <- c()
remaining_gene<-rownames(data)
for(gene1 in remaining_gene){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id
    dist<-data[gene1,gene2]
    if(length(which(duplicated(cluster_id)))){
        same_cluster<-c(same_cluster,dist)
    }else{not_same_cluster<-c(not_same_cluster,dist)}
  }
}
type<-c(rep("co-exp",length(same_cluster)),rep("non co-exp",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data2<-data.frame(type,var)
data2$type<-factor(data2$type,levels=c("non co-exp","co-exp"))
pdf("./00_Figure/Fig2/Fig2M-OR_structure_similarity_distribution.pdf",width=3,height=3)
ggboxplot(data2, x="type", y="var", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("OR RMSD")+
scale_color_manual(values = colors_for_exp_pattern)
dev.off()












