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
After_ISAC_all <- readRDS("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/Unsupervised_ORN_cluster_WNN_add_subcluster.rds")
Before_ISAC <- readRDS("./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")

# FigS2A: Before ISAC
# FigS2A-1: cluster trans dist tree
DefaultAssay(Before_ISAC) <- "integratedRNA_onecluster"
obj <- Before_ISAC
object<- obj
Idents(obj)<- obj$seurat_clusters
DefaultAssay(obj)<-"integratedRNA_onecluster"
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
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

library(ggtree)
pdf("./00_Figure/FigS2/FigS2A-1-Before_ISAC-cluster-ORN-tree-cosine.pdf",width=5,height=14)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()

# FigS2A-2: chemoreceptor gene heatmap
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)
Idents(Before_ISAC)<-factor(Before_ISAC$seurat_clusters,levels=cluster_order)

# change the OR gene name 
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")

DefaultAssay(Before_ISAC)<- "SCT"
all_gene<- rownames(Before_ISAC)

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

Before_ISAC_trans <- RenameGenesSeurat_SCT(Before_ISAC,newnames = all_gene)

DefaultAssay(Before_ISAC)<-"SCT"
p<-DotPlot(Before_ISAC,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled>= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% all_receptor_gene])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot)),gene_order))

gene_order<- dotplot_feature

dotplot_feature<- gene_order
for (i in 1:length(dotplot_feature)){
    if(dotplot_feature[i] %in% ORgene_name_trans$OR_gene){
        tmp=ORgene_name_trans[which(ORgene_name_trans$OR_gene==dotplot_feature[i]),]$last_name;
        dotplot_feature[i]=tmp
    }
}

DefaultAssay(Before_ISAC_trans)<-"SCT"
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
dotplot_feature<- dotplot_feature[which(dotplot_feature %in% rownames(Before_ISAC_trans))]
Before_ISAC_trans<- ScaleData(Before_ISAC_trans,features=dotplot_feature)
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

pdf("./00_Figure/FigS2/FigS2A-2-Before_ISAC-chemoreceptor_gene_heatmap-orderbytree.pdf",width=15, height=20)
#DoHeatmap(Before_ISAC_trans,disp.max = 2,disp.min = -2,size = 4, slot = "scale.data",angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors,myUmapcolors,myUmapcolors), label=FALSE)+ scale_fill_gradientn(colors = solarExtra[4:8])
#DoHeatmap(Before_ISAC_trans,disp.max = 2,disp.min = -2,size = 4, slot = "scale.data",angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors,myUmapcolors,myUmapcolors), label=FALSE)+ scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(Before_ISAC_trans,draw.lines = FALSE,disp.max = 2,disp.min = -2,size = 4, slot = "scale.data",angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors,myUmapcolors,myUmapcolors), label=FALSE)+ scale_fill_gradientn(colors = solarExtra[2:8])
#DoHeatmap(Before_ISAC_trans,disp.max = 1,disp.min = -1,slot = "scale.data",size = 4,angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors,myUmapcolors,myUmapcolors), label=FALSE)+ scale_fill_gradientn(colors = solarExtra[4:8])
dev.off()

pdf("./00_Figure/FigS2/FigS2A-b-remove_nopower-dotplot-orderbytree_lastest.pdf",width=22, height=10)
p<-DotPlot(Before_ISAC,features = gene_order) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p&scale_x_discrete(labels=dotplot_feature)
p<-DotPlot(Before_ISAC,features = gene_order[-1:-4]) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p&scale_x_discrete(labels=dotplot_feature[-1:-4])
dev.off()

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]

dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))



# FigS2B: After ISAC
# FigS2B-1: cluster trans dist tree
DefaultAssay(After_ISAC_all) <- "integratedRNA_onecluster"
obj <- After_ISAC_all
obj$subcluster<- factor(obj$subcluster)
Idents(obj)<- obj$subcluster
object<- obj
DefaultAssay(obj)<-"integratedRNA_onecluster"
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
data.dims <- lapply(X = levels(x = obj), FUN = function(x) {
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

library(ggtree)
pdf("./00_Figure/FigS2/FigS2C-1-After_ISAC_all-cluster-ORN-tree-cosine.pdf",width=12,height=20)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()

# FigS2B-2: chemoreceptor gene heatmap
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)
Idents(After_ISAC_all)<-factor(After_ISAC_all$subcluster,levels=cluster_order)

# change the OR gene name 
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")

DefaultAssay(After_ISAC_all)<- "SCT"
all_gene<- rownames(After_ISAC_all)

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

After_ISAC_all_trans <- RenameGenesSeurat_SCT(After_ISAC_all,newnames = all_gene)
DefaultAssay(After_ISAC_all)<-"SCT"
p<-DotPlot(After_ISAC_all,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled>= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% all_receptor_gene])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot)),gene_order))

gene_order<- dotplot_feature

dotplot_feature<- gene_order
for (i in 1:length(dotplot_feature)){
    if(dotplot_feature[i] %in% ORgene_name_trans$OR_gene){
        tmp=ORgene_name_trans[which(ORgene_name_trans$OR_gene==dotplot_feature[i]),]$last_name;
        dotplot_feature[i]=tmp
    }
}


DefaultAssay(After_ISAC_all_trans)<-"SCT"
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
dotplot_feature<- dotplot_feature[which(dotplot_feature %in% rownames(After_ISAC_all_trans))]
After_ISAC_all_trans<- ScaleData(After_ISAC_all_trans,features=dotplot_feature)

pdf("./00_Figure/FigS2/FigS2B-2-After_ISAC_all-chemoreceptor_gene_heatmap-orderbytree.pdf",width=18, height=20)
DoHeatmap(After_ISAC_all_trans,draw.lines = FALSE,disp.max = 2,disp.min = -2,size = 4, slot = "scale.data",angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors,myUmapcolors,myUmapcolors), label=FALSE)+ scale_fill_gradientn(colors = solarExtra[2:8])
dev.off()

pdf("./00_Figure/FigS2/FigS2C-b-remove_nopower-dotplot-orderbytree_lastest.pdf",width=22, height=22)
p<-DotPlot(After_ISAC_all,features = gene_order) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p&scale_x_discrete(labels=dotplot_feature)
p<-DotPlot(After_ISAC_all,features = gene_order[-1:-4]) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p&scale_x_discrete(labels=dotplot_feature[-1:-4])

dev.off()

















