# Fig3: Cis-Elements Orchestrating Olfactory Receptor Co-Expression 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)

set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C"  , "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# Fig3A:two promoter 

# Umap plot:
library(scCustomize)
library(Nebulosa)

pdf('./00_Figure/Fig3/Fig3A-Or154_163-FeaturePlot.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("LOC102655285"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("36:Or154")
p2<- FeaturePlot(ORN, features = c("LOC107965761"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("36:Or163")
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("LOC102655285","LOC107965761"),
                                custom_palette = BlueAndRed())
p1|p2|p12
dev.off()

obj<- subset(ORN,idents=c("36"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC102655285","LOC107965761")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
colnames(ORN_matrix)<- c("Or154","Or163")
ORN_matrix<- as.data.frame(ORN_matrix)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or163,ORN_matrix$Or154,decreasing=T),]
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/Fig3/Fig3B-coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or163,ORN_matrix$Or154,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or154,ORN_matrix$Or163,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
pheatmap(t(ORN_matrix),
     cluster_cols = T,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# distinguish the coexp and the one OR exp cells
Or154_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or154>0),])
Or163_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or163>0),])
coexp_barcode<- intersect(Or154_barcode,Or163_barcode)
Or154_barcode<- setdiff(Or154_barcode,coexp_barcode)
Or163_barcode<- setdiff(Or163_barcode,coexp_barcode)

# plot the heatmap 
coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or154,coexp_matrix$Or163,decreasing=TRUE),]
Or154_matrix<- ORN_matrix[Or154_barcode,]
Or154_matrix<-Or154_matrix[order(Or154_matrix$Or154,Or154_matrix$Or163,decreasing=TRUE),]
Or163_matrix<- ORN_matrix[Or163_barcode,]
Or163_matrix<-Or163_matrix[order(Or163_matrix$Or154,Or163_matrix$Or163,decreasing=TRUE),]
ORN_matrix<- rbind(coexp_matrix,Or154_matrix,Or163_matrix)
label_pheatmap<- data.frame(label=c(rep("coexp",length(coexp_barcode)),
     rep("only_Or154",length(Or154_barcode)),
     rep("only_Or163",length(Or163_barcode))))
col <- c("#8A9FD1" ,"#C06CAB","#FEE500")

names(col)<-c("only_Or154","only_Or163","coexp")
ann_colors= list(label = col)
rownames(label_pheatmap)<- rownames(ORN_matrix)
pdf("./00_Figure/Fig3/Fig3B-coexp-heatmap-SCT-split_coexp.pdf",height=3,width=12)
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     annotation_col = label_pheatmap,
     annotation_colors = ann_colors,
     show_colnames=F
  )
dev.off()

coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or154,coexp_matrix$Or163,decreasing=TRUE),]
bk <- c(seq(0,3,by=0.1))
pdf("./00_Figure/Fig3/Fig3B-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
pheatmap(t(coexp_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colorRampPalette(colors = c("white","#CC0000"))(length(bk)),
     annotation_legend = TRUE,
     show_rownames=T,breaks=bk,
     #annotation_col = label_pheatmap,
     #annotation_colors = ann_colors,
     show_colnames=F
  )
dev.off()
     
write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or154_Or163/raw_barcode_coexp.txt",row.name=F,col.names=F)


# save the coexp barcode for track 
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or154_Or163/barcode.txt",row.name=F,col.names=F)
sed -i 's/"//g' barcode.txt
# subset the RNA bam
## PBS configure 
#PBS -N get_RNA_bw
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or154_Or163
  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_RNA.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_RNA.bam  > ./barcode_RNA.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_RNA.bedGraph barcode_RNA.norm.bedGraph &> barcode_RNA.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_RNA.norm.bedGraph > barcode_RNA.norm.sorted.bedGraph
  bedGraphToBigWig barcode_RNA.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_RNA.norm.bw
# subset the ATAC bam
## PBS configure 
#PBS -N get_ATAC_bw
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=32G
  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_ATAC.bam  > ./barcode_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_ATAC.bedGraph barcode_ATAC.norm.bedGraph &> barcode_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_ATAC.norm.bedGraph > barcode_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_ATAC.norm.bw 

## PBS configure 
#PBS -N Atacworks
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=32G
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or154_Or163
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate python37
honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
python $atacworks/scripts/main.py denoise \
    --noisybw barcode_ATAC.norm.bw \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    # --regions C234.genome_intervals.bed \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml

# Fig3H:single promoter  (Or25-27) 

# Umap plot:
library(scCustomize)
library(Nebulosa)

pdf('./00_Figure/Fig3/Fig3H-Or25_27-FeaturePlot.pdf', width=18, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("Or25"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)
p2<- FeaturePlot(ORN, features = c("Or26"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)
p3<- FeaturePlot(ORN, features = c("Or27"),cols=c("lightgrey", "#5A4CE1"), max.cutoff =3, order=TRUE,)
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("Or25","Or26","Or27"),
                                custom_palette = BlueAndRed())
p1|p2|p3|p12
dev.off()

obj<- subset(ORN,idents=c("24"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("Or25","Or26","Or27")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
colnames(ORN_matrix)<- c("Or25","Or26","Or27")
ORN_matrix<- as.data.frame(ORN_matrix)
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)
pdf("./00_Figure/Fig3/Fig3I-coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or25,ORN_matrix$Or26,ORN_matrix$Or27,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color =  colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or27,ORN_matrix$Or26,ORN_matrix$Or25,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
pheatmap(t(ORN_matrix),
     cluster_cols = T,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()


# distinguish the coexp and the one OR exp cells
Or27_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or27>0&ORN_matrix$Or26==0&ORN_matrix$Or25==0),])
Or26_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or27==0&ORN_matrix$Or26>0&ORN_matrix$Or25==0),])
Or25_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or25>0),])
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or27>0&ORN_matrix$Or26>0&ORN_matrix$Or25>0),])

# Or25>0 in all cells!
coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or25,coexp_matrix$Or26,coexp_matrix$Or27,decreasing=TRUE),]
bk <- c(seq(0,3,by=0.1))
pdf("./00_Figure/Fig3/Fig3H-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
pheatmap(t(coexp_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colorRampPalette(colors = c("white","#CC0000"))(length(bk)),
     annotation_legend = TRUE,
     show_rownames=T,breaks=bk,
     #annotation_col = label_pheatmap,
     #annotation_colors = ann_colors,
     show_colnames=F
  )
dev.off()



# save the coexp barcode for track 
library(tidyr)
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or25_Or26_Or27/barcode_coexp.txt",row.name=F,col.names=F)

write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or25_Or26_Or27/raw_barcode_coexp.txt",row.name=F,col.names=F)

obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
set.seed(5678)
random_barcode1<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
random_barcode2<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
barcode<- data.frame(barcode=random_barcode1)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or25_Or26_Or27/random_barcode3.txt",row.name=F,col.names=F)
barcode<- data.frame(barcode=random_barcode2)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or25_Or26_Or27/random_barcode4.txt",row.name=F,col.names=F)

sed -i 's/"//g' *.txt

# subset the RNA bam
## PBS configure 
#PBS -N get_RNA_bw
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or25_Or26_Or27
  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes barcode_coexp.txt --cores 20 --out-bam barcode_coexp_RNA.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_coexp_RNA.bam  > ./barcode_coexp_RNA.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_coexp_RNA.bedGraph barcode_coexp_RNA.norm.bedGraph &> barcode_coexp_RNA.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_coexp_RNA.norm.bedGraph > barcode_coexp_RNA.norm.sorted.bedGraph
  bedGraphToBigWig barcode_coexp_RNA.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_coexp_RNA.norm.bw

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes random_barcode1.txt --cores 20 --out-bam random_barcode1_RNA.bam
  bedtools  genomecov  -bg -split -ibam ./random_barcode1_RNA.bam  > ./random_barcode1_RNA.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl random_barcode1_RNA.bedGraph random_barcode1_RNA.norm.bedGraph &> random_barcode1_RNA.norm.bedGraph.log
  sort -k1,1 -k2,2n random_barcode1_RNA.norm.bedGraph > random_barcode1_RNA.norm.sorted.bedGraph
  bedGraphToBigWig random_barcode1_RNA.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt random_barcode1_RNA.norm.bw

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes random_barcode2.txt --cores 20 --out-bam random_barcode2_RNA.bam
  bedtools  genomecov  -bg -split -ibam ./random_barcode2_RNA.bam  > ./random_barcode2_RNA.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl random_barcode2_RNA.bedGraph random_barcode2_RNA.norm.bedGraph &> random_barcode2_RNA.norm.bedGraph.log
  sort -k1,1 -k2,2n random_barcode2_RNA.norm.bedGraph > random_barcode2_RNA.norm.sorted.bedGraph
  bedGraphToBigWig random_barcode2_RNA.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt random_barcode2_RNA.norm.bw

# subset the ATAC bam
## PBS configure 
#PBS -N get_ATAC_bw
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=32G
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or25_Or26_Or27

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode_coexp.txt --cores 20 --out-bam barcode_coexp_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_coexp_ATAC.bam  > ./barcode_coexp_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_coexp_ATAC.bedGraph barcode_coexp_ATAC.norm.bedGraph &> barcode_coexp_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_coexp_ATAC.norm.bedGraph > barcode_coexp_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_coexp_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_coexp_ATAC.norm.bw

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes random_barcode1.txt --cores 20 --out-bam random_barcode1_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./random_barcode1_ATAC.bam  > ./random_barcode1_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl random_barcode1_ATAC.bedGraph random_barcode1_ATAC.norm.bedGraph &> random_barcode1_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n random_barcode1_ATAC.norm.bedGraph > random_barcode1_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig random_barcode1_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt random_barcode1_ATAC.norm.bw

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes random_barcode2.txt --cores 20 --out-bam random_barcode2_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./random_barcode2_ATAC.bam  > ./random_barcode2_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl random_barcode2_ATAC.bedGraph random_barcode2_ATAC.norm.bedGraph &> random_barcode2_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n random_barcode2_ATAC.norm.bedGraph > random_barcode2_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig random_barcode2_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt random_barcode2_ATAC.norm.bw




## PBS configure 
#PBS -N Atacworks
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=32G
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig3/Or154_Or163
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate python37
honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
python $atacworks/scripts/main.py denoise \
    --noisybw barcode_ATAC.norm.bw \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    # --regions C234.genome_intervals.bed \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml






