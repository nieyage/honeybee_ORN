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
DefaultAssay(ORN)<-"raw_RNA"

# Fig3A:two promoter  (Or128 and Or129) 
# Umap plot:
library(scCustomize)
library(Nebulosa)

pdf('./00_Figure/FigS3/FigS3A-Or128_129-FeaturePlot.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("LOC102656904"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("8:Or128")
p2<- FeaturePlot(ORN, features = c("LOC102656221"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("8:Or129F")
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("LOC102656904","LOC102656221"),
                                custom_palette = BlueAndRed())
p1|p2|p12
dev.off()

obj<- subset(ORN,idents=c("8"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC102656904","LOC102656221")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
colnames(ORN_matrix)<- c("Or129F","Or128")
ORN_matrix<- as.data.frame(ORN_matrix)
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/FigS3/FigS3B-coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or128,ORN_matrix$Or129,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or129,ORN_matrix$Or128,decreasing=T),]
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
Or129F_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or129F>0&ORN_matrix$Or128==0),])
Or128_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or129F==0&ORN_matrix$Or128>0),])
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or129F>0&ORN_matrix$Or128>0),])

# plot the heatmap 
coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or128,coexp_matrix$Or129F,decreasing=TRUE),]

Or128_matrix<- ORN_matrix[Or128_barcode,]
Or128_matrix<-Or128_matrix[order(Or128_matrix$Or128,Or128_matrix$Or129F,decreasing=TRUE),]
Or129_matrix<- ORN_matrix[Or129F_barcode,]
Or129_matrix<-Or129_matrix[order(Or129_matrix$Or128,Or129_matrix$Or129F,decreasing=TRUE),]
ORN_matrix<- rbind(coexp_matrix,Or128_matrix,Or129_matrix)
label_pheatmap<- data.frame(label=c(rep("coexp",length(coexp_barcode)),
     rep("only_Or128",length(Or128_barcode)),
     rep("only_Or129",length(Or129F_barcode))))
col <- c("#8A9FD1" ,"#C06CAB","#FEE500")

names(col)<-c("only_Or128","only_Or129","coexp")
ann_colors= list(label = col)
rownames(label_pheatmap)<- rownames(ORN_matrix)
pdf("./00_Figure/FigS3/FigS3B-Or128_129_coexp-heatmap-SCT-split_coexp.pdf",height=3,width=12)
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
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or128,coexp_matrix$Or129F,decreasing=TRUE),]
bk <- c(seq(0,2,by=0.1))
coexp_matrix<-coexp_matrix[,c(2,1)]
pdf("./00_Figure/FigS3/FigS3F-Or128_129-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
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
     
write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or128_Or129/raw_barcode_coexp.txt",row.name=F,col.names=F)



# save the coexp barcode for track 
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or128_Or129/barcode.txt",row.name=F,col.names=F)
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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or128_Or129
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
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or128_Or129

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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or128_Or129
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate python37
honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
python $atacworks/scripts/main.py denoise \
    --noisybw barcode_ATAC.norm.bw \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml

# S3-2：SP example
pdf('./00_Figure/FigS3/FigS3-2-A-Or1_2-FeaturePlot.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("LOC725205"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("17:Or1")
p2<- FeaturePlot(ORN, features = c("LOC100577787"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("17:Or3")
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("LOC725205","LOC100577787"),
                                custom_palette = BlueAndRed())
p1|p2|p12
dev.off()

obj<- subset(ORN,idents=c("17"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC725205","LOC100577787")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
colnames(ORN_matrix)<- c("Or3","Or1")
ORN_matrix<- as.data.frame(ORN_matrix)
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)
ORN_matrix<- ORN_matrix[,c(2,1)]
pdf("./00_Figure/FigS3/FigS3-2-B-OR1_3-coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or1,ORN_matrix$Or3,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or3,ORN_matrix$Or1,decreasing=T),]
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

# coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$cell_group<-as.character(obj$cell_group)
for (i in 1:length(obj$cell_group)){
  if(obj$cell_group[i]!="17"){obj$cell_group[i]="other"}
    }
Idents(obj)<-obj$cell_group
DefaultAssay(obj)<-"peaks_ORN_subcluster"
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
obj_features<- gene
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/FigS3/FigS3C-Or1_3_trackplot-singlepromoter.pdf",width=10,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()
Or1_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or1>0&ORN_matrix$Or3==0),])
Or3_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or1==0&ORN_matrix$Or3>0),])
Or1_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or1>0),])
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or1>0&ORN_matrix$Or3>0),])

coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or1>0&ORN_matrix$Or3>0),])


coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or1,coexp_matrix$Or3,decreasing=TRUE),]
bk <- c(seq(0,2,by=0.1))
pdf("./00_Figure/FigS3/FigS3-Or1_3-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
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
     
write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3/new/raw_barcode_coexp.txt",row.name=F,col.names=F)
# save the coexp barcode for track 
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3/new/barcode.txt",row.name=F,col.names=F)
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
set.seed(110)
random_barcode1<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
random_barcode2<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
barcode<- data.frame(barcode=random_barcode1)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3/new/random_barcode1.txt",row.name=F,col.names=F)
barcode<- data.frame(barcode=random_barcode2)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3/new/random_barcode2.txt",row.name=F,col.names=F)




# save the coexp barcode for track 
barcode<- data.frame(barcode=Or1_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3/barcode.txt",row.name=F,col.names=F)
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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3
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
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or1_Or3

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_ATAC.bam  > ./barcode_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_ATAC.bedGraph barcode_ATAC.norm.bedGraph &> barcode_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_ATAC.norm.bedGraph > barcode_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_ATAC.norm.bw 


DefaultAssay(ORN)<-"raw_RNA"
pdf('./00_Figure/FigS3/FigS3A-Or12_14-FeaturePlot.pdf', width=18, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("Or12"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("19:Or12")
p2<- FeaturePlot(ORN, features = c("Or13"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("19:Or13")
p3<- FeaturePlot(ORN, features = c("Or14"),cols=c("lightgrey", "#5A4CE1"), max.cutoff =3, order=TRUE,)+ggtitle("19:Or14")

p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("Or12","Or13","Or14"),
                                custom_palette = BlueAndRed())
p1|p2|p3|p12
dev.off()
obj<- subset(ORN,idents=c("19"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("Or12","Or13","Or14")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
ORN_matrix<- as.data.frame(ORN_matrix)
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/FigS3/FigS3-Or12_14coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or12,ORN_matrix$Or13,ORN_matrix$Or14,decreasing=T),]
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
Or12_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or12>0&ORN_matrix$Or13==0&ORN_matrix$Or14==0),])
Or13_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or12==0&ORN_matrix$Or13>0&ORN_matrix$Or14==0),])
Or14_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or12==0&ORN_matrix$Or13==0&ORN_matrix$Or14>0),])
Or12_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or12>0),])
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or12>0&ORN_matrix$Or13>0&ORN_matrix$Or14>0),])


coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or12,coexp_matrix$Or13,coexp_matrix$Or14,decreasing=TRUE),]
bk <- c(seq(0,2,by=0.1))
pdf("./00_Figure/FigS3/FigS3_Or12_14-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
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
barcode<- data.frame(barcode=Or12_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14/barcode.txt",row.name=F,col.names=F)
sed -i 's/"//g' barcode.txt



# save the coexp barcode for track 
library(tidyr)
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14/new/barcode_coexp.txt",row.name=F,col.names=F)

write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14/new/raw_barcode_coexp.txt",row.name=F,col.names=F)

obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
set.seed(110)
random_barcode1<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
random_barcode2<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
barcode<- data.frame(barcode=random_barcode1)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14/new/random_barcode1.txt",row.name=F,col.names=F)
barcode<- data.frame(barcode=random_barcode2)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14/new/random_barcode2.txt",row.name=F,col.names=F)

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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14
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
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or12_Or13_Or14

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_ATAC.bam  > ./barcode_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_ATAC.bedGraph barcode_ATAC.norm.bedGraph &> barcode_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_ATAC.norm.bedGraph > barcode_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_ATAC.norm.bw 

# Or169 and Or170
DefaultAssay(ORN)<-"raw_RNA"
pdf('./00_Figure/FigS3/FigS3A-Or169_170-FeaturePlot.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("LOC102654841"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("29:Or169")
p2<- FeaturePlot(ORN, features = c("Or170"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("29:Or170")
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("LOC102654841","Or170"),
                                custom_palette = BlueAndRed())
p1|p2|p12
dev.off()

obj<- subset(ORN,idents=c("29"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC102654841","Or170")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
ORN_matrix<- as.data.frame(ORN_matrix)
colnames(ORN_matrix)<-c("Or169","Or170")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/FigS3/FigS3-Or169_170coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or169,ORN_matrix$Or170,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# distinguish the coexp and the one OR exp cells
Or169_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or169>0&ORN_matrix$Or170==0),])
Or170_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or169==0&ORN_matrix$Or170>0),])
Or169_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or169>0),])
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or169>0&ORN_matrix$Or170>0),])

# save the coexp barcode for track 
barcode<- data.frame(barcode=Or169_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170/barcode.txt",row.name=F,col.names=F)
sed -i 's/"//g' barcode.txt


coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or169,coexp_matrix$Or170,decreasing=TRUE),]
bk <- c(seq(0,2,by=0.1))
pdf("./00_Figure/FigS3/FigS3-Or169-170-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
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
     
write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170/new/raw_barcode_coexp.txt",row.name=F,col.names=F)

obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
set.seed(110)
random_barcode1<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
random_barcode2<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
barcode<- data.frame(barcode=random_barcode1)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170/new/random_barcode1.txt",row.name=F,col.names=F)
barcode<- data.frame(barcode=random_barcode2)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170/new/random_barcode2.txt",row.name=F,col.names=F)
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170/new/barcode_coexp.txt",row.name=F,col.names=F)


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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170
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
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or169_Or170

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_ATAC.bam  > ./barcode_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_ATAC.bedGraph barcode_ATAC.norm.bedGraph &> barcode_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_ATAC.norm.bedGraph > barcode_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_ATAC.norm.bw 


# Or39 and Or40
DefaultAssay(ORN)<-"raw_RNA"
pdf('./00_Figure/FigS3/FigS3A-Or39_40-FeaturePlot.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("LOC100577101"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("10:Or39")
p2<- FeaturePlot(ORN, features = c("LOC100577068"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("10:Or40")
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("LOC100577101","LOC100577068"),
                                custom_palette = BlueAndRed())
p1|p2|p12
dev.off()

obj<- subset(ORN,idents=c("10"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC100577101","LOC100577068")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
ORN_matrix<- as.data.frame(ORN_matrix)
ORN_matrix<- ORN_matrix[,c(2,1)]
colnames(ORN_matrix)<-c("Or39","Or40")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/FigS3/FigS3-Or39_40coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or39,ORN_matrix$Or40,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# distinguish the coexp and the one OR exp cells
Or39_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or39>0&ORN_matrix$Or40==0),])
Or40_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or39==0&ORN_matrix$Or40>0),])
Or39_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or39>0),])
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or39>0&ORN_matrix$Or40>0),])


coexp_matrix<- ORN_matrix[coexp_barcode,]
coexp_matrix<-coexp_matrix[order(coexp_matrix$Or39,coexp_matrix$Or40,decreasing=TRUE),]
bk <- c(seq(0,2,by=0.1))
pdf("./00_Figure/FigS3/FigS3-Or39_40-coexp-heatmap-SCT-only_coexp.pdf",height=3,width=12)
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
     
write.table(coexp_barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40/new/raw_barcode_coexp.txt",row.name=F,col.names=F)
# save the coexp barcode for track 
barcode<- data.frame(barcode=coexp_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40/new/barcode.txt",row.name=F,col.names=F)
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
set.seed(110)
random_barcode1<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
random_barcode2<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
barcode<- data.frame(barcode=random_barcode1)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40/new/random_barcode1.txt",row.name=F,col.names=F)
barcode<- data.frame(barcode=random_barcode2)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40/new/random_barcode2.txt",row.name=F,col.names=F)




# save the coexp barcode for track 
barcode<- data.frame(barcode=Or39_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40/barcode.txt",row.name=F,col.names=F)
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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40
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
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or39_Or40

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_ATAC.bam  > ./barcode_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_ATAC.bedGraph barcode_ATAC.norm.bedGraph &> barcode_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_ATAC.norm.bedGraph > barcode_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_ATAC.norm.bw 

# Or132 and Or133
# Or132 and Or133N
DefaultAssay(ORN)<-"raw_RNA"
pdf('./00_Figure/FigS3/FigS3A-Or132_133N-FeaturePlot.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<- FeaturePlot(ORN, features = c("LOC102653615"),cols=c("lightgrey", "#4DAE49"), max.cutoff =3, order=TRUE,)+ggtitle("6:Or132")
p2<- FeaturePlot(ORN, features = c("LOC102653695"),cols=c("lightgrey", "#E31A1C"), max.cutoff =3, order=TRUE,)+ggtitle("6:Or133N")
p12<- Plot_Density_Joint_Only(seurat_object = ORN, 
                                features = c("LOC102653615","LOC102653695"),
                                custom_palette = BlueAndRed())
p1|p2|p12
dev.off()

obj<- subset(ORN,idents=c("6"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC102653615","LOC102653695")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
ORN_matrix<- as.data.frame(ORN_matrix)
colnames(ORN_matrix)<-c("Or132","Or133N")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/FigS3/FigS3-Or132_133coexp-heatmap-SCT.pdf",height=3,width=12)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or132,ORN_matrix$Or133N,decreasing=T),]
pheatmap(t(ORN_matrix),
     cluster_cols = F,
     cluster_rows = F,
     border=F,
     color = colors_list[10:100],
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# distinguish the coexp and the one OR exp cells
Or132_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or132>0&ORN_matrix$Or133N==0),])
Or133N_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or132==0&ORN_matrix$Or133N>0),])
Or132_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or132>0),])

# save the coexp barcode for track 
barcode<- data.frame(barcode=Or132_barcode)
barcode_sample<-separate(barcode,"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or132_Or133/barcode.txt",row.name=F,col.names=F)
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
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or132_Or133
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
source /public/home/nieyg/.bash_profile
cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/Or132_Or133

  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes barcode.txt --cores 20 --out-bam barcode_ATAC.bam
  bedtools  genomecov  -bg -split -ibam ./barcode_ATAC.bam  > ./barcode_ATAC.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl barcode_ATAC.bedGraph barcode_ATAC.norm.bedGraph &> barcode_ATAC.norm.bedGraph.log
  sort -k1,1 -k2,2n barcode_ATAC.norm.bedGraph > barcode_ATAC.norm.sorted.bedGraph
  bedGraphToBigWig barcode_ATAC.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt barcode_ATAC.norm.bw 





# single OR:Or10 and Or11 


RNAfold --noconv -p < Or11.fa > Or10.res
/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl Or11_ss.ps Or11_dp.ps > Or11_rss.ps
magick convert -density 300 Or11_rss.ps Or11_rss.pdf # 生成二级结构pdf图

mv /md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/Or11_rss.pdf /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/hairpin_loss
mv /md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/Or10_rss.pdf /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/hairpin_loss


>Or10	(((((...........)))))........ ( -0.10)


>Or10
TAATTTTATTCCTTGTCATAATTCGTACCTAAACATCATGAAACTATTCG



RNAfold --noconv -p < Or10.fa > Or10.res
/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl Or10_ss.ps Or10_dp.ps > Or10_rss.ps
magick convert -density 300 Or10_rss.ps Or10_rss.pdf # 生成二级结构pdf图




# coexp_RT:Or64 Or150

RNAfold --noconv -p < Or64.fa > Or64.res
/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl Or64_ss.ps Or64_dp.ps > Or64_rss.ps
magick convert -density 300 Or64_rss.ps Or64_rss2.pdf # 生成二级结构pdf图

RNAfold --noconv -p < Or150.fa > Or150.res
/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl Or150_ss.ps Or150_dp.ps > Or150_rss.ps
magick convert -density 300 Or150_rss.ps Or150_rss.pdf # 生成二级结构pdf图


mv /md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RNAfold/Or150_rss.pdf /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/hairpin_loss
mv /md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RNAfold/Or64_rss.pdf /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/hairpin_loss



# coexp_NRT:Or152(LOC102655559) and Or67(LOC100578045)
paste ID result |sort -k1,1 

>LOC100578045	((((((((............)))))))). ( -1.90)
>LOC100578045
ATATTAAAACAAATATTTTAACGATATAGAAGAATCTAAAAATTTATGTTA

cd /data/R02/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RNAfold/LOC100578045_rss.pdf




RNAfold --noconv -p < Or67.fa > Or67.res
/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl Or67_ss.ps Or67_dp.ps > Or67_rss.ps
magick convert -density 300 Or67_rss.ps Or67_rss.pdf # 生成二级结构pdf图



>LOC102655559	((((((((................))))))))....... ( -4.40)

>LOC102655559
ATAATATAATATATTAATCTTTTTTATTAAATGATGATTAATAATTGATA
RNAfold --noconv -p < Or152.fa > Or152.res
/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl Or152_ss.ps Or152_dp.ps > Or152_rss.ps
magick convert -density 300 Or152_rss.ps Or152_rss.pdf # 生成二级结构pdf图


mv /md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RNAfold/Or152_rss.pdf /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/hairpin_loss
mv /md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RNAfold/Or67_rss.pdf /data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS3/hairpin_loss





Idents(ORN) <- ORN$cell_group
obj <- subset(ORN,idents=c("19"))
DefaultAssay(obj)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
set.seed(110)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj <- subset(ORN,cells=c(obj_barcode,random_barcode))
obj$cell_group<-as.character(obj$cell_group)
for (i in 1:length(obj$cell_group)){
  if(obj$cell_group[i]!="19"){obj$cell_group[i]="other"}
}
Idents(obj)<-obj$cell_group

pdf("./00_Figure/FigS3/Or12_Or13_Or14_fragmentplot_C19.pdf",width=16,height=12)
CoveragePlot(
  object = obj,
  region = "Group2-10013251-10019530",
  window = 150,
  extend.upstream = 700,
  annotation = TRUE,
  extend.downstream = 700,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
dev.off()


CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 30,
    links=F
  )


# 提取 Or10 genebody and downstream TES region and label the hairpin location 

>Or10
TAATTTTATTCCTTGTCATAATTCGTACCTAAACATCATGAAACTATTCG






gtgtataattttttcctccagtatcctttttatttattttatatttaccaaAGAATAGTAATAATTGTAAACTTTATAAACAATACAATGATATTCtggagatttataaattatagtcaTTATgctagaatttataaataataaattctctcacaattatcataaatatatattttcatgtaagatattttttcaaatattttataattttctattataagtttgatttttgaaaaattatcagaatAATCATATTCGAGGTATCGAGTATTTAaagtttacaaaaaaataaattcatcaataaaataattttctataaaatgtaatttatgctAAAGATGTATACAaccaattttcataaaaatttaaattgtttcgcAATCTTTTCGTGTTTTACGATTATTCTCATAACAGTAACAGAATATCCTGTAACATTTCATATtcttaaatacaataaaattgaaatgagtACAGTAATGTTAGCAATTACGAATTATgacaatagaataattttattccttgtCATAATTCGTACCTAAACATCATGAAactattcgttttttttaaatataaacacatagtcaaaatatatttctactttaccaaaaagaagatttagtcttttattttttatttcatttaacttGATCTCAGAAAtctgataaaattttgcagaaaagaaaagaagaaaattatcatgaaaaaattaaaatgttgattTCGCATATCaaacaaatatatgattttagtgtgataaataaagaaatttaatataatattgaaattgaagaaatattacttatttatcattaattgaaaatctgtGAATTTTCatcaggaatatatatatatataccatagctataataatcattgtttcgtattaaattattggattCTTTCATTCTAGAAAAATAGATCTTTAAATCGAAAAACACTAAGATTGCTACACATTTTTTaaggaataatatatactGAATCATGAGATTATTTTTGCGAT





