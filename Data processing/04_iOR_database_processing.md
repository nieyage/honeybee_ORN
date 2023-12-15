## get infomation from the iorbase
1. get the Apis OR sequence and mapping with our dotplot feature 

```
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
dotplot_data<-read.csv("~/project/honeybee/honebee-latest-Version//05_ORN_cluster/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

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
dotplot_feature_OR2<-c(Orco,dotplot_feature)
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature_OR2,]
iOR_fasta<-readAAStringSet("/data/R02/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/iORbase_ApMe_seq.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
total<- c(dotplot_feature_fa,iOR_fasta)

aln <- muscle::muscle(total)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
pdf("./iOR_protein_mapping_our_Data.pdf",width=16,height=16)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none") + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()
# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;
data<-data.frame()
for (i in 1:122){
	row=i;
	for (j in 123:ncol(dist_percent)){
	col=j;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}}
}
data<-data[order(data$rate,decreasing=TRUE),]
ID_trans<-data[!duplicated(data$gene1),]
ID_trans<-data[!duplicated(data$gene2),]
ID_trans<-ID_trans[ID_trans$rate>70,]
write.csv(ID_trans,"iOR_baseID_2_ourID.csv")
```

2. get the Apis OR pdb file 

```
#!/bin/bash

cat iORbase_id.txt | while read line
do
#echo $line
arr=($line)
id=${arr[0]}
echo  https://www.iorbase.com/protein_download/?change_name=$id -O $id.pdb >> url.txt
done

```
# %s/\r//g modify url.txt
sed -i -e 's/^/wget /' url.txt
bash url.txt

3. calculate the structural similarity by RMSD 
***fly coreceptor 
```
ls *.pdb > pdb_list.txt
for file1 in $(<pdb_list.txt)
do
 for file2 in $(<pdb_list.txt)
 do 
  python ../RMSD.py $file1 $file2 >> RMSD_result.txt
 done
done

```

# back 2 R 
```
RMSD <- read.table("RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("pdb_list.txt")
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


```

***fly and honybee coreceptor 
```
ls *.pdb > pdb_list.txt
for file1 in $(<pdb_list.txt)
do
 for file2 in $(<pdb_list.txt)
 do 
  python ../RMSD.py $file1 $file2 >> RMSD_result.txt
 done
done
```

```
RMSD <- read.table("RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("pdb_list.txt")
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
label_pheatmap<- data.frame(species=c(rep("D.melanogaster",5),rep("Apis mellifera",4)))
rownames(label_pheatmap) <- colnames(data)
data[data==0]<- 0.5
data<- log2(data)
ann_colors<-list(
species = c("D.melanogaster"="#7985B6", "Apis mellifera"="#C7B6E1"))

pdf("./honeybee_OR_feature_RMSD.pdf",width=7,height=5)
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


```




*** honebee OR 
```
#!/bin/bash
nLine=`wc -l iOR_baseID_2_ourID.txt | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        ourID=`awk '{print $1}' iOR_baseID_2_ourID.txt | sed -n "${i}p"`
        iORbaseID=`awk '{print $2}' iOR_baseID_2_ourID.txt | sed -n "${i}p"`
        cp ./pdb_file/${iORbaseID}.pdb ./rename_pdb/
        mv ./rename_pdb/${iORbaseID}.pdb ./rename_pdb/${ourID}.pdb
done
cd rename_pdb
ls *.pdb > pdb_list.txt
for file1 in $(<pdb_list.txt)
do
 for file2 in $(<pdb_list.txt)
 do 
  python ../RMSD.py $file1 $file2 >> RMSD_result.txt
 done
done

```

# back 2 R 
```
RMSD <- read.table("RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("pdb_list.txt")
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

dotplot_feature_fa<- all_receptor_gene_fasta[ID,]
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
tree<-as.phylo(clust)
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
groupInfo <- split(Group_info$gene_name, Group_info$seqnames)
data.tree <- groupOTU(tree, groupInfo)

pdf("./OR_sequence_protein_similarity-tree.pdf",width=25,height=10)
ggtree(data.tree,ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab(size=3) + theme(legend.position = "right")
dev.off()
m<-ggtree(data.tree,ladderize = FALSE, branch.length = "none",aes(color=group))
gene_order<-na.omit(m$data[order(m$data$y),]$label)
gene_order<-as.character(gene_order)

# plot heatmap 
library(pheatmap)
data[data==0]<- 0.5
data<- log2(data)
pdf("./honeybee_OR_feature_RMSD.pdf",width=11,height=10)
pheatmap(data[rev(gene_order),rev(gene_order)],
      cluster_cols = F,
      cluster_rows = F,
      color = colorRampPalette(c("#1C214F","#3082BD","#F1ECEC"))(100),
      #annotation_col = barcode_label_pheatmap,
      #annotation_colors = ann_colors,
      #annotation_row = barcode_label_pheatmap,
      annotation_legend = TRUE,
      show_rownames=T,
      show_colnames=T
 )
dev.off()

```
