## phylogenetic tree only honeybee 

1. OR 
# New gtf OR sequence similarity:

```
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(OR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
rownames(Group_info)<-Group_info$gene_name
groupInfo <- split(row.names(Group_info), Group_info$seqnames)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=16,height=16)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()

# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/all_OR_fasta_aa_similarity.csv");

```

2. GR 

```
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(GR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/GR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=10,height=10)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()

# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/GR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/all_GR_fasta_aa_similarity.csv");

```

3. IR 


```
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(IR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/IR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=10,height=10)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()
# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/IR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/all_IR_fasta_aa_similarity.csv");

```

## phylogenetic tree honeybee, fly and  

1. get sequence 

```
sed -i '/#/d;' AaegL5-ChemoreceptorPeptides.fasta 
```



```

library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
AaegL5_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/AaegL5-ChemoreceptorPeptides.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
names_in_myXStringSet <- names(AaegL5_fasta)
#select names from your XStringSet
Aaeg_OR <- AaegL5_fasta[names_in_myXStringSet[grep("AaegOr",names_in_myXStringSet)],]
Aaeg_IR <- AaegL5_fasta[names_in_myXStringSet[grep("AaegIr",names_in_myXStringSet)],]
Aaeg_GR <- AaegL5_fasta[names_in_myXStringSet[grep("AaegGr",names_in_myXStringSet)],]


#get fly OR sequence;
library(biomaRt)
mart <- useMart("ensembl","dmelanogaster_gene_ensembl")
fly<-readRDS("/md01/liyh526/project/Fly-cooperate/7.18run/WNN_ORN_integrated_antenna.rds")
allgene<-rownames(fly)
OR<-grep("^Or[0-9]",allgene,value=T)
IR<-grep("^Ir[0-9]",allgene,value=T)
GR<-grep("^Gr[0-9]",allgene,value=T)

all<-c(OR,IR,GR,"Orco")
OR_ppseqs <- getSequence(id = c(OR,"Orco"),
                          type="external_gene_name",
                          seqType="peptide",
                          mart = mart) 
GR_ppseqs <- getSequence(id = GR,
                          type="external_gene_name",
                          seqType="peptide",
                          mart = mart) 
IR_ppseqs <- getSequence(id = IR,
                          type="external_gene_name",
                          seqType="peptide",
                          mart = mart)                                                     
OR_ppseqs <- OR_ppseqs[OR_ppseqs$peptide!="Sequence unavaliable",]
IR_ppseqs <- IR_ppseqs[IR_ppseqs$peptide!="Sequence unavaliable",]
GR_ppseqs <- GR_ppseqs[GR_ppseqs$peptide!="Sequence unavaliable",]
OR_ppseqs$external_gene_name <- make.unique(OR_ppseqs$external_gene_name)
GR_ppseqs$external_gene_name <- make.unique(GR_ppseqs$external_gene_name)
IR_ppseqs$external_gene_name <- make.unique(IR_ppseqs$external_gene_name)
exportFASTA(OR_ppseqs,"./fly_OR_pep.fasta")
exportFASTA(IR_ppseqs,"./fly_IR_pep.fasta")
exportFASTA(GR_ppseqs,"./fly_GR_pep.fasta")

fly_OR_fasta<-readAAStringSet("./fly_OR_pep.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
fly_GR_fasta<-readAAStringSet("./fly_GR_pep.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
fly_IR_fasta<-readAAStringSet("./fly_IR_pep.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)


All_OR_fasta<- c(OR_fasta,fly_OR_fasta,Aaeg_OR)
All_GR_fasta<- c(GR_fasta,fly_GR_fasta,Aaeg_GR)
All_IR_fasta<- c(IR_fasta,fly_IR_fasta,Aaeg_IR)

# create the chemoreceptor info 
Amel <- names(c(OR_fasta,GR_fasta,IR_fasta))
fly <- names(c(fly_OR_fasta,fly_GR_fasta,fly_IR_fasta))
Aaeg <- names(c(Aaeg_OR,Aaeg_IR,Aaeg_GR))

chemo_info <- data.frame(gene=c(Amel,fly,Aaeg),
species=c(rep("Amel",length(Amel)),rep("Dmel",length(fly)),rep("Aaeg",length(Aaeg))))
rownames(chemo_info) <- chemo_info$gene
groupInfo <- split(chemo_info$gene, chemo_info$species)

# For OR 
aln <- muscle::muscle(All_OR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/three_species_OR_sequence_protein_similarity-tree.pdf",width=20,height=20)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()

# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/cross_species_OR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/cross_species_OR_fasta_aa_similarity.csv");



# For IR 
aln <- muscle::muscle(All_IR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/three_species_IR_sequence_protein_similarity-tree.pdf",width=20,height=20)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()

# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/cross_species_IR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/cross_species_IR_fasta_aa_similarity.csv");


# For GR 
aln <- muscle::muscle(All_GR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/three_species_GR_sequence_protein_similarity-tree.pdf",width=20,height=20)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()
# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/cross_species_GR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/cross_species_GR_fasta_aa_similarity.csv");

```





