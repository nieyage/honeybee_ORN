#Step1: ORN cluster by OR gene 
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

gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']
gene.coords <- gene.coords[gene.coords$gene_name %in% all_receptor_gene]
gene.coords <- gene.coords[gene.coords$type == "transcript"]
gene.coords <- as.data.frame(gene.coords)
strand_p <- gene.coords[gene.coords$strand == "+",]
strand_n <- gene.coords[gene.coords$strand == "-",]

# for strand "+" OR gene 
strand_p <- strand_p[!duplicated(strand_p$gene_name),]
strand_p <- strand_p[order(strand_p$seqnames,strand_p$start,decreasing=F),]
OR_TES_location<- data.frame()
for (i in 1:nrow(strand_p)){
	OR1 <- strand_p[i,]
	if(i!=nrow(strand_p)){
	OR2 <- strand_p[i+1,]
	if(OR1$seqnames==OR2$seqnames){
		OR1_TES_start <- OR1$end+1;
		intersections <- abs(OR1$end-OR2$start);
		if(intersections< 1000){
			OR1_TES_end <- OR2$start-1;
			if(OR1_TES_end<OR1_TES_start){OR1_TES_end=OR1_TES_start+1000}else{OR1_TES_end=OR1_TES_end}
		}else{
			OR1_TES_end <- OR1_TES_start+1000;
		}
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_start,end=OR1_TES_end,strand=OR1$strand)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}else{
		OR1_TES_start <- OR1$end+1;
		OR1_TES_end <- OR1_TES_start+1000;
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_start,end=OR1_TES_end,strand=OR1$strand)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}
}else{
	OR1_TES_start <- OR1$end+1;
	OR1_TES_end <- OR1_TES_start+1000;
	OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_start,end=OR1_TES_end,strand=OR1$strand)
	OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
}
}
OR_TES_location_strand_p<- OR_TES_location;

# for strand "-" OR gene 
strand_n <- strand_n[!duplicated(strand_n$gene_name),]
strand_n <- strand_n[order(strand_n$seqnames,strand_n$start,decreasing=T),]

OR_TES_location<- data.frame()
for (i in 1:nrow(strand_n)){
	OR1 <- strand_n[i,]
	if(i!=nrow(strand_n)){
	OR2 <- strand_n[i+1,]
	if(OR1$seqnames==OR2$seqnames){
		OR1_TES_start <- OR1$start-1;
		intersections <- abs(OR1$start-OR2$end);
		if(intersections< 1000){
			OR1_TES_end <- OR2$end+1;
		}else{
			OR1_TES_end <- OR1_TES_start-1000;
		}
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_end,end=OR1_TES_start,strand=OR1$strand)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}else{
		OR1_TES_start <- OR1$start-1;
		OR1_TES_end <- OR1_TES_start-1000;
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_end,end=OR1_TES_start,strand=OR1$strand)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}
}else{
	OR1_TES_start <- OR1$start-1;
	OR1_TES_end <- OR1_TES_start-1000;
	OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_end,end=OR1_TES_start,strand=OR1$strand)
	OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
}
}
OR_TES_location$strand<- "-"
OR_TES_location<- rbind(OR_TES_location,OR_TES_location_strand_p)
OR_TES_location$seqnames<- as.character(OR_TES_location$seqnames)

# save the OR TES sequence location as bed file (one OR one file)
# read through OR gene 
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
write.csv(as.data.frame(dotplot_feature),"dotplot_feature_TES_type.csv")

OR_type<- read.csv("./13_TES_signal/dotplot_feature_TES_type.csv")
NRT_OR<- OR_type[OR_type$type=="NRT",2]
RT_OR<- OR_type[OR_type$type=="RT",2]
OR_TES_location$V1<- 1;
NRT_OR_TES_location<- OR_TES_location[OR_TES_location$OR%in%NRT_OR,]
RT_OR_TES_location<- OR_TES_location[OR_TES_location$OR%in%RT_OR,]
	write.table(NRT_OR_TES_location[,c(2:4,1,6,5)],paste0("./13_TES_signal/","NRT_OR_TES_location",".bed"),row.names=F,col.names=F,sep="\t")
	write.table(RT_OR_TES_location[,c(2:4,1,6,5)],paste0("./13_TES_signal/","RT_OR_TES_location",".bed"),row.names=F,col.names=F,sep="\t")


OR_TES_location<- OR_TES_location[OR_TES_location$OR%in% OR_type$dotplot_feature,]


for(i in 1:nrow(OR_TES_location)){
	OR_TES_location_subset <- OR_TES_location[i,];
	#OR_TES_location_subset$V1<- 1;
	OR<- OR_TES_location_subset$OR
	write.table(OR_TES_location_subset[,c(2:4,1,5)],paste0("./01_OR_TES_location_bed/",OR,".bed"),row.names=F,col.names=F,sep="\t")
}
write.table(OR_TES_location[,c(2:5,1)],"./OR_TES_location.bed",row.names=F,col.names=F,sep="\t")

# back to shell 
sed -i 's/"//g' *.bed

ls *.bed > OR_promoter.list
for file1 in $(<OR_promoter.list)
do
  bedtools getfasta -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed $file1 -fo ../02_OR_TES_location_fa/$file1.fa -s -name
done




conda activate python27
ls *.fa > sample.list
for sample in $(<sample.list)
do
  jellyfish count -m 5 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer5_count.tsv
  jellyfish count -m 7 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer7_count.tsv  
  jellyfish count -m 9 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer9_count.tsv
  jellyfish count -m 11 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer11_count.tsv
  jellyfish count -m 15 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer15_count.tsv
  jellyfish count -m 19 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer19_count.tsv
  jellyfish count -m 25 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer25_count.tsv
  cat kmer*_count.tsv > $sample.tsv
  rm kmer*_count.tsv
done

rename .bed.fa.tsv .tsv *.tsv
perl -lne 'if ($ARGV=~/(.*).tsv/){print "$1\t$_"}' *.tsv >matrix.count

a=read.csv('/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/02_OR_TES_location_fa/matrix.count',header=F,sep="\t")
colnames(a)=c('OR','kmer','count')
library(reshape2)
counts=dcast(a,formula=kmer~OR)
counts[is.na(counts)]=0
write.table(counts, file="OR_TES_kmer.counts",sep="\t",quote=FALSE,row.names=FALSE)
rownames(counts)<- counts[,1]
counts=counts[,-1]

# make the tree 
library(lsa)
counts<- as.matrix(counts)
cosine_dist <- as.dist(cor(counts))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))

OR_type<- read.csv("./dotplot_feature_TES_type.csv")
RT_OR<- OR_type[OR_type$type=="RT",2]

library(ggtree)
pdf("./dotplot_feature_TES_tree-cosine.pdf",width=12,height=10)
tree <- groupOTU(data.tree, .node=RT_OR)
ggtree(tree,aes(color=group),layout="circular",) + geom_tiplab()+ geom_treescale() 
dev.off()

# DE kmer
library(DESeq2)
counts<- counts[,c(RT_OR,NRT_OR)]
Group <- factor(c(rep("RT",length(RT_OR)),rep("NRT",length(NRT_OR))))
colData <- data.frame(row.names=colnames(counts),Group)
colData
countData<-na.omit(counts)
dds<-DESeqDataSetFromMatrix(countData,colData,formula(~Group)) 
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)
res <- results(dds)
DE_KMER <- subset(res, pvalue < 0.05 )
count=t(scale(t(counts[rownames(DE_KMER),]),scale = T,center = T))
annotation_col = data.frame(colData$Group)
rownames(annotation_col) = factor(rownames(colData))

pdf("DE_KMER-heatmap.pdf",width = 14,height = 8)
pheatmap(count,cluster_cols =F,cluster_rows = T,
        # color = colorRampPalette(colors = c("navy","white","firebrick3")),
         #legend_breaks=seq(-3,3,1),
         annotation_col = annotation_col, 
         #annotation_colors =anncol,
         #annotation_row = annotation_row,
         show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols =T,cluster_rows = T,
        # color = colorRampPalette(colors = c("navy","white","firebrick3")),
         #legend_breaks=seq(-3,3,1),
         annotation_col = annotation_col, 
         #annotation_colors =anncol,
         #annotation_row = annotation_row,
         show_rownames=T,show_colnames=T)

dev.off()

bedtools getfasta -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed OR_TES_location.bed -fo OR_TES_location.fa -s -name
bedtools getfasta -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed NRT_OR_TES_location.bed -fo NRT_OR_TES_location.fa -s -name
bedtools getfasta -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed RT_OR_TES_location.bed -fo RT_OR_TES_location.fa -s -name
nohup meme OR_TES_location.fa -dna -oc ./03_MEME/MEME_all_OR_TES/ -nostatus -mod zoops -nmotifs 20 -minw 6 -maxw 60 &
nohup meme NRT_OR_TES_location.fa -neg ./03_MEME/RT_OR_TES_location.fa -dna -oc ./MEME_NRTvsRT_OR_TES/ -nostatus -mod zoops -nmotifs 20 -minw 6 -maxw 60 -objfun de &
nohup meme RT_OR_TES_location.fa -neg ./03_MEME/NRT_OR_TES_location.fa -dna -oc ./MEME_RTvsNRT_OR_TES/ -nostatus -mod zoops -nmotifs 20 -minw 6 -maxw 60 -objfun de &
nohup meme NRT_OR_TES_location.fa -dna -oc ./03_MEME/MEME_NRT_OR_TES/ -nostatus -mod zoops -nmotifs 20 -minw 6 -maxw 30 &
nohup meme RT_OR_TES_location.fa  -dna -oc ./03_MEME/MEME_RT_OR_TES/ -nostatus -mod zoops -nmotifs 20 -minw 6 -maxw 30 &

findMotifsGenome.pl OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa . -len 6,8,10,12,14,16,18,20,21,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50

nohup findMotifsGenome.pl RT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./HOMER_RT -len 6,8,10,12,14,16,18,20,21,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50 &
nohup findMotifsGenome.pl NRT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./HOMER_NRT -len 6,8,10,12,14,16,18,20,21,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50 &



