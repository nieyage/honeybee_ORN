# get the relationship between the sequence distance and kmer correalation 
library(ggplot2)
library(dplyr)
set.seed(1234)
# sequence distance
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
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
#Fig3A OR correlation  heatmap and sequence tree 
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
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
# corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature,]
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)

sequence_dist<- as.matrix(sdist)

# kmer correalation 
# get the promoter nucleotide sequence of dotplot OR 
ls *.bed > OR_promoter.list
for file1 in $(<OR_promoter.list)
do
  bedtools getfasta -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed $file1 -fo $file1.fa
done

conda activate python27
ls *.fa > sample.list
for sample in $(<sample.list)
do
print ${sample}
  jellyfish count -m 5 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer5_count.tsv
  jellyfish count -m 7 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer7_count.tsv  
  jellyfish count -m 9 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer9_count.tsv
  jellyfish count -m 11 -s 10G -t 12 -o mer_counts.jf $sample
  jellyfish dump  -c -t mer_counts.jf > kmer11_count.tsv
  cat kmer5_count.tsv kmer7_count.tsv kmer9_count.tsv kmer11_count.tsv > $sample.tsv
  rm kmer5_count.tsv kmer7_count.tsv kmer9_count.tsv kmer11_count.tsv
done

rename _promoter.fa.tsv .tsv *.tsv
perl -lne 'if ($ARGV=~/(.*).tsv/){print "$1\t$_"}' *.tsv >matrix.count

a=read.csv('/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/matrix.count',header=F,sep="\t")
colnames(a)=c('OR','kmer','count')
library(reshape2)
counts=dcast(a,formula=kmer~OR)
counts[is.na(counts)]=0
write.table(counts, file="/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/OR_kmer.counts",sep="\t",quote=FALSE,row.names=FALSE)
rownames(counts)<- counts[,1]
counts=counts[,-1]

# k=5,7,9,11
#data_cor <- counts

# k=5 
data_cor<- counts[which(str_length(rownames(counts))==5),]
#data_cor<- data_cor[which(rowSums(data_cor)>5),]
OR_kmer_cor <- cor(data_cor)

OR_kmer_cor <-OR_kmer_cor[which(rownames(OR_kmer_cor)%in%dotplot_feature),which(colnames(OR_kmer_cor)%in%dotplot_feature)]

OR_kmer_cor_long<- melt(OR_kmer_cor)
sequence_dist_long<- melt(sequence_dist)
sequence_dist_long$sample<- paste(sequence_dist_long$Var1,sequence_dist_long$Var2,sep="_")
OR_kmer_cor_long$sample<- paste(OR_kmer_cor_long$Var1,OR_kmer_cor_long$Var2,sep="_")
data <- merge(sequence_dist_long,OR_kmer_cor_long,by="sample")
data<- data[,c(1,4,7)]
colnames(data)<- c("OR_pair","sequence_dist","OR_kmer_cor")
data<- data[!duplicated(data),]
data<- data[data$sequence_dist!=0,]
data<- data[data$sequence_dist<300,]
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/kmer5_sequence_distm300_OR_kmer.pdf",width=6,height=4)
 ggplot(data, aes(sequence_dist,OR_kmer_cor))+ 
  geom_point()+ 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  theme_classic()+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
dev.off()

# file 300 dist 

# density plot 
data_kmer_freq<- rowSums(counts)
data<-as.data.frame(data_kmer_freq);
library(ggplot2)
pdf("./data_kmer_freq.pdf",width=10,height=4)
ggplot(data, aes(x=data_kmer_freq)) + xlab("gene_CV")+
              geom_density(alpha=.25) + theme_classic() 
dev.off()

data_cor<- counts
data_cor<- data_cor[which(rowSums(data_cor)>20),]
data_kmer_freq<- rowSums(data_cor)
data<-as.data.frame(data_kmer_freq);
library(ggplot2)
pdf("./data_kmer_freq.pdf",width=10,height=4)
ggplot(data, aes(x=data_kmer_freq)) + xlab("gene_CV")+
              geom_density(alpha=.25) + theme_classic() 
dev.off()

OR_kmer_cor <- cor(data_cor)
OR_kmer_cor <-OR_kmer_cor[which(rownames(OR_kmer_cor)%in%dotplot_feature),which(colnames(OR_kmer_cor)%in%dotplot_feature)]
OR_kmer_cor_long<- melt(OR_kmer_cor)
sequence_dist_long<- melt(sequence_dist)
sequence_dist_long$sample<- paste(sequence_dist_long$Var1,sequence_dist_long$Var2,sep="_")
OR_kmer_cor_long$sample<- paste(OR_kmer_cor_long$Var1,OR_kmer_cor_long$Var2,sep="_")
data <- merge(sequence_dist_long,OR_kmer_cor_long,by="sample")
data<- data[,c(1,4,7)]
colnames(data)<- c("OR_pair","sequence_dist","OR_kmer_cor")
data<- data[!duplicated(data),]
data<- data[data$sequence_dist!=0,]
data<- data[data$sequence_dist<300,]
library(ggpmisc)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/filter_20_kmer_sequence_dist300_OR_kmer.pdf",width=6,height=4)
 ggplot(data, aes(sequence_dist,OR_kmer_cor))+ 
  geom_point()+ 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  theme_classic()+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
dev.off()


