# Homer for NRT OR 

# nohup findMotifsGenome.pl RT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./HOMER_RT -len 5,6,8,10,12 -S 100 -size given -p 6 -oligo &
# nohup findMotifsGenome.pl NRT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./HOMER_NRT -len 5,6,8,10,12 -S 100 -size given -p 6 -oligo &
# nohup findMotifsGenome.pl OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./HOMER_all -len 5,6,8,10,12 -S 100 -size given -p 6 -oligo &
# 
# Optimize motifs ("-opt <motif file>")
# Instead of looking for novel de novo motifs, HOMER will instead try to optimize the motif supplied.  This is cool when trying to change the length of a motif, or find a very long version of a given motif.  For example, if you specify "-opt <file>" and "-len 50", it will try to expand the motif to 50bp and optimize it.
# 
# nohup findMotifsGenome.pl NRT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./HOMER_NRTOptimize -len 5,6,8,10,12 -S 100 -size given -p 6 -oligo -opt /data/R02/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/HOMER_NRT/homerResults/motif4.motif &
# 
# 
# 
# 
# # Finding Instance of Specific Motifs
# # get the oligoT motif 
# # Scanning sequence for motifs
# findMotifsGenome.pl ERalpha.peaks hg18 MotifOutputDirectory/ -find motif1.motif > outputfile.txt
# annotatePeaks.pl ERalpha.peaks hg18 -m motif1.motif >输出文件.txt
# 
# findMotifsGenome.pl NRT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa ./MotifOutputDirectory/ -find /data/R02/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/HOMER_NRT/homerResults/motif4.motif  > outputfile_motif_find.txt
# annotatePeaks.pl NRT_OR_TES_location.bed /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -m /data/R02/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/HOMER_NRT/homerResults/motif4.motif  > outputfile_motif.txt

# in R 

# 定义函数来提取上游序列
# for the strand +
max_length = 30;
extract_upstream_sequence_strand_p <- function(dna_sequence, pattern, max_length) {
  # 查找`ployT`尾巴的位置
  match <- gregexpr(pattern, dna_sequence)
  match <- as.integer(match[[1]])
  if (length(match) == 0) {
    return(NULL)  # 如果没有匹配到`ployT`尾巴，返回NULL
  }
  upstream_sequence <- c()
  if (length(match)==1) {
  # 计算`ployT`上游序列的起始位置
  start <- max(match - max_length, 1)
  # 提取上游序列
  upstream_sequence <- substring(dna_sequence, start, match - 1)
  return(upstream_sequence)
  }
  if (length(match)>1) {
    for(i in 1:length(match)){
      start <- max(match[i] - max_length, 1)
      upstream_sequence_subset <- substring(dna_sequence, start, match[i] - 1)
      upstream_sequence <- c(upstream_sequence_subset,upstream_sequence)
    }
    return(upstream_sequence)
  }
}

# for the strand -
extract_upstream_sequence_strand_n <- function(dna_sequence, pattern, max_length) {
  # 查找`ployT`尾巴的位置
  match <- gregexpr(pattern, dna_sequence)
  match <- as.integer(match[[1]])
  if (length(match) == 0) {
    return(NULL)  # 如果没有匹配到`ployT`尾巴，返回NULL
  }
  upstream_sequence <- c()
  pattern_length <- length(strsplit(pattern,split="")[[1]])
  sequence_length <- length(strsplit(dna_sequence,split="")[[1]])
  if (length(match)==1) {
  # 计算`ployT`上游序列的起始位置
  start <- match+pattern_length
  end <- min(start+max_length, sequence_length)
  # 提取上游序列
  upstream_sequence <- substring(dna_sequence, start, end)
  return(upstream_sequence)
  }
  if (length(match)>1) {
    for(i in 1:length(match)){
      start <- match+pattern_length
      end <- min(start+max_length, sequence_length)
      upstream_sequence_subset <- substring(dna_sequence, start, end)
      upstream_sequence <- c(upstream_sequence_subset,upstream_sequence)
    }
    return(upstream_sequence)
  }
}

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
strand_p_gene <- strand_p$gene_name
strand_n_gene <- strand_n$gene_name

library(Biostrings)
# for NRT_OR
NRT_OR_fasta<-readDNAStringSet("/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/NRT_OR_TES_location.fa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
upstream_sequence <- data.frame()
pattern <- "TTTTT"
max_length <- 50
for(i in 1:53){
  name<- strsplit(names(NRT_OR_fasta[i]),split="::")[[1]][1];
  sequence <- toString(NRT_OR_fasta[i])
  if(name %in% strand_p_gene){
      upstream_sequence_subset <- extract_upstream_sequence_strand_p(sequence, pattern,max_length);
      tmp<- data.frame(name,sequence=upstream_sequence_subset)
      upstream_sequence <- rbind(upstream_sequence,tmp)
  }
  if(name %in% strand_n_gene){
      upstream_sequence_subset <- extract_upstream_sequence_strand_n(sequence, pattern,max_length);
      tmp<- data.frame(name,sequence=upstream_sequence_subset)
      upstream_sequence <- rbind(upstream_sequence,tmp)
  };
}
upstream_sequence_last <- upstream_sequence[upstream_sequence$sequence!="",]
upstream_sequence_last <- upstream_sequence_last[!duplicated(upstream_sequence_last),]
# 将数据框转换为FASTA格式的对象
fasta_obj <- DNAStringSet(as.character(upstream_sequence_last$sequence))
names(fasta_obj)<- upstream_sequence_last$name
# 使用writeXStringSet函数保存为FASTA文件
writeXStringSet(fasta_obj, file = "/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/NRT_OR_TES_upstream50bp.fa")

# for RT_OR
RT_OR_fasta<-readDNAStringSet("/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RT_OR_TES_location.fa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
upstream_sequence <- data.frame()
pattern <- "TTTTT"
max_length <- 50
for(i in 1:14){
  name<- strsplit(names(RT_OR_fasta[i]),split="::")[[1]][1];
  sequence <- toString(RT_OR_fasta[i])
  if(name %in% strand_p_gene){
      upstream_sequence_subset <- extract_upstream_sequence_strand_p(sequence, pattern,max_length);
      tmp<- data.frame(name,sequence=upstream_sequence_subset)
      upstream_sequence <- rbind(upstream_sequence,tmp)
  }
  if(name %in% strand_n_gene){
      upstream_sequence_subset <- extract_upstream_sequence_strand_n(sequence, pattern,max_length);
      tmp<- data.frame(name,sequence=upstream_sequence_subset)
      upstream_sequence <- rbind(upstream_sequence,tmp)
  };
}
upstream_sequence_last <- upstream_sequence[upstream_sequence$sequence!="",]
upstream_sequence_last <- upstream_sequence_last[!duplicated(upstream_sequence_last),]
# 将数据框转换为FASTA格式的对象
fasta_obj <- DNAStringSet(as.character(upstream_sequence_last$sequence))
names(fasta_obj)<- upstream_sequence_last$name
# 使用writeXStringSet函数保存为FASTA文件
writeXStringSet(fasta_obj, file = "/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/RT_OR_TES_upstream50bp.fa")

# MEME denovo motif enrichment for NRT and RT and NRT vs RT 
nohup meme NRT_OR_TES_upstream30bp.fa \
-neg ./RT_OR_TES_upstream30bp.fa \
-dna -oc ./04_MEME_for_upstream30bp/MEME_NRTvsRT_OR_TES_up30/ \
-nostatus -mod zoops -nmotifs 20 -minw 5 -maxw 20 -objfun de &

nohup meme RT_OR_TES_upstream30bp.fa \
-neg ./NRT_OR_TES_upstream30bp.fa \
-dna -oc ./04_MEME_for_upstream30bp/MEME_RTvsNRT_OR_TES_up30/ \
-nostatus -mod zoops -nmotifs 20 -w 5 -minw 5 -maxw 20 -objfun de &

nohup meme NRT_OR_TES_upstream30bp.fa -dna -oc ./04_MEME_for_upstream30bp/MEME_NRT_OR_TES_up30/  -nostatus -mod zoops -nmotifs 20 -minw 5 -maxw 20 &
nohup meme RT_OR_TES_upstream30bp.fa  -dna -oc ./04_MEME_for_upstream30bp/MEME_RT_OR_TES_up30/ -nostatus -mod zoops -nmotifs 20 -minw 5 -maxw 20 &

# for 50bp

nohup meme NRT_OR_TES_upstream50bp.fa \
-neg ./RT_OR_TES_upstream50bp.fa \
-dna -oc ./05_MEME_for_upstream50bp/MEME_NRTvsRT_OR_TES_up50/ \
-nostatus -mod zoops -nmotifs 20 -minw 5 -maxw 20 -objfun de &

nohup meme RT_OR_TES_upstream50bp.fa \
-neg ./NRT_OR_TES_upstream50bp.fa \
-dna -oc ./05_MEME_for_upstream50bp/MEME_RTvsNRT_OR_TES_up50/ \
-nostatus -mod zoops -nmotifs 20 -w 5 -minw 5 -maxw 20 -objfun de &

nohup meme NRT_OR_TES_upstream50bp.fa -dna -oc ./05_MEME_for_upstream50bp/MEME_NRT_OR_TES_up50/  -nostatus -mod zoops -nmotifs 20 -minw 5 -maxw 20 &
nohup meme RT_OR_TES_upstream50bp.fa  -dna -oc ./05_MEME_for_upstream50bp/MEME_RT_OR_TES_up50/ -nostatus -mod zoops -nmotifs 20 -minw 5 -maxw 20 &





# RT vs NRT Intergenic region length
# for strand "+" OR gene 
strand_p <- strand_p[!duplicated(strand_p$gene_name),]
strand_p <- strand_p[order(strand_p$seqnames,strand_p$start,decreasing=F),]
strand_p_intersections<- data.frame()
for (i in 1:nrow(strand_p)){
  OR1 <- strand_p[i,]
  if(i!=nrow(strand_p)){
  OR2 <- strand_p[i+1,]
  if(OR1$seqnames==OR2$seqnames){
    intersections <- OR2$start-OR1$end
    tmp<- data.frame(OR=OR1$gene_name,intersections=intersections)
    strand_p_intersections<- rbind(strand_p_intersections,tmp)
  }
  }
}
# for strand "-" OR gene 
strand_n <- strand_n[!duplicated(strand_n$gene_name),]
strand_n <- strand_n[order(strand_n$seqnames,strand_n$start,decreasing=T),]
strand_n_intersections<- data.frame()
for (i in 1:nrow(strand_n)){
  OR1 <- strand_n[i,]
  if(i!=nrow(strand_n)){
  OR2 <- strand_n[i+1,]
  if(OR1$seqnames==OR2$seqnames){
    intersections <- -(OR2$start-OR1$end)
    tmp<- data.frame(OR=OR1$gene_name,intersections=intersections)
    strand_n_intersections<- rbind(strand_n_intersections,tmp)
  }
}
}

intersections_last<- rbind(strand_n_intersections,strand_p_intersections)

# add the RT info 
OR_type<- read.csv("../13_TES_signal/dotplot_feature_TES_type.csv")
intersections_last$type <- OR_type[match(intersections_last$OR,OR_type$dotplot_feature),]$type
intersections_last <- na.omit(intersections_last)

intersections_last$type<-factor(intersections_last$type,levels=c("NRT","RT"))
library(ggpubr)
library(cowplot)

pdf("./13_TES_signal/NRTvsRT_Intergenic_region_length_distribution.pdf",width=3,height=3)
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,2000)+
stat_compare_means()+theme(legend.position="none")+ylab("Intergenic region length")
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,20000)+
stat_compare_means()+theme(legend.position="none")+ylab("Intergenic region length")
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("Intergenic region length")

dev.off()


# for all OR gene 


library(Biostrings)
# for NRT_OR
OR_fasta<-readDNAStringSet("/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/NRT_OR_TES_location.fa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
upstream_sequence <- data.frame()
pattern <- "TTTTT"
max_length <- 50
for(i in 1:53){
  name<- strsplit(names(NRT_OR_fasta[i]),split="::")[[1]][1];
  sequence <- toString(NRT_OR_fasta[i])
  if(name %in% strand_p_gene){
      upstream_sequence_subset <- extract_upstream_sequence_strand_p(sequence, pattern,max_length);
      tmp<- data.frame(name,sequence=upstream_sequence_subset)
      upstream_sequence <- rbind(upstream_sequence,tmp)
  }
  if(name %in% strand_n_gene){
      upstream_sequence_subset <- extract_upstream_sequence_strand_n(sequence, pattern,max_length);
      tmp<- data.frame(name,sequence=upstream_sequence_subset)
      upstream_sequence <- rbind(upstream_sequence,tmp)
  };
}
upstream_sequence_last <- upstream_sequence[upstream_sequence$sequence!="",]
upstream_sequence_last <- upstream_sequence_last[!duplicated(upstream_sequence_last),]
# 将数据框转换为FASTA格式的对象
fasta_obj <- DNAStringSet(as.character(upstream_sequence_last$sequence))
names(fasta_obj)<- upstream_sequence_last$name
# 使用writeXStringSet函数保存为FASTA文件
writeXStringSet(fasta_obj, file = "/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/NRT_OR_TES_upstream50bp.fa")




