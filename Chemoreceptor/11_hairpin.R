RNAfold -j6 --noconv --auto-id NRT_OR_TES_upstream50bp.fa > NRT_OR_TES_upstream50bp.out
RNAfold -j6 --noconv --auto-id RT_OR_TES_upstream50bp.fa > RT_OR_TES_upstream50bp.out


sort -t ">" -k 2,2 -o sorted_sequences.fasta NRT_OR_TES_upstream50bp.fa

grep ">" RT_OR_TES_upstream50bp.fa>RT_ID
grep -E -o "\(.+?\)" RT_OR_TES_upstream50bp.out >RTresult
paste RT_ID RTresult |sort -k1,1 


R

data<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig7/dotplot_feature_TES_type_hairpin2.csv")
# The number of hairpin in OR TES region 
my_comparisons <- list( c("NRT", "RT") )
library(ggpubr)
pdf("./00_Figure/Fig7/Fig7M_hairpin_number_NRT_vs_RT.pdf",width=4,height=4)
ggboxplot(data, x="type", y="Hairpin_number", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("# of hairpin")
dev.off()

RNAfold --noconv -p < LOC100578045.fa > LOC100578045.res

/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl LOC100578045_ss.ps LOC100578045_dp.ps > LOC100578045_rss.ps
magick convert -density 300 LOC100578045_rss.ps LOC100578045_rss.pdf # 生成二级结构pdf图




OR_fasta<-readDNAStringSet("/data/R02/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/02_OR_TES_location_fa/all_OR_TES.fa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
upstream_sequence <- data.frame()
pattern <- "TTTTT"
max_length <- 50
for(i in 1:177){
  name<- strsplit(names(OR_fasta[i]),split="::")[[1]][1];
  sequence <- toString(OR_fasta[i])
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
writeXStringSet(fasta_obj, file = "/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/OR_TES_upstream50bp.fa")


RNAfold -j6 --noconv --auto-id OR_TES_upstream50bp.fa > OR_TES_upstream50bp.out

grep ">" OR_TES_upstream50bp.fa>OR_ID
grep -E -o "\(.+?\)" OR_TES_upstream50bp.out >all_RNAfold_result
paste OR_ID all_RNAfold_result |sort -k1,1 






