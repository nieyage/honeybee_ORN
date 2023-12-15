class calculate_genome:
    def __init__(self,gtf,fasta):
        self.gtf = gtf
        self.fasta = fasta
        
    def chromoLen(self):
        chrLen = {}
        for rec in SeqIO.parse(self.fasta,'fasta'):
            chrLen[rec.id] = len(rec.seq)
            
        return chrLen
    
    def gcContent(self,gc_window):
        self.gc_window = gc_window
        
        gc_content = {'chr_id':[],
                     'bin_start':[],
                     'gc':[]}
        chr_len = self.chromoLen()
        
        #print(chr_len)
        for rec in SeqIO.parse(self.fasta,'fasta'):
            for i in range(0,chr_len[rec.id],self.gc_window):
                #print(rec.id)
                gc_content['chr_id'].append(rec.id)
                gc_content['bin_start'].append(i)
                gc_content['gc'].append(round(GC(rec.seq[i:i+self.gc_window]),2))
        return pd.DataFrame(gc_content)
    
    def geneDensity(self,gene_window):
        self.gene_window = gene_window
        final_df = []
        df = pd.read_table(self.gtf,header=None,comment="#",sep="\t",
                           usecols=[0,2,3,4],
                           names="Chromosome Feature Start End".split())
        #df.columns = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()
        df = df[df.Feature=="transcript"]
        chrLen = self.chromoLen()
        for chr_id in chrLen.keys():
            print(chr_id)
            df1 = df[df.Chromosome==chr_id]
            gene_start = [int(a) for a in df1.Start]
            gene_start.insert(0,0)
            gene_start.append(round(chrLen[chr_id]/self.gene_window)*self.gene_window+self.gene_window)
            #print(gene_start)
            bin_start = [int(a.left) for a in pd.cut(gene_start,bins=round(chrLen[chr_id]/self.gene_window)+1).value_counts().index]
            bin_start[0] = 0
            gene_count = list(pd.cut(gene_start,bins=round(chrLen[chr_id]/self.gene_window)+1).value_counts().values)
            #print(len(bin_start))
            #print(len(gene_count))
            #print("OK")
            final_df.append(pd.DataFrame({'chr_id':chr_id,'bin_start':bin_start,'gene_count':gene_count}))
            
        return pd.concat(final_df)
        


# get OR gene gtf in R 
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
library(reshape2)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']
gene.coords <- gene.coords[gene.coords$gene_name %in% OR_gene]
gene.coords <- gene.coords[gene.coords$type == "transcript"]
gene.coords <- as.data.frame(gene.coords)
write.table(gene.coords, file = "OR_gene_filtered.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# 导入必要的库
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC

# 创建 calculate_genome 实例并传入 GTF 和 FASTA 文件的路径
honeybee = calculate_genome('/md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/OR_gene_filtered.gtf', '/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa')

# 计算染色体长度
chromosome_lengths = honeybee.chromoLen()
print(chromosome_lengths)

# 计算 GC 内容
gc_content_df = honeybee.gcContent(gc_window=10)
print(gc_content_df)

# 计算基因密度
gene_density_df = honeybee.geneDensity(gene_window=5000)
print(gene_density_df)

df = pd.DataFrame(chromosome_lengths.items(), columns=['Chromosome', 'Length'])

# Save the DataFrame to a CSV file
df.to_csv('honeybee_chromoLen.csv', index=False)

# 保存 GC 内容数据到 CSV 文件
gc_content_df.to_csv('honeybee_gc_content.csv', index=False)

# 保存基因密度数据到 CSV 文件
gene_density_df.to_csv('OR_gene_density.csv', index=False)



OR_number<- as.data.frame(table((as.data.frame(gene.coords)[,1])))
OR_number<- OR_number[order(OR_number$Freq,decreasing=TRUE),]
OR_chrom<- as.character(OR_number[OR_number$Freq>5,]$Var1)
df<-read.csv("honeybee_chromoLen.csv")
gc<-read.csv("honeybee_gc_content.csv")
genedensity<-read.csv("OR_gene_density.csv")
colnames(df)<- c("chr_id","chr_len")
df<- df[df$chr_id%in%OR_chrom,]
gc<- gc[gc$chr_id%in%OR_chrom,]
genedensity<- genedensity[genedensity$chr_id%in%OR_chrom,]

gc_content_df<- data.frame(seq_ID=gc$chr_id,seq_start=gc$bin_start,  seq_end=gc$bin_start+999, Abundance=gc$gc/sum(gc$gc))
write.csv(gc_content_df,"gc_content_df_density.csv")


OR_gene_df<- data.frame(seq_ID=gene.coords$seqnames,seq_start=gene.coords$start,  seq_end=gene.coords$end, Abundance=1)

write.csv(OR_gene_df,"OR_gene_df.csv")

gene.coords[gene.coords$gene_name%in%c("Or25","Or26","Or27"),]

  gene.coords[gene.coords$gene_name%in%c("Or9","Or10","Or11"),],]
   seqnames    start      end width strand     source       type score phase
47   Group2 10000277 10002223  1947      + BestRefSeq transcript    NA    NA
48   Group2 10004087 10005752  1666      + BestRefSeq transcript    NA    NA
49   Group2 10008722 10010579  1858      + BestRefSeq transcript    NA    NA
    transcript_id gene_id gene_name   gene_biotype
47 NM_001242960.1     Or9       Or9 protein_coding
48 NM_001242961.1    Or10      Or10 protein_coding
49 NM_001242962.1    Or11      Or11 protein_coding

  gene.coords[gene.coords$gene_name%in%c("Or25","Or26","Or27"),],]
   seqnames    start      end width strand     source       type score phase
55   Group2 10046322 10048685  2364      +     Gnomon transcript    NA    NA
56   Group2 10049035 10051106  2072      + BestRefSeq transcript    NA    NA
57   Group2 10051858 10053414  1557      + BestRefSeq transcript    NA    NA
    transcript_id gene_id gene_name   gene_biotype
55 XM_016917965.2    Or25      Or25 protein_coding
56 NM_001242968.1    Or26      Or26 protein_coding
57 NM_001242969.1    Or27      Or27 protein_coding

# 
https://www.omicstudio.cn/tool/50






data<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")
data$gene_exp_pattern<- dotplot_data[match(data$OR_gene,dotplot_data$features.plot),]$gene_type
write.csv(data,"./OR_gene_df_gene_exp_pattern.csv")
# UCSC show region chrLG2:9976746-10160000

sed 's/\(chr[0-9]\+\)/"\1"/g' OR_name_map_output_ucsc.gtf > OR_name_map_output_ucscmodified.gtf


# For OR9-11 chrLG2:10000000-10013000

# For Or25-27 chrLG2:10046200-10055000


library(ggplot2)
# 读取GTF文件，假设GTF文件名为"yourfile.gtf"
gtf_data <- read.csv("/Users/fraya/Documents/project/honeybee/chemoreceptor/OR_gene_df_gene_exp_pattern_Group2.csv", header = FALSE)

# 指定要绘制的区域范围，这里假设绘制chr1的1-10000范围内的基因
chr <- "Group2"
start_pos <- 9976746
end_pos <- 10160000
head(gtf_data)
# 筛选符合区域范围的基因
genes <- subset(gtf_data, V2 == chr & V3 >= start_pos & V4 <= end_pos)

gg <- ggplot(data = genes, aes(x = (V4 + V3) / 2, y = V5, shape = V5,color=V6)) +
  geom_jitter(height = 100, width = 0, size = 5) +
  #geom_point(size = 5) +
  scale_shape_manual(values = c(24, 25)) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "") +
  theme(axis.text = element_blank(), axis.title = element_text(size = 12))
gg <- gg + geom_text(aes(label = V1), vjust = -0.7, size = 4,  nudge_y = 100)

# 打印图形
print(gg)

genes <- subset(gtf_data, V1 == chr & V4 >= start_pos & V5 <= end_pos & V3 == "gene")
library(ggplot2)
library(gggenes)
BiocManager::install("gggenes")

ggplot(genes, aes(xmin = V3, xmax = V4, 
                          y = V2, fill = V6, label = V1, forward = 1)) +
  geom_gene_arrow() +
  facet_wrap(~ V2, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  geom_gene_label(align = "centre",min.size = 2,reflow=T)



# gene type for dotplot OR 
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
dotplot_data<-read.csv("../05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
coexp_OR <- unique(dotplot_data[dotplot_data$id%in% multiOR_cluster,]$features.plot)
data<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/dotplot_feature_TES_type.csv")
RT_OR<- data[data$type=="RT",2]
dotplot_data$gene_type <- ifelse(dotplot_data$features.plot%in% coexp_OR,"coexp(NRT)","single_OR")
dotplot_data[which(dotplot_data$features.plot%in% RT_OR),]$gene_type="coexp(RT)"

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

# Step1: Intergenic region length

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

intersections_last$type <- dotplot_data[match(intersections_last$OR,dotplot_data$features.plot),]$gene_type
intersections_last <- na.omit(intersections_last)

intersections_last$type<-factor(intersections_last$type,levels=c("single_OR","coexp(NRT)","coexp(RT)"))
library(ggpubr)
library(cowplot)
my_comparisons <- list( c("single_OR", "coexp(NRT)"), c("single_OR", "coexp(RT)"), c("coexp(NRT)", "coexp(RT)") )

pdf("./singlevsNRTvsRT_Intergenic_region_length_distribution.pdf",width=6,height=3)
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,2000)+
stat_compare_means(comparisons = my_comparisons, label.y = c(18, 22, 26))+theme(legend.position="none")+ylab("Intergenic region length")
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,20000)+
stat_compare_means(comparisons = my_comparisons, label.y = c(18, 22, 26))+theme(legend.position="none")+ylab("Intergenic region length")
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("Intergenic region length")
dev.off()

# Step2: dN/dS(compare with single exp OR) coexp (NRT vs RT) 
library(ggpubr)
library(cowplot)
kaks_data<- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/result_dir_muscle2/OR_pair_kaks.out")
colnames(kaks_data)<- c("OR_pair","method","KA","Ks","Ka/Ks","p_value","length","S-sites","N-sites","Folo-sites","Substitutions","S-sub","N-sub","Folo-S-sub","Folo-N-sub","Divergence-Time","Substitution-Rate-Ratio","GC","ML_score","AICc","Akaike-Weight","Model")
dNdS_data <- kaks_data[,c(1,5)]
library(tidyr)
df <- kaks_data %>%
  separate(OR_pair, into = c("OR1", "OR2"), sep = "-")

coexp_OR_RT<- dotplot_data[dotplot_data$gene_type=="coexp(RT)",]$features.plot
coexp_OR_NRT<- dotplot_data[dotplot_data$gene_type=="coexp(NRT)",]$features.plot
single_OR <- dotplot_data[dotplot_data$gene_type=="single_OR",]$features.plot
coexp_OR_RT_vs_single_OR <- df[df$OR1%in% coexp_OR_RT,]
coexp_OR_RT_vs_single_OR <- coexp_OR_RT_vs_single_OR[coexp_OR_RT_vs_single_OR$OR2%in% single_OR,]
coexp_OR_NRT_vs_single_OR <- df[df$OR1%in% coexp_OR_NRT,]
coexp_OR_NRT_vs_single_OR <- coexp_OR_NRT_vs_single_OR[coexp_OR_NRT_vs_single_OR$OR2%in% single_OR,]

data<- data.frame(type=c(rep("coexp(NRT)",length(coexp_OR_NRT_vs_single_OR[,6])),
  rep("coexp(RT)",length(coexp_OR_RT_vs_single_OR[,6]))),
dNdS=c(coexp_OR_NRT_vs_single_OR[,6],coexp_OR_RT_vs_single_OR[,6])
)
my_comparisons <- list( c("coexp(NRT)", "coexp(RT)") )

pdf("./NRTvsRT_dNdS_compare_singleOR.pdf",width=3,height=3)
ggboxplot(data, x="type", y="dNdS", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,1)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("dN/dS")
dev.off()

# Step3: ploy T number and hairpin number
a=read.csv('/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/02_OR_TES_location_fa/matrix.count',header=F,sep="\t")
colnames(a)=c('OR','kmer','count')
library(reshape2)
counts=dcast(a,formula=kmer~OR)
counts[is.na(counts)]=0
write.table(counts, file="OR_TES_kmer.counts",sep="\t",quote=FALSE,row.names=FALSE)
rownames(counts)<- counts[,1]
counts=counts[,-1]

polyT_sequences <- character(0)

# 生成包含5到15个T的polyT字符串并存储在字符向量中
for (length_T in 5:15) {
  polyT_sequence <- paste(rep("T", length_T), collapse = "")
  polyT_sequences <- c(polyT_sequences, polyT_sequence)
}
polyT_data<- data.frame(OR=colnames(counts),polyT_count=colSums(counts[which(rownames(counts)%in%polyT_sequences),]))
polyT_data$gene_type<- dotplot_data[match(polyT_data$OR,dotplot_data$features.plot),]$gene_type
my_comparisons <- list( c("single_OR", "coexp(NRT)"), c("single_OR", "coexp(RT)"), c("coexp(NRT)", "coexp(RT)") )

pdf("./NRTvsRT_polyT_count.pdf",width=4,height=4)
ggboxplot(polyT_data, x="gene_type", y="polyT_count", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("# of polyT")
dev.off()








