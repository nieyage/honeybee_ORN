
# Step1: parpare the OR gene cds and pep file (honeybee and ant)

honeybee_cds=/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/all_receptor_gene_cds.fa
honeybee_pep=/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/all_receptor_gene_pep.aa

# other bee sequence 

awk -F, '{if(NR>1) print ">"$1"\n"$2}' ant_OR_IR_GR_cds.csv > ant_OR_IR_GR_cds.fasta
awk -F, '{if(NR>1) print ">"$1"\n"$2}' ant_OR_IR_GR_pep.csv > ant_OR_IR_GR_pep.fasta

H_saltator_cds=/md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/ant_OR_IR_GR_cds.fasta
H_saltator_pep=/md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/ant_OR_IR_GR_pep.fasta


# Step2: homologous gene screening
cd /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/02_homologous_gene

###先使用makeblastdb建库，如果是老版本的blast的话，则使用formatdb
# 双向最优配对序列的序列名
makeblastdb -in $H_saltator_pep -dbtype prot
###使用blast进行alignment，得到一个表格，输出格式为m6
blastp -query $honeybee_pep -db $H_saltator_pep -evalue 1e-5 -max_target_seqs 1 -num_threads 2 -out honeybee_H_saltator_blastp_out.m6 -outfmt 6
###而后使用shell指令进行简单的排序跟取值，如果存在多个物种比对，则需要写更为复杂的脚本去取双向最优序列
cut -f1-2 honeybee_H_saltator_blastp_out.m6 | sort | uniq > honeybee_H_saltator.homolog

# Step3: remove the gap by ParaAT;
###此处用到的cds跟pep为这两个物种直接cat起来
###ParaAT需要全局调用从能正常使用，同时要确保clustalw已经安装，不然会报错，可以直接用conda装
###-p是指定线程数的文件，直接生成一个文件里面填个数字即可
cat $honeybee_cds $H_saltator_cds > /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/honeybee_H_saltator_cds.fa
cat $honeybee_pep $H_saltator_pep > /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/honeybee_H_saltator_pep.aa

ParaAT.pl -h ./02_homologous_gene/honeybee_H_saltator.homolog \
-n /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/honeybee_H_saltator_cds.fa \
-a /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/honeybee_H_saltator_pep.aa \
-m clustalw2 -p proc -f paml -o ./03_honeybee_H_saltator_paml

# Step4: CODEML in PAML for dN/dS
# install the PAML
wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
tar xf paml4.9j.tgz
rm bin/*.exe 
cd src 
make -f Makefile 
# rm *.o 
mv baseml basemlg codeml pamp evolver yn00 chi2 ../bin 

cat 03_honeybee_H_saltator_paml/*.paml > ./05_PAML/honeybee_H_saltator.codon

cd /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/05_PAML

ls -l |grep "^-"|wc -l
219
# tree file 
vi test.trees
  2  1

(1,2);
# config file 
cp /md01/nieyg/software/paml4.9j/codonml.ctl .

/md01/nieyg/software/paml4.9j/bin/codeml codonml.ctl


# Step5: calculate  the dNdS by kaks_calculate
conda activate kaks_calculate

ParaAT.pl -h ./02_homologous_gene/honeybee_H_saltator.homolog \
-n /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/honeybee_H_saltator_cds.fa \
-a /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/01_input_file/honeybee_H_saltator_pep.aa \
-m muscle -p proc -f axt -g -k -o ./04_honeybee_H_saltator_kaks_calculate

cat *kaks > merge.kaks
awk '!(NR%2)' merge.kaks > honeybee_H_saltator_kaks_calculate.out


# Step6: dN/dS(compare with single exp OR and coexp (NRT and RT) 

library(ggpubr)
library(cowplot)
kaks_data<- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/04_honeybee_H_saltator_kaks_calculate/honeybee_H_saltator_kaks_calculate.out")
colnames(kaks_data)<- c("OR_pair","method","KA","Ks","Ka/Ks","p_value","length","S-sites","N-sites","Folo-sites","Substitutions","S-sub","N-sub","Folo-S-sub","Folo-N-sub","Divergence-Time","Substitution-Rate-Ratio","GC","ML_score","AICc","Akaike-Weight","Model")
dNdS_data <- kaks_data[,c(1,5)]
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
coexp_OR <- unique(dotplot_data[dotplot_data$id%in% multiOR_cluster,]$features.plot)
data<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/dotplot_feature_TES_type.csv")
RT_OR<- data[data$type=="RT",2]
dotplot_data$gene_type <- ifelse(dotplot_data$features.plot%in% coexp_OR,"coexp(NRT)","single_OR")
dotplot_data[which(dotplot_data$features.plot%in% RT_OR),]$gene_type="coexp(RT)"

write.csv(dotplot_data,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_gene_exptype.csv")
library(tidyr)
df <- kaks_data %>%
  separate(OR_pair, into = c("OR1", "OR2"), sep = "-")

coexp_OR_RT<- dotplot_data[dotplot_data$gene_type=="coexp(RT)",]$features.plot
coexp_OR_NRT<- dotplot_data[dotplot_data$gene_type=="coexp(NRT)",]$features.plot
single_OR <- dotplot_data[dotplot_data$gene_type=="single_OR",]$features.plot

single_OR<- df[df$OR1%in% single_OR,]
coexp_OR_RT<- df[df$OR1%in% coexp_OR_RT,]
coexp_OR_NRT<- df[df$OR1%in% coexp_OR_NRT,]

data<- data.frame(type=c(rep("coexp_NRT",length(coexp_OR_NRT[,6])),
  rep("coexp_RT",length(coexp_OR_RT[,6])),
  rep("single_OR",length(single_OR[,6]))),
dNdS=c(coexp_OR_NRT[,6],coexp_OR_RT[,6],single_OR[,6]))
my_comparisons <- list( c("coexp_NRT", "coexp_RT"),c("coexp_NRT", "single_OR"),c("coexp_RT", "single_OR") )

pdf("./00_Figure/Fig7/dNdS_compare_split_RTandNRT_ant.pdf",width=5,height=5)
ggboxplot(data, x="type", y="dNdS", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,1)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("dN/dS")
dev.off()

data<- data.frame(type=c(rep("coexp_OR",length(coexp_OR_NRT[,6])+length(coexp_OR_RT[,6])),
  rep("single_OR",length(single_OR[,6]))),
dNdS=c(coexp_OR_NRT[,6],coexp_OR_RT[,6],single_OR[,6]))
my_comparisons <- list( c("coexp_OR", "single_OR") )

pdf("./00_Figure/Fig7/dNdS_compare_single_vs_coexp_ant.pdf",width=4,height=5)
ggboxplot(data, x="type", y="dNdS", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,1)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("dN/dS")
dev.off()

#!!!!!!! the results of GBE paper

data <- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/GBE_bee_dNdS_results.csv")

data <- na.omit(data)

OR_gene<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")
OR=unique(dotplot_data$features.plot)
OR_name<- OR_gene[match(OR,OR_gene$OR_gene),4]
OR_name[39]<- "AmGr9"
OR_name[61]<- "AmOr105"
#(base) grep "LOC100576522" OR_gene_mapping.txt | sort -k3,3nr -t$'\t'|less
#(base) grep "LOC107966050" OR_gene_mapping.txt | sort -k3,3nr -t$'\t'|less

OR_name[27]<- "AmOr162P" # pident 72.805
OR_name[19]<- "AmOr119" # pident 44.444

OR_type<- dotplot_data[match(OR,dotplot_data$features.plot),]$gene_type
dNdS_data<- data.frame(OR,OR_name,OR_type)

dNdS_data$OR_name <- gsub("[A-Za-z]*$", "", dNdS_data$OR_name)
data$Branch<- gsub("[A-Za-z]*$", "", data$Branch)
dNdS_data$dNdS<- data[match(dNdS_data$OR_name,data$Branch),]$Mean_dNdS

dNdS_data$exp_pattern<- dNdS_data$OR_type
dNdS_data$exp_pattern<- gsub("coexp\\(NRT\\)","coexp_NRT",dNdS_data$exp_pattern)
dNdS_data$exp_pattern<- gsub("coexp\\(RT\\)","coexp_RT",dNdS_data$exp_pattern)
my_comparisons <- list( c("coexp_NRT", "coexp_RT"),c("coexp_NRT", "single_OR"),c("coexp_RT", "single_OR") )

pdf("./00_Figure/Fig7/dNdS_compare_split_RTandNRT_bee.pdf",width=5,height=5)
ggboxplot(dNdS_data, x="exp_pattern", y="dNdS", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,1)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("dN/dS")
dev.off()

dNdS_data$exp_pattern<- dNdS_data$OR_type
dNdS_data$exp_pattern<- gsub("coexp\\(NRT\\)","coexp_OR",dNdS_data$exp_pattern)
dNdS_data$exp_pattern<- gsub("coexp\\(RT\\)","coexp_OR",dNdS_data$exp_pattern)

my_comparisons <- list( c("coexp_OR", "single_OR") )

pdf("./00_Figure/Fig7/dNdS_compare_single_vs_coexp_bee.pdf",width=4,height=5)
ggboxplot(dNdS_data, x="exp_pattern", y="dNdS", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,1)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("dN/dS")
dev.off()




